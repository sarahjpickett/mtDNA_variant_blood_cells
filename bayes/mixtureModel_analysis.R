###
# READ READ READ READ READ
#
# if JAGS not installed on computer: can install at following link
# https://mcmc-jags.sourceforge.io
# then can install package in R with the following command
# install.packages("rjags")
# If interested this is a blog about getting started with JAGS, rjags and bayesian modelling
# https://www.r-bloggers.com/2012/04/getting-started-with-jags-rjags-and-bayesian-modelling/
#
## - Install/update packages with following commands
# install.packages("rjags")
# install.packages("parallel")
#
###

library("rjags")
library("parallel")

dir.create("./OUTPUT", showWarnings = FALSE)

# load data from Sarah's repo
# The working directory should be Sarah's repo
sc_data = read.csv("../../output/all_cleaned_sc_data.csv", header=TRUE)

# sets data to be between 0
sc_data[,"HET"] = sc_data[,"HET"]/100
sc_data[,"cell_type"] = sub(" ", "", sc_data[,"cell_type"])
sc_data[,"cell_type"] = sub("/", "_", sc_data[,"cell_type"])

# identifying cell types and patients
cell_types = unique(sc_data$cell_type)
pts = unique(sc_data$Patient)

# Model that can be read by JAGS
modelstring = "
model {
  for(i in 1:N){
    # calculate log density for both components
    dens[i,1] = dnorm(Y[i], mu, tau) / ( phi((1-mu)*sqrt(tau))-phi(-mu*sqrt(tau)) ) # truncated normal [0,1] density
    dens[i,2] = 1 # uniform density

    # The ones trick
    prob[i] = (pi_infl*dens[i,1] + (1-pi_infl)*dens[i,2]) / Const
    ones[i] ~ dbern( prob[i] )
  }
  # priors
  pi ~ dbeta(alpha_pi, beta_pi)
  pi0 ~ dbeta(alpha_pi0, beta_pi0)
  z ~ dbern(pi0)
  pi_infl = ifelse(z==1, pi, 0)

  mu ~ dunif(lower_mu,upper_mu)
  sd ~ dgamma(alpha_sd, beta_sd)
  tau = 1/sd^2

  # prediction
  z_pred ~ dbern(pi_infl)
  comp_pred = 2 - z_pred
  pred[1] ~ dnorm(mu, tau) T(0,1)
  pred[2] ~ dunif(0, 1)
  Ypred = pred[comp_pred]
}
"

## number of samples from the posterior distribution
Npost = 5000
## level of thinning
# this removes correlation between successive samples from the posterior
# so we have independent draws from the posterior beliefs on parameters
MCMCthin = 100
## burn-in period 
# the number of samples allowed for the posterior MCMC chain to reach its th
# posterior distribution. If inference output looks more like a stock price
# rather than a caterpillar, increase this first.
MCMCBurnin = 5000

## function to perform inference and return prior and posterior values
inference = function(cellpat_data){
  if( nrow(cellpat_data)>0 ){
    type = unique(cellpat_data$cell_type)
    if( type == "8+ EM/TEMRA") type = "8+ EM-TEMRA"
    fp_post = paste0("./OUTPUT/", type,"_POST.csv")
    fp_prior = paste0("./OUTPUT/", type,"_PRIOR.csv")
    if(!file.exists(fp_post) | !file.exists(fp_prior)){ 
    # check that the output files do not already exist
      
      # define the parameter values to be passed to the model
      data = list(Y=cellpat_data$HET, ones=rep(1,nrow(cellpat_data)),
                  alpha_pi=1, beta_pi=1,
                  alpha_pi0=1, beta_pi0=1,
                  lower_mu=0, upper_mu=0.2,
                  alpha_sd=1, beta_sd=5,
                  N=nrow(cellpat_data), Const=1e2)
      # the Const parameter is part of the ones trick and ensures prob[i] is less than 1
      # its value value does not affect the output. If inference fails try increasing

      # data to sample from the priors
      pred = data
      pred$Y = NULL
      pred$N = 0
      pred$ones = NULL

      # where inference is performed, using rjags package
      model = jags.model(textConnection(modelstring), data=data, n.chains=1)
      update(model, n.iter=MCMCBurnin)
      output = coda.samples(model=model, variable.names=c("mu","sd","pi","pi0","pi_infl","Ypred","pred"),
                                  n.iter=Npost*MCMCthin, thin=MCMCthin)
      # samples from the priors
      model_prior = jags.model(textConnection(modelstring), data=pred, n.chains=1)
      output_prior = coda.samples(model=model_prior, variable.names=c("mu","sd","pi","pi0","pi_infl","Ypred","pred"),
                                        n.iter=Npost, thin=1)

      # saves the posterior samples in a data frame
      post_df = as.data.frame(output[[1]])
      prior_df = as.data.frame(output_prior[[1]])
      
      ## For each posterior and prior draw we calculate the point at which 
      # the two densities cross and the percentage of blood cells below this 
      # threshold
      threshold_post = double(Npost)
      threshold_prior = double(Npost)
      percent_below_post = double(Npost)
      percent_below_prior = double(Npost)
      for(i in 1:Npost){
        mu_post = post_df[i,"mu"]
        sd_post = post_df[i,"sd"]
        h_post =  1 - post_df[i,"pi_infl"]

        mu_prior = prior_df[i,"mu"]
        sd_prior = prior_df[i,"sd"]
        h_prior =  1 - prior_df[i,"pi_infl"]
        b=1
        a=0

        trunc_const_post = pnorm(b,mu_post,sd_post)-pnorm(a,mu_post,sd_post)
        trunc_const_prior = pnorm(b,mu_prior,sd_prior)-pnorm(a,mu_prior,sd_prior)

        thr_post = sqrt(-2*sd_post^2*log( sqrt(2*pi*sd_post^2)*h_post/(1-h_post)*trunc_const_post ))+mu_post
        thr_prior = sqrt(-2*sd_prior^2*log(sqrt(2*pi*sd_prior^2)*h_prior/(1-h_prior)*trunc_const_prior ))+mu_prior

        threshold_post[i] = ifelse(thr_post>0 & thr_post<1, thr_post, NA)
        threshold_prior[i] = ifelse(thr_prior>0 & thr_prior<1, thr_prior, NA)

        percent_below_post[i] = sum(cellpat_data[["HET"]]<threshold_post[i])/nrow(cellpat_data)
        percent_below_prior[i] = sum(cellpat_data[["HET"]]<threshold_prior[i])/nrow(cellpat_data)
      }
      # add the newly calculated variables to the data frames  
      post_df[["threshold"]] = threshold_post
      post_df[["percent_below"]]= percent_below_post
      prior_df[["threshold"]] = threshold_prior
      prior_df[["percent_below"]] = percent_below_prior
      
      # return prior and posterior data frames in list
      return(list("post"=post_df, "prior"=prior_df))
    }
  }
}

# create a list of data to be passed to the inference function
# for each cell and patient combination
data_list = list()
for(pat in pts){
  for(type in cell_types){
    data_list[[paste0(type,"_",pat)]] = sc_data[sc_data$cell_type==type & sc_data$Patient==pat, ]
  }
}

#data_list[["Bnaive_R013"]] = sc_data[sc_data$cell_type=="Bnaive" & sc_data$Patient=="R013",]

## perform inference across cores
# default to a one less core than is present on your machine 
# can change this by altering the ncores parameter
ncores = detectCores() - 1 
cl = makeCluster(ncores)
{
  clusterExport(cl, c("inference", "Npost", "MCMCthin", "MCMCBurnin", "modelstring"))
  clusterEvalQ(cl, {
    library("rjags")
  })
  output = parLapply(cl, data_list, inference)
}
stopCluster(cl)

# save the output as a .csv file
for(cellpat in names(output)){
  if(!is.null(output[[cellpat]])){
    fp_post = paste0("./OUTPUT/",cellpat,"_POST.csv")
    fp_prior = paste0("./OUTPUT/",cellpat,"_PRIOR.csv")
    write.csv(output[[paste(cellpat)]][["post"]], file=fp_post)
    write.csv(output[[paste(cellpat)]][["prior"]], file=fp_prior)
  }
}







