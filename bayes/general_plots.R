library("HDInterval")

# creates folders to save output in 
dir.create("./PDF", showWarnings = FALSE)
dir.create("./PDF/cellpat_bool", showWarnings = FALSE)

myBlack = function(alpha) rgb(0,0,0, alpha)
myGrey = function(alpha) rgb(66/255,66/255,66/255, alpha)
myLightBlue = function(alpha) rgb(0,191/255,1, alpha)
myBlue = function(alpha) rgb(0,0,128/255, alpha)
myCerulean = function(alpha) rgb(42/255, 82/255, 190/255, alpha)
myRed = function(alpha) rgb(1,0,0, alpha)
myPurple = function(alpha) rgb(148/255,0/255,211/255, alpha)
myGreen = function(alpha) rgb(0,100/255,0, alpha)
myPink = function(alpha) rgb(1,105/255,180/255, alpha)
myYellow = function(alpha) rgb(1,215/255,0, alpha)

myCols = function(alpha, labels){
  col_vec = c(myLightBlue(alpha), myBlue(alpha), myRed(alpha), myPink(alpha), myGreen(alpha), myYellow(alpha))
  names(col_vec) = labels
  return(col_vec)
}

rtnorm <- function(n, mu, sigma, low, high) {
  # find quantiles that correspond the the given low and high levels.
  p_low <- pnorm(low, mu, sigma)
  p_high <- pnorm(high, mu, sigma)
  
  # draw quantiles uniformly between the limits and pass these
  # to the relevant quantile function.
  qnorm(runif(n, p_low, p_high), mu, sigma)
}

rbern = function(n, p){
  rbinom(n, size=1, p)
}

# function to plot the prior and posterior distribution of each parameter
# must have the 
priorpost_plot = function(prior, post){
  param_names = colnames(prior[[1]])
  for(param in param_names){
    traceplot(post[,param], main=paste("Traceplot of", param))
    
    prior_dens = density(prior[,param][[1]])
    post_dens = density(post[,param][[1]])
    x_lim = range(c(prior_dens$x, post_dens$x))
    y_max= max(c(prior_dens$y, post_dens$y))
    
    plot(1, type='n', xlim=c(0,1), ylim=c(0,y_max), xlab="", ylab="", 
         main=paste("Density of ", param))
    lines(prior_dens, lwd=2, col=myGreen)
    lines(post_dens, lwd=2, col=myBlack)
    legend("topright", lty=1, lwd=2, col=c(myGreen, myBlack), legend=c("prior", "post"))
  }
}

mcmc_plot = function(posterior, prior, title){
  
  for(param in colnames(posterior)[!colnames(posterior)%in%c("threshold", "percent_below")]){
    param_vec = posterior[,paste(param)]
    plot(ts(param_vec), ylab=paste(param), xlab="Index", 
         main=paste(param,"Trace Plot"), cex.lab=1.2, cex.main=1.4 )
    if(length(unique(param_vec))>1) acf(param_vec, lag=20, main="", cex.lab=1.2)
    else plot(NA, xlim=c(0,20), ylim=c(0,1), ylab="ACF", xlab="Lag", cex.lab=1.2)
    
    plot(density(param_vec[!is.na(param_vec)]), lwd=2, col=myBlue(0.5),
         main=paste(param,"Posterior Desntiy"), cex.lab=1.2, cex.main=1.4 )
    lines(density(prior[, paste(param)]), lwd=2, col=myGreen(0.5), )
    
    title(main=title, outer=TRUE, line=-2, cex.main=1.6)
  }
}

# load data from Sarah's repo
# The working directory should be Sarah's repo
raw_data = read.csv("../../output/all_cleaned_sc_data.csv", header=TRUE)

# sets data to be between 0 and 1
sc_data = raw_data
sc_data[,"HET"] = raw_data[,"HET"]/100
sc_data[,"cell_type"] = sub(" ", "", sc_data[,"cell_type"])
sc_data[,"cell_type"] = sub("/", "_", sc_data[,"cell_type"])


# identifying cell types and patients
cell_types_all = unique(sc_data$cell_type)
cell_types = cell_types_all[!(cell_types_all %in% c("4+TEMRA","8+EM_TEMRA"))]

ocell_types = c("34+Progenitor", "Monocyte","4+naive","4+CM","4+EM", "4+TEMRA",
                "8+naive","8+CM","8+EM","8+TEMRA","Bnaive","Bmemory")

pts = unique(sc_data$Patient)
ages = vector("numeric")
for(i in 1:length(pts)){
  ages[paste(pts[i])] = unique(sc_data[sc_data$Patient==pts[i], "age"])
}
names(ages) = pts

pat_cells = list()
for(pat in pts[order(ages)]){
  ctypes = unique(sc_data[sc_data$Patient==pat, "cell_type"])
  pat_cells[[pat]] = ocell_types[ocell_types %in% ctypes]
}

cellpat_post = list()
cellpat_prior = list()
for(pat in pts){
  for(type in ocell_types){
    fp_post = paste0("./OUTPUT/cellpat_bool/",type,"_",pat,"_POST.csv")
    fp_prior = paste0("./OUTPUT/cellpat_bool/",type,"_",pat,"_PRIOR.csv")
    if(file.exists(fp_post)){
      cellpat_post[[paste0(pat,"_",type)]] = read.csv(file=fp_post, header=TRUE)
    }
    if(file.exists(fp_prior)){
      cellpat_prior[[paste0(pat,"_",type)]] = read.csv(file=fp_prior, header=TRUE)
    }
  }
}

for(cellpat in names(cellpat_post)){
  df_post = cellpat_post[[paste(cellpat)]]
  df_prior = cellpat_prior[[paste(cellpat)]]
  cellpat_post[[paste(cellpat)]] = df_post[,!(names(df_post)%in%c("X"))]
  cellpat_prior[[paste(cellpat)]] = df_prior[,!(names(df_prior)%in%c("X"))]
}


temp_vec = vector("numeric", length=length(ocell_types))
names(temp_vec) = ocell_types
median_pi = matrix(NA, nrow=length(ocell_types), ncol=length(pts), 
                   dimnames=list(ocell_types, pts))
for(pat in pts){
  for(cell in ocell_types){
    froot = paste(pat, cell, sep="_")
    if( !is.null(cellpat_post[[froot]][,"pi"])){
      temp_vec[cell] = median(cellpat_post[[froot]][,"pi_infl"])
    } else {
      temp_vec[cell] = NA
    }
  }
  median_pi[,pat] = temp_vec
}

write.table(median_pi, file="median_piProp.csv", sep=", ")

transparency = 0.005

pdf("./PDF/data.pdf", width=12, height=8)
{
  for(pat in pts[order(ages)]){
    par(mfrow=c(2,2))
    for(cell in ocell_types){
      if( length(data_list[[pat]][[cell]])>0){
        hist(data_list[[pat]][[cell]], breaks=seq(0,1,0.01), freq=FALSE,
             col=myGrey(0.4), border=myGrey(0.4), 
             xlab="Mutation Load", ylab="Density", 
             main=paste(pat, cell,"\nAge:",ages[pat])
        )
      }
    }
  } 
}
dev.off()

pdf("./PDF/cellpat_bool/MCMC.pdf", width=12, height=8)
{
  for(cellpat in names(cellpat_post)){
    par(mfrow=c(3,3))
    mcmc_plot(cellpat_post[[cellpat]], cellpat_prior[[cellpat]], title=cellpat)
  }
}
dev.off()

pdf("./PDF/cellpat_bool/model_posteriors.pdf", height=12, width=14.5)
{
  for(cellpat in names(cellpat_post)){
    post = cellpat_post[[paste(cellpat)]]
    prior = cellpat_prior[[paste(cellpat)]]
    
    cellpat_split = strsplit(cellpat, split="_")
    pat = cellpat_split[[1]][1]
    type = cellpat_split[[1]][2]
    if( length(cellpat_split[[1]])!=2){
      type = paste0(cellpat_split[[1]][2],"_",cellpat_split[[1]][3] )
    }
    cellpat_data = sc_data[sc_data$Patient==pat & sc_data$cell_type==type, ]
    
    par(mfrow=c(3,2))
    ymax = max(c(hist(post[,"mu"], plot=FALSE, breaks=seq(0,0.2,0.002))$density, hist(prior[,"mu"], plot=FALSE, breaks=seq(0,0.2,0.002))$density))
    hist(prior[,"mu"], border=myRed(0.3), col=myRed(0.3), 
         main=expression(mu~"Prior & Posterior"), xlab=expression(mu),
         xlim=c(0,0.2), ylim=c(0,ymax), cex.main=1.2, freq=FALSE, breaks=seq(0,0.2,0.002))
    hist(post[,"mu"], col=myBlue(0.3), border=myBlue(0.3),
         add=TRUE, freq=FALSE, breaks=seq(0,0.2,0.002))
    legend("topright", lty=1, lwd=2, col=c(myBlue(0.9), myRed(0.9)), legend=c("Post","Prior"))
    
    xmax = max(c(post[,"sd"]), prior[,"sd"]) + 0.003
    ymax = max(c(hist(post[,"sd"], plot=FALSE, breaks=seq(from=0,to=xmax,by=0.003))$density, hist(prior[,"sd"], plot=FALSE, breaks=seq(0,xmax,by=0.003))$density))
    hist(prior[,"sd"], freq=FALSE, border=myRed(0.3), col=myRed(0.3), 
         breaks=seq(0,xmax,by=0.003), main=expression(sigma~"Prior & Posterior"),
         xlim=c(0,0.3), ylim=c(0,ymax),
         cex.main=1.2, xlab=expression(sigma))
    hist(post[,"sd"], freq=FALSE, breaks=seq(0,xmax,by=0.003), 
         border=myBlue(0.3), col=myBlue(0.3), add=TRUE)
    legend("topright", lty=1, lwd=2, col=c(myBlue(0.9), myRed(0.9)), legend=c("Post","Prior"))
    
    # ymax = max(c(hist(post[,"pi0"], plot=FALSE, breaks=seq(0,1,0.01))$density, hist(prior[,"pi0"], plot=FALSE, breaks=seq(0,1,0.01))$density))
    # hist(prior[,"pi0"], border=myRed(0.3), col=myRed(0.3), xlim=c(0,1), ylim=c(0,ymax), 
    #      main=expression(pi[0]~"Prior & Posterior"), xlab=expression(pi[0]), cex.main=1.2,
    #     freq=FALSE, breaks=seq(0,1,0.01))
    # hist(post[,"pi0"], lwd=2, col=myBlue(0.3), border=myBlue(0.3), 
    #      add=TRUE, freq=FALSE, breaks=seq(0,1,0.01))
    # legend("topright", lty=1, lwd=2, col=c(myBlue(0.6), myRed(0.6)), legend=c("Post","Prior"))
    # 
    # ymax = max(c(hist(post[,"pi"], plot=FALSE, breaks=seq(0,1,0.01))$density, hist(prior[,"pi"], plot=FALSE, breaks=seq(0,1,0.01))$density))
    # hist(prior[,"pi"], border=myRed(0.3), col=myRed(0.3), xlim=c(0,1), ylim=c(0,ymax), 
    #      main=expression(pi~"Prior & Posterior"), xlab=expression(pi), cex.main=1.2,
    #      freq=FALSE, breaks=seq(0,1,0.01))
    # hist(post[,"pi"], lwd=2, col=myBlue(0.3), border=myBlue(0.3), 
    #      add=TRUE, freq=FALSE, breaks=seq(0,1,0.01))
    # legend("topright", lty=1, lwd=2, col=c(myBlue(0.6), myRed(0.6)), legend=c("Post","Prior"))
    # 
    ymax = max(c(hist(post[,"pi_infl"], plot=FALSE, breaks=seq(0,1,0.01))$density, hist(prior[,"pi_infl"], plot=FALSE, breaks=seq(0,1,0.01))$density))
    hist(prior[,"pi_infl"], border=myRed(0.3), col=myRed(0.3), xlim=c(0,1), ylim=c(0,ymax), 
         main=expression(pi~"Prior & Posterior"), xlab=expression(tilde(pi)),
         cex.main=1.2,  freq=FALSE, breaks=seq(0,1,0.01))
    hist(post[,"pi_infl"], lwd=2, col=myBlue(0.3), border=myBlue(0.3), 
         add=TRUE, freq=FALSE, breaks=seq(0,1,0.01))
    legend("topright", lty=1, lwd=2, col=c(myBlue(0.6), myRed(0.6)), legend=c("Post","Prior"))
    
    ndx_post = !is.na(post[,"threshold"])
    ndx_prior = !is.na(prior[,"threshold"])
    {
    if( sum(ndx_post)>0 & sum(ndx_prior)>0 ){
      hist(post[ndx_post,"threshold"], border=myBlue(0.3), col=myBlue(0.3),
           main="threshold Prior & Post", freq=FALSE, breaks=seq(0,1,0.01),
           xlab=paste("Prior:",sum(ndx_prior),"Post:",sum(ndx_post)),
           xlim=c(0,1), cex.main=1.2)
      hist(prior[ndx_prior,"threshold"], border=myRed(0.3), col=myRed(0.3),
           breaks=seq(0,1,0.01), add=TRUE, freq=FALSE )
    } else if(sum(ndx_post)>0 & sum(ndx_prior)==0) {
      hist(post[ndx_post,"threshold"], border=myBlue(0.3), col=myBlue(0.3),
           main="threshold Prior & Post", freq=FALSE, breaks=seq(0,1,0.01),
           xlab=paste("Prior:",sum(ndx_prior),"Post:",sum(ndx_post)),
           xlim=c(0,1), cex.main=1.2)
    } else if(sum(ndx_post)==0 & sum(ndx_prior)>0) {
      hist(prior[ndx_prior,"threshold"], border=myRed(0.3), col=myRed(0.3),
           main="threshold Prior & Post", freq=FALSE, breaks=seq(0,1,0.01),
           xlab=paste("Prior:",sum(ndx_prior),"Post:",sum(ndx_post)),
           xlim=c(0,1), cex.main=1.2)
    } else {
      hist(rep(-99,1e2), freq=FALSE, xlim=c(0,1), ylim=c(0,1),
           main="threshold Prior & Post",
           xlab=paste("Prior:",sum(ndx_prior),"Post:",sum(ndx_post)),
          cex.main=1.2)
    }
    }
    legend("topright", lty=1, lwd=2, col=c(myBlue(0.9), myRed(0.9)), legend=c("Post","Prior"))

    # compone_dens = density(post[,"pred.1."])
    # data_dens = hist(cellpat_data$HET, breaks=seq(0,1,0.01), plot=FALSE)$density
    # ymax = max(c(data_dens, compone_dens$y))
    # hist(cellpat_data$HET, breaks=seq(0,1,0.01), cex.main=1.2, xlim=c(0,1), ylim=c(0,ymax),
    #      xlab="Mutation Load", col=myBlack(0.2), border=myBlack(0.2), 
    #      main="Component Posterior Predicitve", freq=FALSE)
    # lines(density(post[,"pred.1."]), lwd=2, col=myBlue(0.9))
    # lines(density(post[,"pred.2."]), lwd=2, col=myRed(0.9))
    # legend("topright", lty=1, lwd=2, col=c(myBlue(0.9), myRed(0.9)), legend=c("Spike", "Slab"))
    # 
    
    hist(cellpat_data$HET, breaks=seq(0,1,0.01), freq=FALSE, 
         main="Predictive Posterior and threshold quantiles", xlab="Mutation Load",
         col=myBlack(0.2), border=myBlack(0.2))
    hist(post[,"Ypred"], breaks=seq(0,1,by=0.01), add=TRUE, freq=FALSE,
         col=myPurple(0.4), border=myPurple(0.4))
    thresh_qnts = quantile(post[!is.na(post[,"threshold"]),"threshold"], c(0.025,0.5,0.975))
    abline(v=thresh_qnts, lty=c(2,1,2), col=myGreen(0.9), lwd=c(1,2,1))
    legend("topright", lty=1, lwd=2, col=c(myBlue(0.9), myGreen(0.9)), legend=c("Pred", "thresh"))
    
    title(main=paste(cellpat), line=-2, outer=TRUE, cex.main=2)
  }
}
dev.off()

pdf("./PDF/cellpat_bool/pi_prior.pdf", width=12, height=8)
{
  n=1e4
  pi0 = runif(n)
  z0 = rbinom(n,1,pi0)
  pis = vector("numeric", n)
  for(i in 1:n){
    pis[i] = ifelse(z0[i]==0, runif(1), 0)
  }
  hist(pis, breaks=seq(0,1,0.01), freq=FALSE, 
       col=myGrey(0.6), border=myGrey(0.6), 
       xlab=expression(pi), main=expression(pi~"Prior Distribution"))
}
dev.off()

spike_prop = list()
for(pat in pts){
  spike_prop[[pat]] = list()
  spike_prop[[pat]]$Prior = cellpat_prior[[paste0(pat, "_Monocyte")]][,"pi_infl"]
  for(type in ocell_types){
    cellpat = paste0(pat,"_",type)
    if(!is.null(cellpat_post[[cellpat]])){
      spike_prop[[pat]][[type]] = cellpat_post[[cellpat]][,"pi_infl"]
    }
  }
}
pdf("./PDF/cellpat_bool/spike_prop.pdf", height=8, width=12)
{
  for(pat in names(spike_prop)){
    cells = gsub("genitor", "", gsub("cyte", "", names(spike_prop[[pat]])))
    cell_splits = c(grep("Prior", cells), grep("Mo", cells), max(grep("4", cells)), max(grep("8", cells))) + c(0.5,0.5,0.5,0.5)
    par(mfrow=c(1,1))
    stripchart(spike_prop[[pat]], vertical=TRUE, method="jitter", pch=20, col=c(myGreen(transparency), rep(myBlue(transparency), length(spike_prop[[pat]])-1)), 
               xlab="Cell Type", ylab="Proportion from spike", main=paste(pat, "Age:", ages[pat]), 
               ylim=c(0,1), cex.lab=1.4, cex.main=1.4, 
               group.names = cells)
    abline(v=cell_splits, lwd=6, lty=2, col=myGrey(0.7))
  }
}
dev.off()

jiggle = 0.4
pi_age_list = list()
{
  for(cell in ocell_types){
    pi_age_list[[cell]] = list()
    for(pat in pts){
      pi_post = cellpat_post[[paste0(pat,"_",cell)]][,"pi_infl"]
      N = length(pi_post)
      pi_age_list[[cell]][[pat]][["HET"]] = pi_post
      pi_age_list[[cell]][[pat]][["age"]] = rep(ages[pat], N) + rnorm(N,0,jiggle)
    }
  }
}

pdf("./PDF/cellpat_bool/spike_prop_withAge.pdf", width=10, height=8)
{
  for(cell in names(pi_age_list)){
    plot(1, type='n', xlim=c(10,60), ylim=c(0,1), 
         xlab="Age", ylab="Spike proportion", main=paste0(cell,"\n Spike Proportion"), 
         cex.lab=1.2, cex.main=1.4)
    i = 1
    for(pat in pts){
      points(pi_age_list[[cell]][[pat]]$age, pi_age_list[[cell]][[pat]]$HET, 
             pch=20, col=myCols(transparency, pts)[pat])
      i = i + 1
    }
    legend("topright", legend=pts, pch=20, col=myCols(0.9, pts) )
    
  }
}
dev.off()

# list of cell type pairs we wish to compare
# constructed in a dumb way...
comp_pairs = vector("character")
{
  i = 1
  for(cell1 in c("34+Progenitor", "Monocyte", "4+naive", "4+CM", "4+EM", "8+naive","8+CM","8+EM", "Bnaive")){
    if(cell1=="34+Progenitor"){
      for(cell2 in ocell_types[!cell_types%in%c("34+Progenitor")]){
        comp_pairs[i] = paste0(cell1, "_", cell2)
        i = i + 1
      }
    }
    if(cell1=="4+naive"){
      for(cell2 in c("4+CM", "4+EM", "4+TEMRA")){
        comp_pairs[i] = paste0(cell1, "_", cell2)
        i = i + 1
      }
    }
    if(cell1=="4+CM"){
      for(cell2 in c("4+EM")){
        comp_pairs[i] = paste0(cell1, "_", cell2)
        i = i + 1
      }
    }
    if(cell1=="4+EM"){
      comp_pairs[i] = paste0(cell1, "_4+TEMRA")
      i = i + 1
    }
    if(cell1=="8+naive"){
      for(cell2 in c("8+CM", "8+EM", "8+TEMRA")){
        comp_pairs[i] = paste0(cell1, "_", cell2)
        i = i + 1
      }
    }
    if(cell1=="8+CM"){
      for(cell2 in c("8+EM", "8+TEMRA")){
        comp_pairs[i] = paste0(cell1,"_",cell2)
        i = i + 1
      }
    }
    if(cell1=="8+EM"){
      comp_pairs[i] = paste0(cell1, "_8+TEMRA")
      i = i + 1
    }
    if(cell1 == "Bnaive"){
      comp_pairs[i] = paste0(cell1, "_Bmemory")
      i = i + 1
    }
  }
}

# construct lists
picomp_list = list()
picomp_list_qnts = list()
for(cellpair in comp_pairs){
  pair = strsplit(cellpair, split="_")[[1]]
  cell_one = pair[1]
  cell_two = pair[2]
  
  picomp_list[[cellpair]] = list()
  for(pat in pts[order(ages)]){
    cellone_pi = cellpat_post[[paste0(pat,"_",cell_one)]][,"pi_infl"]
    celltwo_pi = cellpat_post[[paste0(pat,"_",cell_two)]][,"pi_infl"]
    diff = cellone_pi - celltwo_pi
    picomp_list[[cellpair]][[pat]] = diff
    picomp_list_qnts[[cellpair]]$med[pat] = round(quantile(diff, c(0.5)),3)
    picomp_list_qnts[[cellpair]]$lower[pat] = round(hdi(diff, credMass=0.95)[1], 3)
    picomp_list_qnts[[cellpair]]$upper[pat] = round(hdi(diff, credMass=0.95)[2], 3)
    prob = sum(diff>0)/length(diff)
    nullprob = 1/length(diff)
    picomp_list_qnts[[cellpair]]$prob_val[pat] = ifelse(prob==0, nullprob, prob)
    picomp_list_qnts[[cellpair]]$prob_char[pat] = ifelse(prob==0, paste("<",nullprob), paste("=",prob))
  }
}

# for each patient a list of combinations in the form
# for(pat in pts){
#   for(celltwo in celltwos){
#     list[[pat]][[cellone]] = *vector of celltwos* 
cellpairs = list() # list to store (desired) output
{
  for(pat in pts[order(ages)]){ # cycle through patients
    cellpairs[[pat]] = list() # create patient list
    
    for(cellone in pat_cells[[pat]][!grepl("34+", pat_cells[[pat]])]){ # cycle through all permutations of comparison
      # even if the cells aren't recorded in patient
      tmp = vector("numeric")
      for(celltwo in pat_cells[[pat]]){
        if(paste0(cellone, "_", celltwo) %in% comp_pairs){
          tmp = c(tmp, celltwo)
        }
      }
      if(length(tmp)>0) cellpairs[[pat]][[cellone]] = tmp
    }
  }
}

sig_col = function( sig, alpha=0.7 ){
  if(is.null(sig))  return(NULL)
  else if(sig) return( myRed(alpha) ) 
  else if(!sig) return( myCerulean(alpha) ) 
}


# for all cell combos tested within one patient 
# a list indicating whether that difference is significant 
picomp_sig = list()
for(pat in pts){
  temp_list = list()
  for(cellone in names(cellpairs[[pat]])){
    for(celltwo in cellpairs[[pat]][[cellone]]){
      cellpair = paste0(cellone,"_",celltwo)
      if(cellone %in% pat_cells[[pat]] & celltwo %in% pat_cells[[pat]]){
        lower = picomp_list_qnts[[cellpair]]$lower[pat]
        upper = picomp_list_qnts[[cellpair]]$upper[pat]
        temp_list[[cellpair]] = !(lower<=0 & 0<=upper)
        
        diff = picomp_list[[cellpair]][[pat]]
        n = length(diff)
        prob = sum(diff>0)/n
        
        # temp_list[[cellpair]] = (prob<0.05 | prob>0.95)
      }
    }
  }
  picomp_sig[[pat]] = temp_list
}

pdf("./PDF/cellpat_bool/spikeprop_branch_comp.pdf", width=12, height=8)
{
  for(pat in pts[order(ages)]){
    spike_prop[[pat]]$Prior = NULL
    # plot all posterior draws for pi for each patient in cell "age" order
    cells = names(spike_prop[[pat]])
    cell_splits = c(grep("Mo", cells, fixed=TRUE), max(grep("4", cells, fixed=TRUE)), max(grep("8", cells, fixed=TRUE))) + c(0.5,0.5,0.5)
    
    par(mar=c(12,5,9,3), mgp=c(3,1,0), xpd=TRUE)
    groupnames = sub("ory","", sub("genitor", "", sub("cyte","", names(spike_prop[[pat]]))))
    stripchart(spike_prop[[pat]], vertical=TRUE, method="jitter", pch=20, col=myBlack(0.002), 
               xlab="", ylab="Proportion from spike",
               ylim=c(-0,1), cex.lab=1.4, group.names=groupnames)
    text(x=length(cells)/2, y=1.45, labels=paste(pat, "Age:", ages[pat]), cex=1.4)
    text(x=length(cells)/2, y=-0.55, labels="Cell Type", cex=1.4)
    for(i in seq_along(cell_splits)){
      lines(x=rep(cell_splits[i], 2), y=c(0,1), lty="dashed", col=myGrey(0.7), lwd=4) # adds lines to split by cell group
    }
    
    cmp_list = cellpairs[[pat]] # list of cell comparisons made
    
    for(cellone in names(cmp_list)[!grepl("34+", names(cmp_list))]){ # add the comparison lines
      
      cell_index = 1:length(spike_prop[[pat]] )
      names(cell_index) = names(spike_prop[[pat]])
      
      cellone_top = quantile(spike_prop[[pat]][[cellone]], 0.99)
      cellone_bottom = quantile(spike_prop[[pat]][[cellone]], 0.01)
      ytop = 1.15
      ybottom = -0.35
      side = 0.07
      for(celltwo in cmp_list[[cellone]]){
        
        if(grepl("naive", cellone)){
          signif = picomp_sig[[pat]][[paste0(cellone,"_",celltwo)]]
          lines(c(cell_index[cellone], cell_index[cellone]), c(ytop-side, ytop), 
                lwd=5, col=sig_col(signif))
          
          lines(c(cell_index[cellone], cell_index[celltwo]), c(ytop, ytop), 
                lwd=5, col=sig_col(signif))
          
          lines(c(cell_index[celltwo], cell_index[celltwo]), c(ytop, ytop-side ), 
                lwd=5, col=sig_col(signif)) 
          ytop = ytop + 0.1
        } else {
          ybottom = ifelse( cellone=="8+EM" | cellone=="4+EM", -0.25, ybottom)
          ybottom = ifelse( cellone=="4+CM" && length(cellpairs[[pat]][[cellone]])==1, -0.25, ybottom)
          
          signif = picomp_sig[[pat]][[paste0(cellone,"_",celltwo)]]
          
          celltwo_top = quantile(spike_prop[[pat]][[celltwo]], 0.99)
          lines(c(cell_index[cellone], cell_index[cellone]), c(ybottom+side, ybottom), 
                lwd=5, col=sig_col(signif))
          
          lines(c(cell_index[cellone], cell_index[celltwo]), c(ybottom, ybottom), 
                lwd=5, col=sig_col(signif))
          
          lines(c(cell_index[celltwo], cell_index[celltwo]), c(ybottom, ybottom+side), 
                lwd=5, col=sig_col(signif))
          ybottom = ybottom - 0.1
        }
      }
    }
  }
}
dev.off()

comp_pairs2 = c("34+Progenitor_Monocyte", "Monocyte_Bnaive", "Monocyte_Bmemory")
picomp_list2 = list()
picomp_list_qnts2 = list()
for(cellpair in comp_pairs2){
  pair = strsplit(cellpair, split="_")[[1]]
  cell_one = pair[1]
  cell_two = pair[2]
  
  picomp_list2[[cellpair]] = list()
  for(pat in pts[order(ages)]){
    cellone_pi = cellpat_post[[paste0(pat,"_",cell_one)]][,"pi_infl"]
    celltwo_pi = cellpat_post[[paste0(pat,"_",cell_two)]][,"pi_infl"]
    diff = cellone_pi - celltwo_pi
    picomp_list2[[cellpair]][[pat]] = diff
    picomp_list_qnts2[[cellpair]]$med[pat] = round(quantile(diff, c(0.5)),3)
    picomp_list_qnts2[[cellpair]]$lower[pat] = round(hdi(diff, credMass=0.95)[1], 3)
    picomp_list_qnts2[[cellpair]]$upper[pat] = round(hdi(diff, credMass=0.95)[2], 3)
    prob = sum(diff>0)/length(diff)
    nullprob = 1/length(diff)
    picomp_list_qnts2[[cellpair]]$prob_val[pat] = ifelse(prob==0, nullprob, prob)
    picomp_list_qnts2[[cellpair]]$prob_char[pat] = ifelse(prob==0, paste("<",nullprob), paste("=",prob))
  }
}

cellpairs2 = list() # list to store (desired) output
{
  for(pat in pts[order(ages)]){ # cycle through patients
    cellpairs2[[pat]] = list() # create patient list
    
    for(cellone in pat_cells[[pat]]){ # cycle through all permutations of comparison
      # even if the cells aren't recorded in patient
      tmp = vector("numeric")
      for(celltwo in pat_cells[[pat]]){
        if(paste0(cellone, "_", celltwo) %in% comp_pairs2){
          tmp = c(tmp, celltwo)
        }
      }
      if(length(tmp)>0) cellpairs2[[pat]][[cellone]] = tmp
    }
  }
}

picomp_sig2 = list()
for(pat in pts){
  temp_list = list()
  for(cellone in names(cellpairs2[[pat]])){
    for(celltwo in cellpairs2[[pat]][[cellone]]){
      cellpair = paste0(cellone,"_",celltwo)
      if(cellone %in% pat_cells[[pat]] & celltwo %in% pat_cells[[pat]]){
        lower = picomp_list_qnts2[[cellpair]]$lower[pat]
        upper = picomp_list_qnts2[[cellpair]]$upper[pat]
        temp_list[[cellpair]] = !(lower<=0 & 0<=upper)
        
        diff = picomp_list2[[cellpair]][[pat]]
        n = length(diff)
        prob = sum(diff>0)/n
        
        # temp_list[[cellpair]] = (prob<0.05 | prob>0.95)
      }
    }
  }
  picomp_sig2[[pat]] = temp_list
}

pdf("./PDF/cellpat_bool/something_robust.pdf", width=12, height=8)
{ 
  for(pat in pts[order(ages)]){
    spike_prop[[pat]]$Prior = NULL
    # plot all posterior draws for pi for each patient in cell "age" order
    cells = names(spike_prop[[pat]])
    cell_splits = c(grep("Mo", cells, fixed=TRUE), max(grep("4", cells, fixed=TRUE)), max(grep("8", cells, fixed=TRUE))) + c(0.5,0.5,0.5)
    
    par(mar=c(12,5,9,3), mgp=c(3,1,0), xpd=TRUE)
    groupnames = sub("ory","", sub("genitor", "", sub("cyte","", names(spike_prop[[pat]]))))
    stripchart(spike_prop[[pat]], vertical=TRUE, method="jitter", pch=20, col=myBlack(0.002), 
               xlab="", ylab="Proportion from spike",
               ylim=c(-0,1), cex.lab=1.4, group.names=groupnames)
    text(x=length(cells)/2, y=1.45, labels=paste(pat, "Age:", ages[pat]), cex=1.4)
    text(x=length(cells)/2, y=-0.55, labels="Cell Type", cex=1.4)
    for(i in seq_along(cell_splits)){
      lines(x=rep(cell_splits[i], 2), y=c(0,1), lty="dashed", col=myGrey(0.7), lwd=4) # adds lines to split by cell group
    }
    
    cmp_list = cellpairs2[[pat]] # list of cell comparisons made
    
    for(cellone in names(cmp_list)){ # add the comparison lines
      
      cell_index = 1:length(spike_prop[[pat]] )
      names(cell_index) = names(spike_prop[[pat]])
      
      cellone_top = quantile(spike_prop[[pat]][[cellone]], 0.99)
      cellone_bottom = quantile(spike_prop[[pat]][[cellone]], 0.01)
      ytop = 1.15
      ybottom = -0.35
      side = 0.07
      
      sig1 = picomp_sig2[[pat]][["34+Progenitor_Monocyte"]]
      lines(c(cell_index["34+Progenitor"], cell_index["34+Progenitor"]), c(ybottom+side, ybottom), 
            lwd=5, col=sig_col(sig1))
      lines(c(cell_index["34+Progenitor"], cell_index["Monocyte"]), c(ybottom, ybottom), 
            lwd=5, col=sig_col(sig1))
      lines(c(cell_index["Monocyte"], cell_index["Monocyte"]), c(ybottom, ybottom+side), 
            lwd=5, col=sig_col(sig1)) 
      
      for(celltwo in cmp_list[["Monocyte"]]){
          signif = picomp_sig2[[pat]][[paste0(cellone,"_",celltwo)]]
          
          celltwo_top = quantile(spike_prop[[pat]][[celltwo]], 0.99)
          lines(c(cell_index[cellone], cell_index[cellone]), c(ytop-side, ytop), 
                lwd=5, col=sig_col(signif))
          
          lines(c(cell_index[cellone], cell_index[celltwo]), c(ytop, ytop), 
                lwd=5, col=sig_col(signif))
          
          lines(c(cell_index[celltwo], cell_index[celltwo]), c(ytop, ytop-side), 
                lwd=5, col=sig_col(signif))
          ytop = ytop + 0.1
      }
    }
  }
}
dev.off()

#####
## IS THE SPIKE DIFFERENT FROM ZERO?
#####
sig_col = function( sig, alpha=0.7 ){
  if(is.null(sig))  return(NULL)
  else if(sig) return( myRed(alpha) ) 
  else if(!sig) return( myCerulean(alpha) ) 
}

lengths = function(list){
  vec = vector("numeric")
  for(name in names(list)){
    vec[name] = length(list[[name]])
  }
  vec
}

spike_list = list()
spike_qnts = list()
for(pat in pts[order(ages)]){
  spike_list[[pat]] = list()
  for(cell in pat_cells[[pat]]){
    pipost = cellpat_post[[paste(pat,cell,sep="_")]][,"pi_infl"]
    if(!is.null(pipost)){
      Z_spike = rbern(n=length(pipost), p=pipost)
      spike_cellpat = cellpat_post[[paste(pat, cell, sep="_")]][Z_spike==1,"pred.1."]
      spike_list[[pat]][[cell]] = spike_cellpat
      spike_qnts[[pat]][[cell]]$hdi = round(hdi(spike_cellpat, credMass=0.90), 3)
      spike_qnts[[pat]][[cell]]$qnt = round(quantile(spike_cellpat, 0.05), 3)
      spike_qnts[[pat]][[cell]]$sig_hdi = 0<spike_qnts[[pat]][[cell]]$hdi[1]
      spike_qnts[[pat]][[cell]]$sig_qnt = 0<spike_qnts[[pat]][[cell]]$qnt
      spike_qnts[[pat]][[cell]]$Ntotal = length(pipost)
      spike_qnts[[pat]][[cell]]$Nspike = sum(Z_spike)
    } else {
      spike_list[[pat]][[cell]] = NA
      spike_qnts[[pat]][[cell]] = NULL
    }
  }
}

pdf("./PDF/cellpat_bool/spike_post_qnts.pdf", width=12, height=8)
{
  par(mar=c(6,5,5,3))
  for(pat in pts[order(ages)]){
    cellnames = gsub("cyte", "", gsub( "genitor", "", names(spike_list[[pat]])))
    cell_index = 1:length(spike_list[[pat]])
    names(cell_index) = names(spike_list[[pat]])
    stripchart(spike_list[[pat]], vertical=TRUE, pch=20, cex=0.7, 
               col=myBlack(transparency), method="jitter", ylim=c(0,1), 
               main=paste(pat,"Age:", ages[pat]), ylab="Predictive posterior of spike component",
               xlab="Cell Type", cex.lab=1.6, cex.main=1.8,
               group.names=cellnames)
    
    for(cell in ocell_types){
      if(!is.null(spike_qnts[[pat]][[cell]])){
        qnts = spike_qnts[[pat]][[cell]]
        arrows(x0=cell_index[cell], y0=qnts$hdi[1],
               x1=cell_index[cell], y1=qnts$hdi[2],
               angle=90, lty=1, lwd=3, length=0.12,
               col=sig_col(qnts$sig_qnt, alpha=qnts$Nspike/qnts$Ntotal), 
               code=3)
      }
    }
    text(cell_index, y=rep(1, length(cell_index)), 
         labels=paste0("N=",lengths(spike_list[[pat]])) )
  }
}
dev.off()

