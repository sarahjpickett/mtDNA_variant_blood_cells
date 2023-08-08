library("HDInterval")
library("coda")
myBlack = function(alpha) rgb(0,0,0, alpha)
myGrey = function(alpha) rgb(66/255,66/255,66/255, alpha)
myRed = function(alpha) rgb(1,0,0, alpha)
myBlue = function(alpha) rgb(0,0,1, alpha)

# creates folders to save output in 
dir.create("./PDF", showWarnings = FALSE)

raw_data = read.csv("../../output/all_cleaned_sc_data.csv", header=TRUE)

# order cell types
ocell_types = c("34+Progenitor", "Monocyte","4+naive","4+CM","4+EM", "4+TEMRA",
                "8+naive","8+CM","8+EM","8+TEMRA","Bnaive","Bmemory")

# sets data to be between 0 and 1
sc_data = raw_data
sc_data[,"HET"] = raw_data[,"HET"]/100
sc_data[,"cell_type"] = sub(" ", "", sc_data[,"cell_type"])
sc_data[,"cell_type"] = sub("/", "_", sc_data[,"cell_type"])

# identifying cell types and patients
cell_types_all = unique(sc_data$cell_type)
cell_types = cell_types_all[!(cell_types_all %in% c("4+TEMRA","8+EM_TEMRA"))]

# vector of cell types in order, as desired
ocell_types = c("34+Progenitor", "Monocyte","4+naive","4+CM","4+EM", "4+TEMRA",
                "8+naive","8+CM","8+EM","8+TEMRA","Bnaive","Bmemory")

# vector of cell type pairs we wish to compare
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

# change the patient labels for the plots
pts_change = c("P10", "P15", "P17", "P18", "P19", "P22")
names(pts_change) = c("R013", "R025", "R027", "R028", "R044", "R0149")
pts = unique(sc_data$Patient)
# vector of patient ages
ages = vector("numeric")
for(i in 1:length(pts)){
  ages[paste(pts[i])] = unique(sc_data[sc_data$Patient==pts[i], "age"])
}
names(ages) = pts

# list of cell types present for each patient
pat_cells = list()
for(pat in pts[order(ages)]){
  ctypes = unique(sc_data[sc_data$Patient==pat, "cell_type"])
  pat_cells[[pat]] = ocell_types[ocell_types %in% ctypes]
}

# load inference output for each patient - cell type  combination and store
# in a list
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

# remove unwanted columns
for(cellpat in names(cellpat_post)){
  df_post = cellpat_post[[paste(cellpat)]]
  df_prior = cellpat_prior[[paste(cellpat)]]
  cellpat_post[[paste(cellpat)]] = df_post[,!(names(df_post)%in%c("X"))]
  cellpat_prior[[paste(cellpat)]] = df_prior[,!(names(df_prior)%in%c("X"))]
}
# matrices to store lower and upper quantiles and medians for each cell+patient
# for the posterior proportion of near zero mutation load
lq_pi = matrix(NA, nrow=length(ocell_types), ncol=length(pts), 
               dimnames=list(ocell_types, pts))
median_pi = lq_pi
uq_pi = lq_pi
# calculate and store quantiles and medians
for(pat in pts){
  temp_vec = vector("numeric", length=length(ocell_types))
  names(temp_vec) = ocell_types
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

# save the medians as a csv file
write.table(round(median_pi, 4), file="median_nearZeroProportion.csv", sep=", ")


# list of cell type pairs we wish to compare
# constructed in an inefficient way
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

# construct lists of difference of distributions and high density intervals 
# to test for significant differences
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

# function to return colour to indicate significance
sig_col = function( sig, alpha=1.0 ){
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

# function to create colour name and transparency level
col_alpha = function(color_name, alpha=1.0) {
  rgb.val = col2rgb(color_name)
  t.col = rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = alpha*255 )
  return(t.col)
}

transparency = 0.015
# create vector of colours for each cell type
col_vec = vector()
col_vec["34+Progenitor"] = col_alpha("#BC3C29FF", transparency)
col_vec["Monocyte"] = col_alpha("#0072B5FF", transparency)
col_vec["4+naive"] = col_alpha("#E18727FF", transparency)
col_vec["4+CM"] = col_alpha("#E18727FF", transparency)
col_vec["4+EM"] = col_alpha("#E18727FF", transparency)
col_vec["4+TEMRA"] = col_alpha("#E18727FF", transparency)
col_vec["8+naive"] = col_alpha("#20854EFF", transparency)
col_vec["8+CM"] = col_alpha("#20854EFF", transparency)
col_vec["8+EM"] = col_alpha("#20854EFF", transparency)
col_vec["8+TEMRA"] = col_alpha("#20854EFF", transparency)
col_vec["Bnaive"] = col_alpha("#7876B1FF", transparency)
col_vec["Bmemory"] = col_alpha("#7876B1FF", transparency)

# list of lists for the near zero proportion of mutation load
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

# plot strip chart of proportion of cells with a near zero mutation load
png("./PDF/nearZero_m3242AG_mode.png", width=12, height=15, unit="in", res=300)
{
  par(mfrow=c(3,2), mar=c(9,9,8,3), mgp=c(3,1,0),  xpd=TRUE, 
      cex.main=2, cex.lab=2, cex.axis=1.5)
  
  for(pat in pts[order(ages)]){
    spike_prop[[pat]]$Prior = NULL
    # plot all posterior draws for pi for each patient in cell "age" order
    cells = names(spike_prop[[pat]])
    cell_splits = c(grep("Mo", cells, fixed=TRUE), max(grep("4", cells, fixed=TRUE)), max(grep("8", cells, fixed=TRUE))) + c(0.5,0.5,0.5)

    groupnames = names(spike_prop[[pat]])
    stripchart(spike_prop[[pat]], vertical=TRUE, method="jitter", pch=20, col=col_vec[groupnames], 
               xlab="", ylab="",
               ylim=c(-0,1), group.names=groupnames, yaxt="n", xaxt="n")
    
    for(i in 1:length(groupnames)){
      if( (pat=="R044" && groupnames[i]=="8+naive") | (pat=="R013" && (groupnames[i]=="Bmemory"||groupnames[i]=="Bnaive")) ){
        hdi_obj = hdi(density(spike_prop[[pat]][[groupnames[i]]]), allowSplit=TRUE)

        for(j in 1:nrow(hdi_obj)){
          den = density(spike_prop[[pat]][[groupnames[i]]])
          hdi_region = hdi_obj[j,]
          peak = den$x[hdi_region[1]<=den$x & den$x<=hdi_region[2]][ which.max(den$y[hdi_region[1]<=den$x & den$x<=hdi_region[2]]) ]

          lines(c(i-0.2, i+0.2), rep(peak,2), col=myBlack(1.0), lwd=1)
          # lines(c(i-0.15, i+0.15), rep(hdi_region[1],2),
          #       col=myBlack(0.5), lwd=1, lty=3)
          # lines(c(i-0.15, i+0.15), rep(hdi_region[2],2),
          #       col=myBlack(0.5), lwd=1, lty=3)
        }

      } else {

        hdi_region = hdi( spike_prop[[pat]][[groupnames[i]]] )
        den = density(spike_prop[[pat]][[groupnames[i]]])

        if( hdi_region[1] == hdi_region[2] ){
          peak = hdi_region[1]
        } else {
          peak = den$x[ hdi_region[1]<= den$x & den$x<=hdi_region[2]][ which.max(den$y[hdi_region[1]<=den$x & den$x<=hdi_region[2]]) ]
        }

        lines(c(i-0.2, i+0.2), rep(peak,2),
              col=myBlack(1.0), lwd=1)
        # lines(c(i-0.15, i+0.15), rep(hdi_region[1],2),
        #       col=myBlack(0.5), lwd=1, lty=3)
        # lines(c(i-0.15, i+0.15), rep(hdi_region[2],2),
        #       col=myBlack(0.5), lwd=1, lty=3)
      }
    }

    axis(2, at=0:10/10, cex=1.5)
    axis(1, at=1:length(spike_prop[[pat]]), labels=FALSE)
    if( pat=="R0149" ){
      text(x=1:length(spike_prop[[pat]]), y=-0.24, labels=groupnames, srt=35, adj=1, cex=1.5 )
    } else if( pat=="R027" ){
      text(x=1:length(spike_prop[[pat]]), y=-0.14, labels=groupnames, srt=35, adj=1, cex=1.5 )
    } else {
      text(x=1:length(spike_prop[[pat]]), y=-0.19, labels=groupnames, srt=35, adj=1, cex=1.5 )
    }
    
    text(x=-1.5, y=0.5, labels="Proportion of cells with near zero\nlevel of m.3243A>G", 
         cex=2, srt=90)
    text(x=length(cells)/2, y=1.3, labels=paste(pts_change[pat], "Age:", ages[pat]), cex=2)
    # text(x=length(cells)/2, y=-0.5, labels="Cell Type", cex=2)
    
    for(i in seq_along(cell_splits)){
      lines(x=rep(cell_splits[i], 2), y=c(0,1), lty="dashed", col=myGrey(0.7), lwd=3) # adds lines to split by cell group
    }
    
    cmp_list = cellpairs[[pat]] # list of cell comparisons made
    
    for(cellone in names(cmp_list)[!grepl("34+", names(cmp_list))]){ # add the comparison lines
      
      cell_index = 1:length(spike_prop[[pat]] )
      names(cell_index) = names(spike_prop[[pat]])
      
      cellone_top = quantile(spike_prop[[pat]][[cellone]], 0.99)
      cellone_bottom = quantile(spike_prop[[pat]][[cellone]], 0.01)
      ytop = 1.1
      ybottom = -0.14
      side = 0.03
      for(celltwo in cmp_list[[cellone]]){
        
        if(grepl("naive", cellone)){
          signif = picomp_sig[[pat]][[paste0(cellone,"_",celltwo)]]
          if(signif){
            lines(c(cell_index[cellone], cell_index[cellone]), c(ytop-side, ytop), 
                  lwd=2, col=sig_col(signif))
            
            lines(c(cell_index[cellone], cell_index[celltwo]), c(ytop, ytop), 
                  lwd=2, col=sig_col(signif))
            
            lines(c(cell_index[celltwo], cell_index[celltwo]), c(ytop, ytop-side ), 
                  lwd=2, col=sig_col(signif)) 
            ytop = ytop + 0.05
          }
        } else {
          ybottom = ifelse( cellone=="8+EM" | cellone=="4+EM", -0.14, ybottom)
          ybottom = ifelse( cellone=="4+CM" && length(cellpairs[[pat]][[cellone]])==1, -0.14, ybottom)
          
          signif = picomp_sig[[pat]][[paste0(cellone,"_",celltwo)]]
          if(signif){
            celltwo_top = quantile(spike_prop[[pat]][[celltwo]], 0.99)
            lines(c(cell_index[cellone], cell_index[cellone]), c(ybottom+side, ybottom), 
                  lwd=2, col=sig_col(signif))
            
            lines(c(cell_index[cellone], cell_index[celltwo]), c(ybottom, ybottom), 
                  lwd=2, col=sig_col(signif))
            
            lines(c(cell_index[celltwo], cell_index[celltwo]), c(ybottom, ybottom+side), 
                  lwd=2, col=sig_col(signif))
            ybottom = ybottom - 0.05
          }
        }
      }
    }
  }
}
dev.off()

# matrices to store KDE peaks and 95% HDI quantiles
peaks = matrix(NA, nrow=length(ocell_types), ncol=length(pts))
colnames(peaks) = pts[order(ages)]
rownames(peaks) = ocell_types
lows = peaks
highs = peaks

for(pat in pts[order(ages)]){
  for(cell in ocell_types){
    if( !is.null(spike_prop[[pat]][[cell]]) ){
      if( !(pat=="R044"&&cell=="8+naive") && !(pat=="R013"&&cell=="Bmemory") && !(pat=="R013"&&cell=="Bnaive") ){
        hdi = hdi( spike_prop[[pat]][[cell]] )
        dens = density(spike_prop[[pat]][[cell]])
        hdi_dens = hdi( dens )
        
        peaks[cell, pat] = dens$x[which.max(dens$y)]
        lows[cell, pat] = hdi[1]
        highs[cell, pat] = hdi[2]
      }
    }
  }
}

# plot posterior densities 
pdf("./PDF/nearZeroProp_posterior.pdf", width=14, height=10)
{
  for(pat in pts[order(ages)]){
    par(mfrow=c(3,2), mar=c(6,6,6,3))
    for(cell in ocell_types){
      if( !is.null(spike_prop[[pat]][[cell]]) ){
        dens = density(spike_prop[[pat]][[cell]])
        hist(spike_prop[[pat]][[cell]], freq=FALSE, breaks=seq(0,1,0.005),
             col=myGrey(0.2), border=myGrey(0.2), 
             xlab="Near zero proportion", ylab="Density", 
             main=paste(pat, cell), cex.main=2, cex.lab=2, cex.axis=1.5)
        lines(dens, lwd=2)
        
        abline(v=c(lows[cell, pat], highs[cell, pat]), 
               lwd=2, col=myBlue(0.5), lty=2)
        
        abline(v=median_pi[cell,pat], 
               lwd=2, col=myBlack(0.5), lty=1)
      }
    }
  }
}
dev.off()



