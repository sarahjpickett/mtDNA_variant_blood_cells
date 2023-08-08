library(ggplot2)
library(plyr)
library(tidyverse)
library(ggpubr)
library(ggsci)

my_theme = theme_classic() + theme(text = element_text(size = 20))


all_all = read.csv(file="output/all_cleaned_sc_data.csv", stringsAsFactors = F)
all_all = subset(all_all, cell_type!="8+ EM/TEMRA")


all_all$major_cell_type = substr(all_all$cell_type, 1, 3)
all_all$major_cell_type = ifelse(substr(all_all$major_cell_type, 1, 1)=="B", "B", all_all$major_cell_type)

summaries = ddply(all_all, .(Patient, cell_type), summarise, median = median(HET, na.rm=T), mean = mean(HET, na.rm=T))

all_all = merge(all_all, summaries, all=T)
all_all$cell_type = factor(all_all$cell_type, levels=c("34+ Progenitor", "Monocyte", "4+ naive", "4+ CM", "4+ EM","4+ TEMRA",
                                                       "8+ naive", "8+ CM", "8+ EM", "8+ TEMRA", "B naive", "B memory"))

pnums = read.csv(file="pnums.csv", stringsAsFactors = F)
names(pnums)[1] = "Patient"
all_all = merge(all_all, pnums, all.x=T)

## add in a dummy variable for P17 cd34+ cells and other missing cell types

temp = data.frame("Patient"="R027", "pnums"="P17", "cell_type"=c("34+ Progenitor", "4+ CM"), "major_cell_type"=c("34+", "4+"))
all_all= merge(all_all, temp, all=T)
temp = data.frame("Patient"="R025", "pnums"="P15", "cell_type"="8+ naive", "major_cell_type"="8+")
all_all= merge(all_all, temp, all=T)
temp = data.frame("Patient"=c("R149", "R044", "R013", "R025", "R027"), "pnums"=c("P22", "P19", "P10", "P15", "P17"), "cell_type"="4+ TEMRA", "major_cell_type"="4+")
all_all= merge(all_all, temp, all=T)
temp = data.frame("Patient"="R028", "pnums"="P18", "cell_type"="4+ CM", "major_cell_type"="4+")
all_all= merge(all_all, temp, all=T)

# order of major_cell_type

unique(all_all$major_cell_type)
all_all$major_cell_type = factor(all_all$major_cell_type, levels=c("34+", "Mon", "4+ ", "8+ ", "B"))

all_all = all_all[order(all_all$cell_type),]

## just dotplot ####
#pdf("output/manuscript/combined_dotplot3.pdf", width=12, height=8)

for(patient in unique(all_all$pnums)){
  tiff(file=paste("output/manuscript/dotplot_", patient, ".tiff", sep=""), width=550, height=480)
  
  plot = ggplot(subset(all_all, pnums==patient), aes(x=cell_type, y=HET)) + my_theme +
    #geom_dotplot(binaxis="y", stackdir="center", dotsize = 1.5,  method="histodot", binwidth=0.5, stackratio = 0.5, alpha=0.8, aes(fill=major_cell_type, colour=major_cell_type))+
    geom_dotplot(binaxis="y", stackdir="center", dotsize = 1.2,  method="histodot", binwidth=1, stackratio = 0.15, alpha=0.6, aes(fill=major_cell_type, colour=major_cell_type))+
    #geom_dotplot(binaxis="y", stackdir="center", dotsize = 1,  method="histodot", binwidth=1, stackratio = 0.1, alpha=0.4, aes(fill=major_cell_type, colour=major_cell_type))+
    ggtitle(paste(patient, ": Age = ", unique(subset(all_all, pnums==patient & !is.na(age))$age), sep=""))+
    coord_cartesian(ylim=c(4.5,100)) +
    ylab("m.3243A>G level (%)") +
    xlab("") +
    scale_color_nejm(name = "", na.translate = T, drop = FALSE) +
    scale_fill_nejm(name = "", na.translate = T, drop = FALSE) +
    ylim(c(0,100)) +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    theme(legend.position = "none") +
    geom_vline(xintercept = c(2,3,7,11) - 0.5, col='grey', lwd=0.5, linetype="dotted")
  print(plot)
  dev.off()
}