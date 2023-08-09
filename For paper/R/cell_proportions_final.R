library(tidyverse)
library(ggsignif)
library(ggpubr)
library("ggsci")
options(scipen = 100) # prevent scientific notation

my_theme = theme_classic() + theme(text = element_text(size = 20))

pnums = read.csv(file="For paper/output/pnums.csv", stringsAsFactors = F)
names(pnums) = c("id", "pnums")

### load data ###

dat_l = read.csv(file = "For paper/input/proportions_data.csv", stringsAsFactors = FALSE)

##### fig1 (old s3a figure) #################################################

dat_s3a = subset(dat_l, (cell=="perc_T" | cell=="perc_B" |cell=="perc_NK" |cell=="perc_DR"))
head(dat_s3a)
table(dat_s3a$cell)

dat_s3a$cell = plyr::revalue(dat_s3a$cell, c(
  "perc_T" = "T cells",                 
  "perc_B" = "B cells",                 
  "perc_NK" = "NK cells",                
  "perc_DR" = "DR+ cells"              
))

dat_s3a$cell = factor(dat_s3a$cell, levels = c("T cells", "B cells", "NK cells", "DR+ cells"))

# numbers of patients and controls
head(dat_s3a)
dat_s3a[which(dat_s3a$cell=="T cells"),c("type", "cell", "value")]

s3a = ggplot(dat_s3a, aes(x=cell, y=value)) + my_theme + scale_colour_nejm(name=NULL) +
  geom_boxplot(outlier.shape=NA, aes(colour=type)) + 
  geom_point(position=position_dodge(width=0.75), alpha=0.6, aes(colour=type))+
  scale_y_continuous("% of PBMC", breaks=c(0,25,50,75,100), limits = c(0,100)) +
  xlab("")


# S3b - ratios of cells

dat_s3b = subset(dat_l, cell=="ratio_4_8" | cell=="ratio_naive_mat_cd4" |cell=="ratio_naive_mat_cd8" |cell=="ratio_naive_mat_B" | cell=="ratio_mat_naive_nk")
head(dat_s3b)
table(dat_s3b$cell)

dat_s3b$cell = plyr::revalue(dat_s3b$cell, c(
  "ratio_4_8"= "CD4+:CD8+\nT cells",
  "ratio_naive_mat_cd4" = "Naïve:Memory\nCD4+",
  "ratio_naive_mat_cd8" = "Naïve:Memory\nCD8+",
  "ratio_naive_mat_B" = "Naïve:Memory\nB cells", 
  "ratio_mat_naive_nk" = "Mature:Immature\nNK cells" 
))

dat_s3b$cell = factor(dat_s3b$cell, levels = c("CD4+:CD8+\nT cells", "Naïve:Memory\nCD4+", "Naïve:Memory\nCD8+", "Naïve:Memory\nB cells", "Mature:Immature\nNK cells"))


# split this up into two figures due to differences in scale

dat_s3bi = subset(dat_s3b, cell!="Mature:Immature\nNK cells" & cell!="Naïve:Memory\nB cells")
dat_s3bi$cell = factor(dat_s3bi$cell, levels=unique(dat_s3bi$cell))
dat_s3bii = subset(dat_s3b, cell=="Mature:Immature\nNK cells" | cell=="Naïve:Memory\nB cells" )
dat_s3bii$cell = factor(dat_s3bii$cell, levels=c("Naïve:Memory\nB cells","Mature:Immature\nNK cells"))

s3bii = ggplot(dat_s3bii, aes(x=cell, y=value)) + my_theme + scale_colour_nejm(name=NULL) +
  geom_boxplot(outlier.shape=NA, aes(colour=type)) + 
  geom_point(position=position_dodge(width=0.75), alpha=0.6, aes(colour=type))+
  scale_y_continuous("Ratio")  + #, breaks=c(0,25,50,75,100), limits = c(4.7,100)) +
  xlab("") +
  #stat_compare_means(aes(group = type), label = "p.signif", method="wilcox.test") 
  #stat_compare_means(aes(group = type), label = "p.format", method="wilcox.test") +
  geom_hline(yintercept=1, linetype="dashed") +
  theme(legend.position = "none")

s3bi = ggplot(dat_s3bi, aes(x=cell, y=value)) + my_theme + scale_colour_nejm(name=NULL) +
  geom_boxplot(outlier.shape=NA, aes(colour=type)) + 
  geom_point(position=position_dodge(width=0.75), alpha=0.6, aes(colour=type))+
  scale_y_continuous("Ratio")  + #, breaks=c(0,25,50,75,100), limits = c(4.7,100)) +
  xlab("") +
  #stat_compare_means(aes(group = type), label = "p.signif", method="wilcox.test") 
  #stat_compare_means(aes(group = type), label = "p.format", method="wilcox.test") +
  geom_hline(yintercept=1, linetype="dashed") +
  theme(legend.position = "none")

# S3c - % of DR+ cells

dat_s3c = subset(dat_l, cell=="perc_14" | cell=="perc_16" |cell=="perc_mDC" |cell=="perc_pDC" | cell=="perc_cd34")
head(dat_s3c)
table(dat_s3c$cell)

dat_s3c$cell = plyr::revalue(dat_s3c$cell, c(
  "perc_14" = "CD14+ Monocytes" ,
  "perc_16" = "CD16+ Monocytes",
  "perc_mDC" = "CD11c+ mDC",
  "perc_pDC" = "CD123+ pDC",
  "perc_cd34" = "CD34+ cells"
))

dat_s3c$cell = factor(dat_s3c$cell, levels = c("CD14+ Monocytes", "CD16+ Monocytes", "CD11c+ mDC", "CD123+ pDC", "CD34+ cells"))

# split this up into two figures due to differences in scale

dat_s3ci = subset(dat_s3c, cell=="CD14+ Monocytes" | cell=="CD16+ Monocytes")
dat_s3ci$cell = factor(dat_s3ci$cell, levels=unique(dat_s3ci$cell))
dat_s3cii = subset(dat_s3c, cell!="CD14+ Monocytes" & cell!="CD16+ Monocytes")
dat_s3cii$cell = factor(dat_s3cii$cell, levels=unique(dat_s3cii$cell))

s3ci = ggplot(dat_s3ci, aes(x=cell, y=value)) + my_theme + scale_colour_nejm(name=NULL) +
  geom_boxplot(outlier.shape=NA, aes(colour=type)) + 
  geom_point(position=position_dodge(width=0.75), alpha=0.6, aes(colour=type))+
  #  geom_jitter(width=0.1, height=0, aes(colour=type, group=type)) +
  scale_y_continuous("% of HLA DR+ Cells", breaks=c(0,25,50,75,100), limits = c(0,100)) +
  xlab("") +
  #stat_compare_means(aes(group = type), label = "p.signif", method="wilcox.test") 
  #stat_compare_means(aes(group = type), label = "p.format", method="wilcox.test") +
  theme(legend.position = "none")  

s3cii = ggplot(dat_s3cii, aes(x=cell, y=value)) + my_theme + scale_colour_nejm(name=NULL) +
  geom_boxplot(outlier.shape=NA, aes(colour=type)) + 
  geom_point(position=position_dodge(width=0.75), alpha=0.6, aes(colour=type))+
  #  geom_jitter(width=0.1, height=0, aes(colour=type, group=type)) +
  scale_y_continuous("% of HLA DR+ Cells", breaks=c(0,5,10,15,20), limits = c(0,15)) +
  xlab("") +
  #stat_compare_means(aes(group = type), label = "p.signif", method="wilcox.test") 
  #stat_compare_means(aes(group = type), label = "p.format", method="wilcox.test") +
  theme(legend.position = "none")  

##### stats ####

# make a dataframe to store results of lms

s3_stats = data.frame(group1 = "Control", group2="Patient", cell = c(levels(dat_s3a$cell), levels(dat_s3b$cell), levels(dat_s3c$cell)), beta=NA, se=NA, p=NA)

# T cells and HLA DR+ cells

mod = lm(value ~ type + age, subset(dat_s3a, cell=="T cells"))
summary(mod)
s3_stats[which(s3_stats$cell=="T cells"), c(4:6)] = summary(mod)$coefficients[2,c(1,2,4)]

mod = lm(value ~ type + age, subset(dat_s3a, cell=="HLA DR+ cells"))
summary(mod)
s3_stats[which(s3_stats$cell=="DR+ cells"), c(4:6)] = summary(mod)$coefficients[2,c(1,2,4)]

mod = lm(value ~ type + age, subset(dat_s3a, cell=="B cells"))
summary(mod)
s3_stats[which(s3_stats$cell=="B cells"), c(4:6)] = summary(mod)$coefficients[2,c(1,2,4)]

mod = lm(value ~ type + age, subset(dat_s3a, cell=="NK cells"))
summary(mod)
s3_stats[which(s3_stats$cell=="NK cells"), c(4:6)] = summary(mod)$coefficients[2,c(1,2,4)]

# age not significant and same results still signif even when age in model. How to present this?

## ratios

mod = lm(value ~ type + age, subset(dat_s3b, cell=="CD4+:CD8+\nT cells"))
summary(mod)
s3_stats[which(s3_stats$cell=="CD4+:CD8+\nT cells"), c(4:6)] = summary(mod)$coefficients[2,c(1,2,4)]

mod = lm(value ~ type + age, subset(dat_s3b, cell=="Naïve:Mature\nCD4+"))
summary(mod)
s3_stats[which(s3_stats$cell=="Naïve:Memory\nCD4+"), c(4:6)] = summary(mod)$coefficients[2,c(1,2,4)]

mod = lm(value ~ type + age, subset(dat_s3b, cell=="Naïve:Mature\nCD8+"))
summary(mod)
s3_stats[which(s3_stats$cell=="Naïve:Memory\nCD8+"), c(4:6)] = summary(mod)$coefficients[2,c(1,2,4)]

mod = lm(value ~ type, subset(dat_s3b, cell=="Naïve:Memory\nB cells"))
summary(mod)
s3_stats[which(s3_stats$cell=="Naïve:Memory\nB cells"), c(4:6)] = summary(mod)$coefficients[2,c(1,2,4)]

mod = lm(value ~ type + age, subset(dat_s3b, cell=="Mature:Immature\nNK cells"))
summary(mod)
s3_stats[which(s3_stats$cell=="Mature:Immature\nNK cells"), c(4:6)] = summary(mod)$coefficients[2,c(1,2,4)]


# proportions

mod = lm(value ~ type + age, subset(dat_s3c, cell=="CD14+ Monocytes"))
summary(mod)
s3_stats[which(s3_stats$cell=="CD14+ Monocytes"), c(4:6)] = summary(mod)$coefficients[2,c(1,2,4)]

mod = lm(value ~ type + age, subset(dat_s3c, cell=="CD16+ Monocytes"))
summary(mod)
s3_stats[which(s3_stats$cell=="CD16+ Monocytes"), c(4:6)] = summary(mod)$coefficients[2,c(1,2,4)]

mod = lm(value ~ type + age, subset(dat_s3c, cell=="CD11c+ mDC"))
summary(mod)
s3_stats[which(s3_stats$cell=="CD11c+ mDC"), c(4:6)] = summary(mod)$coefficients[2,c(1,2,4)]

mod = lm(value ~ type + age, subset(dat_s3c, cell=="CD123+ pDC"))
summary(mod)
s3_stats[which(s3_stats$cell=="CD123+ pDC"), c(4:6)] = summary(mod)$coefficients[2,c(1,2,4)]

mod = lm(value ~ type + age, subset(dat_s3c, cell=="CD34+ cells"))
summary(mod)
s3_stats[which(s3_stats$cell=="CD34+ cells"), c(4:6)] = summary(mod)$coefficients[2,c(1,2,4)]

### add these p values and betas to graphs.....

s3_stats
# FDR corrected ps
s3_stats$p.adj = round(p.adjust(s3_stats$p, method="fdr"),4)
s3_stats$p = round(s3_stats$p, 4)
s3_stats$p.label = paste("FDR = ",s3_stats$p.adj, sep="")
s3_stats$p.label_2 = ifelse(s3_stats$p.adj > 0.05, "FDR = ns", s3_stats$p.label)
s3_stats$p.label_unadj = ifelse(s3_stats$p==0, "p < 0.001", ifelse(
  s3_stats$p<0.05, paste("p = ",s3_stats$p, sep=""), "p = ns"))
#s3_stats$p.label_unadj_2 = ifelse(s3_stats$p > 0.05, "ns", s3_stats$p.label_unadj)
s3_stats$p.label_all = paste(s3_stats$p.label_unadj, s3_stats$p.label_2, sep="\n")
s3_stats$y.position = c(rep(105,4), rep(8.2,3), rep(54,2), rep(105,2), rep(15.75,3))
s3_stats$y.position.all = 0.9*s3_stats$y.position

# both p labels
s3a = s3a + stat_pvalue_manual(subset(s3_stats, cell%in%levels(dat_s3a$cell)), label="p.label_all", tip.length = 0.01, x="cell", y.position = "y.position.all")
s3bi = s3bi + stat_pvalue_manual(subset(s3_stats, cell%in%levels(dat_s3bi$cell)), label="p.label_all", tip.length = 0.01, x="cell", y.position = "y.position.all")
s3bii = s3bii + stat_pvalue_manual(subset(s3_stats, cell%in%levels(dat_s3bii$cell)), label="p.label_all", tip.length = 0.01, x="cell", y.position = "y.position.all")
s3ci = s3ci + stat_pvalue_manual(subset(s3_stats, cell%in%levels(dat_s3ci$cell)), label="p.label_all", tip.length = 0.01, x="cell", y.position = "y.position.all")
s3cii = s3cii + stat_pvalue_manual(subset(s3_stats, cell%in%levels(dat_s3cii$cell)), label="p.label_all", tip.length = 0.01, x="cell", y.position = "y.position.all")


# write out plots to files

#tiff("output/sarah/proportions_t_b_nk_dr.tiff", width=800, height=400)
#s3a
#dev.off()

#tiff("output/sarah/ratios1.tiff", width=400, height=400)
#s3bi
#dev.off()

#tiff("output/sarah/ratios2.tiff", width=400, height=400)
#s3bii
#dev.off()

#tiff("output/sarah/proportions_of_dr_cells1.tiff", width=400, height=400)
#s3ci
#dev.off()

#tiff("output/sarah/proportions_of_dr_cells2.tiff", width=400, height=400)
#s3cii
#dev.off()

