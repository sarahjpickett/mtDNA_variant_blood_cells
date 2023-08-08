library(tidyverse)
library(ggsignif)
library(ggpubr)
library("ggsci")
options(scipen = 100) # prevent scientific notation

my_theme = theme_classic() + theme(text = element_text(size = 20))

pnums = read.csv(file="output/sarah/pnums.csv", stringsAsFactors = F)
names(pnums) = c("id", "pnums")

### load data and tidy ####

props = read.csv(file="input/200505 Heteroplasmy Paper All Data_IGF_final_051022.csv", stringsAsFactors = F)
head(props)

# remove controls that are not part of the study (OR x2 and me!) and R005 (sample with not enough blood and no het data available)
props = subset(props, id!="Sample 1" & id!="Sample 2" & id!="Control OR" & id!="R005")

# also remove sample Control 8 and A2546 due to strange flow data
props = subset(props, id!="Control 8" & id!="A2546")

# remove non-3243 patients "R0286"           "R0179"           "R0181"           "R0308" 
props = subset(props, id!="R0286" & id!="R0179" & id!="R0181" & id!="R0308")

# remove R017 (no data)  
props = subset(props, id!="R017")

# calculate a few more ratios...

names(props)

# some sense checks to make sure I understand data. 
sum(round(props$raw_cd4 / props$raw_cd8, 3) != round(props$X4_8_ratio, 3), na.rm = T)

sum(round(props$raw_4_naive / props$raw_cd4 * 100, 3) != round(props$perc_4_naive_of_4,3), na.rm = T)
props[round(props$raw_4_naive / props$raw_cd4 * 100, 2) != round(props$perc_4_naive_of_4,2),]


sum(round((props$raw_4_em +  props$raw_4_cm + props$raw_4_temra) / props$raw_cd4 * 100, 3) != round(props$perc_4_memory_of_4,3), na.rm = T)

sum(round((props$raw_8_em +  props$raw_8_cm + props$raw_8_temra) / props$raw_cd8 * 100, 3) != round(props$perc_8_memory_of_8,3), na.rm = T)

sum(round(props$raw_B1 / props$raw_B * 100, 3) != round(props$perc_B1, 3), na.rm=T) # why doesn't this work out?


# just select raw data columns then recalculate percentages etc....

dat = props[,c(1:36)]

dat=props
names(dat)

# calculate percentages and ratios

dat$cd4_naive_memory_ratio = dat$raw_4_naive / (dat$raw_4_em +  dat$raw_4_cm + dat$raw_4_temra)

dat$cd8_naive_memory_ratio = dat$raw_8_naive / (dat$raw_8_em +  dat$raw_8_cm + dat$raw_8_temra)

dat$perc_T = (dat$raw_T / (dat$raw_T + dat$raw_B + dat$raw_NK + dat$raw_DR)) * 100 # ok
round(dat$perc_T,4) == round(props$perc_T,4)

dat$perc_B = (dat$raw_B / (dat$raw_T + dat$raw_B + dat$raw_NK + dat$raw_DR)) * 100 # ok
round(dat$perc_B,4) == round(props$perc_B,4)

dat$perc_NK = (dat$raw_NK / (dat$raw_T + dat$raw_B + dat$raw_NK + dat$raw_DR)) * 100 # ok
round(dat$perc_NK,4) == round(props$perc_NK,4) # this doesn't match because there is a mistake in the spreadsheet when calculating total number of cells. 
# I have not corrected this but will use these corrected data in downstream analyses


dat$perc_DR = (dat$raw_DR / (dat$raw_T + dat$raw_B + dat$raw_NK + dat$raw_DR)) * 100 # ok
round(dat$perc_DR,4) == round(props$perc_DR,4)


dat$perc_cd4 = dat$raw_cd4 / dat$raw_T * 100
dat$perc_cd8 = dat$raw_cd8 / dat$raw_T * 100

dat$ratio_4_8 = dat$raw_cd4 / dat$raw_cd8
round(dat$ratio_4_8,4) == round(props$X4_8_ratio,4) # check

dat$ratio_naive_mat_cd4 = dat$raw_4_naive / (dat$raw_4_cm + dat$raw_4_em + dat$raw_4_temra)
dat$ratio_naive_mat_cd8 = dat$raw_8_naive / (dat$raw_8_cm + dat$raw_8_em + dat$raw_8_temra)
dat$ratio_naive_mat_nk = dat$raw_56_nk / dat$raw_16_nk # 16 high are mature cells
dat$ratio_mat_naive_nk = dat$raw_16_nk / dat$raw_56_nk # 16 high are mature cells
dat$ratio_naive_mat_B = dat$raw_B2 / dat$raw_B4

# as a % of HLA DR+ cells
# monocytes
dat$perc_14 = dat$raw_14 / dat$raw_DR *100 
dat$perc_16 = dat$raw_16 / dat$raw_DR *100
# DCs
dat$perc_mDC = dat$raw_mDC / dat$raw_DR *100 
dat$perc_pDC = dat$raw_pDC / dat$raw_DR *100
# CD34+
dat$perc_cd34 = dat$raw_cd34 / dat$raw_DR *100

# t cells subsets as a % of t cell subsets....
# quick check....
dat$raw_4_em +  dat$raw_4_cm + dat$raw_4_temra + dat$raw_4_naive -dat$raw_cd4 

dat$perc_naive_cd4 = dat$raw_4_naive / dat$raw_cd4 * 100
dat$perc_mature_cd4 = 100 - dat$perc_naive_cd4
dat$perc_sorted_mature_cd4 = (dat$raw_4_temra+dat$raw_4_em) / dat$raw_cd4 * 100

dat$perc_naive_cd8 = dat$raw_8_naive / dat$raw_cd8 * 100
dat$perc_mature_cd8 = 100 - dat$perc_naive_cd8
dat$perc_sorted_mature_cd8 = (dat$raw_8_temra+dat$raw_8_em) / dat$raw_cd8 * 100

# and B cells

dat$perc_naive_B = dat$raw_B2 / dat$raw_B *100
dat$perc_mature_B = dat$raw_B4 / dat$raw_B *100

# and NK cells

dat$perc_16_NK = dat$raw_16_nk / dat$raw_NK *100
dat$perc_56_NK = dat$raw_56_nk / dat$raw_NK *100

# change to long format

to_gather = names(dat)
to_gather = to_gather[-c(1:8)]

dat_l = gather(dat, cell, value, to_gather)
dim(dat_l)
head(dat_l)


## add pnums to dat
dat = merge(dat, pnums, all.x=T)
dat[,c(1,91)]

### draw graphs ans perform analysis ####


##### PLUMB1 (old s3a figure) #################################################

#dat_s3a = subset(dat_l, (cell=="perc_T" | cell=="perc_B" |cell=="perc_NK" |cell=="perc_DR") & id!="R002" & id!="R009") # this also removes the two patients with lymphocyte counts outside of the normal range
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



### Ratios vs age 

head(dat_s3b)
ggplot(subset(dat_s3b, cell=="CD4+:CD8+\nT cells"), aes(x=age, y=value)) + my_theme +
  geom_point(aes(colour=type)) +
  geom_smooth(aes(colour=type, fill=type), method="lm") +
  ylab("CD4+:CD8+ Ratio") +
  geom_abline(intercept=1, slope=0, linetype="dotted")

summary(lm(value ~ age+type, subset(dat_s3b, cell=="CD4+:CD8+\nT cells")))
summary(lm(value ~ age*type, subset(dat_s3b, cell=="CD4+:CD8+\nT cells")))
summary(lm(value ~ age, subset(dat_s3b, cell=="CD4+:CD8+\nT cells" & type=="Control")))
summary(lm(value ~ age, subset(dat_s3b, cell=="CD4+:CD8+\nT cells" & type=="Patient")))
# CHECK!!!

## get het and cn per cell type and merge in

head(het_clean)

het_cell_type = spread(het_clean[,c("patient_id", "cell_type", "het")], cell_type, het)
names(het_cell_type)[1] = "id" 
names(het_cell_type)[2:length(names(het_cell_type))] = paste("het_", names(het_cell_type)[2:length(names(het_cell_type))], sep="")
het_cell_type$id = as.character(het_cell_type$id)


cn_cell_type = spread(het_clean[,c("patient_id", "cell_type", "normalised_cn")], cell_type, normalised_cn)
names(cn_cell_type)[1] = "id" 
names(cn_cell_type)[2:length(names(cn_cell_type))] = paste("cn_", names(cn_cell_type)[2:length(names(cn_cell_type))], sep="")
cn_cell_type$id = as.character(cn_cell_type$id)

head(dat)

dat_original = dat
dat = merge(dat_original, cn_cell_type, all.x=T)
dat = merge(dat, het_cell_type, all.x=T)

head(dat)

# now graph cell-specific het against cell proportions - COME BACK TO THIS AFTER HOLIDAY!!!

ggplot(dat, aes(x=`het_CD4+ Naive T`, y=perc_T)) + my_theme +
  geom_point() +
  geom_smooth(method="lm")

ggplot(dat, aes(x=`het_CD8+ Naive T`, y=perc_T)) + my_theme +
  geom_point() +
  geom_smooth(method="lm")

ggplot(dat, aes(x=`het_CD8+ Naive T`, y=ratio_naive_mat_cd8)) + my_theme +
  geom_point() +
  geom_smooth(method="lm") +
  stat_cor(aes(label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                                                                       ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
                                                                       sep = "*', '*")), label.x = 0, label.y = 8, key_glyph = draw_key_rect)


ggplot(dat, aes(x=`het_CD8+ Memory T`, y=ratio_naive_mat_cd8)) + my_theme +
  geom_point() +
  geom_smooth(method="lm") +
  stat_cor(aes(label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                             ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
                             sep = "*', '*")), label.x = 0, label.y = 8, key_glyph = draw_key_rect)


ggplot(dat, aes(x=`het_CD4+ Naive T`, y=ratio_naive_mat_cd4)) + my_theme +
  geom_point() +
  geom_smooth(method="lm")+
  stat_cor(aes(label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                             ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
                             sep = "*', '*")), label.x = 0, label.y = 4, key_glyph = draw_key_rect)


ggplot(dat, aes(x=`het_CD4+ Memory T`, y=ratio_naive_mat_cd4)) + my_theme +
  geom_point() +
  geom_smooth(method="lm")+
  stat_cor(aes(label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                             ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
                             sep = "*', '*")), label.x = 0, label.y = 4, key_glyph = draw_key_rect)


ggplot(dat, aes(x=`het_Naive B`, y=ratio_naive_mat_B)) + my_theme +
  geom_point() +
  geom_smooth(method="lm")+
  stat_cor(aes(label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                             ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
                             sep = "*', '*")), label.x = 0, label.y = 20, key_glyph = draw_key_rect)


ggplot(dat, aes(x=`het_Memory B`, y=ratio_naive_mat_B)) + my_theme +
  geom_point() +
  geom_smooth(method="lm")+
  stat_cor(aes(label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                             ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
                             sep = "*', '*")), label.x = 0, label.y = 20, key_glyph = draw_key_rect)


# is this just an age effect?
ggplot(dat, aes(x=age, y=ratio_naive_mat_cd4)) + my_theme +
  geom_point() +
  geom_smooth(method="lm")+
  stat_cor(aes(label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                             ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
                             sep = "*', '*")), label.x = 0, label.y = 20, key_glyph = draw_key_rect)

summary(lm(ratio_naive_mat_cd4 ~ age, subset(dat, type=="Control")))
summary(lm(ratio_naive_mat_cd4 ~ age, subset(dat, type=="Patient")))
summary(lm(ratio_naive_mat_cd4 ~ age, dat))
summary(lm(ratio_naive_mat_cd4 ~ age + `het_CD4+ Naive T`, dat))
summary(lm(ratio_naive_mat_cd4 ~ age*type - type, dat))

ggplot(dat, aes(x=type, y=age)) + my_theme +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height=0, width=0.2)
wilcox.test(age~type, dat)

cd4_ratio_age = ggplot(dat, aes(x=age, y=ratio_naive_mat_cd4)) + my_theme +
  geom_point(aes(colour=type)) +
  geom_smooth(method="lm", aes(colour=type, fill=type))+
  stat_cor(aes(colour=type, label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                             ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
                             sep = "*', '*")), label.x = 58, label.y = c(4,3.8), key_glyph = draw_key_rect)+
  xlab("Age") + ylab("Naive:Memory CD4+ Ratio") + 
  scale_color_npg(name = "", na.translate = F) +
  scale_fill_npg(name = "", na.translate = F) +
  theme(legend.position = "none")

cd8_ratio_age = ggplot(dat, aes(x=age, y=ratio_naive_mat_cd8)) + my_theme +
  geom_point(aes(colour=type)) +
  geom_smooth(method="lm", aes(colour=type, fill=type))+
  stat_cor(aes(colour=type, label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                                          ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
                                          sep = "*', '*")), label.x = 58, label.y = c(7.6, 7.2), key_glyph = draw_key_rect) +
  xlab("Age") + ylab("Naive:Memory CD8+ Ratio") + 
  scale_color_npg(name = "", na.translate = F) +
  scale_fill_npg(name = "", na.translate = F) +
  theme(legend.position = "none")

b_ratio_age = ggplot(dat, aes(x=age, y=ratio_naive_mat_B)) + my_theme +
  geom_point(aes(colour=type)) +
  geom_smooth(method="lm", aes(colour=type, fill=type))+
  stat_cor(aes(colour=type, label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                                          ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
                                          sep = "*', '*")), label.x = 58, label.y = c(24.5,23.25), key_glyph = draw_key_rect)+
  xlab("Age") + ylab("Naive:Memory B Cell Ratio") + 
  scale_color_npg(name = "", na.translate = F) +
  scale_fill_npg(name = "", na.translate = F) +
  theme(legend.position = "none")


nk_ratio_age = ggplot(dat, aes(x=age, y=ratio_mat_naive_nk)) + my_theme +
  geom_point(aes(colour=type)) +
  geom_smooth(method="lm", aes(colour=type, fill=type))+
  stat_cor(aes(colour=type, label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                                          ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
                                          sep = "*', '*")), label.x = 58, label.y = c(50,47.5), key_glyph = draw_key_rect) +
  xlab("Age") + ylab("Mature:Immature NK Cell Ratio") + 
  scale_color_npg(name = "", na.translate = F) +
  scale_fill_npg(name = "", na.translate = F) +
  theme(legend.position = "none")


# Extract legends so they can be printed separately
leg_portrait = get_legend(ggplot(dat, aes(x=age, y=age)) + my_theme +
  geom_smooth(aes(colour=type, fill=type)) + 
  scale_color_npg(name = "", na.translate = F) +
  scale_fill_npg(name = "", na.translate = F))
as_ggplot(leg_portrait)

leg_landscape = get_legend(ggplot(dat, aes(x=age, y=age)) + my_theme +
                            geom_smooth(aes(colour=type, fill=type)) + 
                            scale_color_npg(name = "", na.translate = F) +
                            scale_fill_npg(name = "", na.translate = F)+
                            theme(legend.position = "top"))
as_ggplot(leg_landscape)

#tiff("output/sarah/control_patient_legend_portrait.tiff", width=400, height=400)
#as_ggplot(leg_portrait)
#dev.off()

#tiff("output/sarah/control_patient_legend_landscape.tiff", width=400, height=400)
#as_ggplot(leg_landscape)
#dev.off()

#tiff("output/sarah/cd4_ratio_age.tiff", width=400, height=400)
#cd4_ratio_age
#dev.off()
#tiff("output/sarah/cd8_ratio_age.tiff", width=400, height=400)
#cd8_ratio_age
#dev.off()
#tiff("output/sarah/b_ratio_age.tiff", width=400, height=400)
#b_ratio_age
#dev.off()
#tiff("output/sarah/nk_ratio_age.tiff", width=400, height=400)
#nk_ratio_age
#dev.off()
