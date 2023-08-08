
library(forcats)
library(nlme)
library(dplyr)
library(ggpubr)
library("ggsci")
library("gridExtra")
options(scipen = 100) # prevent scientific notation

my_theme = theme_classic() + theme(text = element_text(size = 20))


#load scripts needed to run this

source("For paper/R/heteroplasmy_data_loading.R")


# write out pnums

write.csv(file="For paper/output/pnums.csv", quote=F, row.names=F, unique(all[which(all$type=="Patient"), c("patient_id", "pnums")]))

### FIGURE 1 #####

# het vs. age for cell_types

het_clean$maturity2 = ifelse(het_clean$maturity=="CD34+", NA, 
                             ifelse(het_clean$maturity=="", NA, as.character(het_clean$maturity)))
het_clean$maturity2 = factor(het_clean$maturity2, levels=c("Naive/Immature", "Memory/Mature"))

temp1 = unique(het_clean$broad_cell_type)
temp2 = c("CD34+", "CD4+ T", "CD8+ T", "Monocytes", "B Cells", "Natural Killer", "Dendritic",     
"Platelets", "Granulocytes", "Whole Blood") 
temp3 = data.frame(broad_cell_type = temp1, verbose_broad_cell_type = temp2)

het_clean = merge(het_clean, temp3, all.x=T)

het_clean$verbose_broad_cell_type = factor(het_clean$verbose_broad_cell_type, levels=c("CD34+", "B Cells", "CD4+ T", "CD8+ T", "Natural Killer", "Monocytes", "Dendritic",     
                                                                          "Granulocytes", "Platelets", "Whole Blood"))

cell_het_age = ggplot(het_clean, aes(x=age, y=het)) + my_theme +
  xlab("Age") + ylab("m.3243A>G Level") +
  geom_point(aes(colour=maturity2)) + 
  geom_point(data=subset(het_clean, is.na(maturity2))) +
  geom_smooth(method="lm", aes(colour=maturity2, fill=maturity2), alpha=0.2) +
  geom_smooth(method="lm", data=subset(het_clean, is.na(maturity2)), alpha=0.2) +
  facet_wrap(~verbose_broad_cell_type, nrow=2) +
  scale_color_npg(name = "", na.translate = F) +
  scale_fill_npg(name = "", na.translate = F) +
  #stat_cor(size=4, method = "pearson", label.x =40, label.y = 60, p.accuracy = 0.001, r.accuracy = 0.01)
  stat_cor( aes(colour = maturity2, label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                                                ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
                                                sep = "*', '*")), label.x = c(30, 30), label.y = c(65, 50), key_glyph = draw_key_rect) +
  stat_cor( data=subset(het_clean, is.na(maturity2)),aes(label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                                               ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
                                               sep = "*', '*")), label.x = 30, label.y = 65, key_glyph = draw_key_rect) +
  theme(legend.position="right")
  
#make all points grey and highlight 6 patients 

scpatients = subset(het_clean, patient_id == "R027" | patient_id == "R025" | patient_id == "R044" | patient_id == "R0149" | patient_id == "R013" | patient_id == "R028" )


#remove whole blood from this graph and age coloured boxplot

minusWB = subset(het_clean, cell_type != "Whole Blood")

scpatients1 = subset(scpatients, cell_type != "Whole Blood" )

highlight = ggplot(minusWB, aes(x=cell_type, y=het, group = pnums)) + my_theme +
  geom_point(aes(), colour = "gray") + geom_line(aes(), colour = "gray") +
  theme(axis.text.x = element_text(angle = 45, hjust=1,size = 11)) +
  ylab("m.3243A>G level") +
  xlab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_line(data = scpatients1, aes(group = pnums, colour = pnums), size = 1) +
  geom_point(data = scpatients1, aes(group = pnums, colour = pnums)) +
  theme(legend.title = element_blank()) +
  scale_color_npg(name = "Patient")+
  ylim(0, 75)



final_het = ggplot(minusWB, aes(x=cell_type, y=het)) + my_theme +
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(colour=age), height=0, width=0.1) +
  theme(axis.text.x = element_text(angle = 45, hjust=1,size = 15)) +
  ylab("m.3243A>G level") +
  xlab("") +
  scale_colour_gradient(high="red", low="royalblue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position="none") +
  ylim(0, 75)


#make separate dot plot of just whole blood values 

WB = subset(het_clean, cell_type == "Whole Blood")

#wb_scatter = ggplot(WB, aes(x = age, y = het)) + theme_classic() + 
#  geom_point(size=2) +
#  theme(axis.text.x = element_text(angle = 45, hjust=1,size = 11)) +
#  ylab("m.3243A>G level") +
#  xlab("Age") +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


WBPLOT = ggplot(WB, aes(x=cell_type, y=het)) + my_theme +
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(colour=age), height=0, width=0.1) +
  theme(axis.text.x = element_text(angle = 45, hjust=1,size = 15)) +
  ylab("") +
  xlab("") +
  scale_colour_gradient(high="red", low="royalblue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylim(0, 75)


# combine plots into a single figure

final_het_wb = ggarrange(final_het, WBPLOT, align="h", widths=c(4,1.5))
highlight2 = ggarrange(highlight, NA, align="h", widths=c(4.41,0.54))

#pdf("output/sarah/figure2.pdf")
#ggarrange(cell_het_age,
#  ggarrange(final_het_wb, highlight2, labels=c("b)","c)"), nrow=2),
#labels=c("a)","", ""), nrow=1)
#dev.off()

#patchwork

# write out plots and create multi-panel figure in ppt

tiff(filename="For paper/output/final_het_wb.tiff", width=550, height=480)
final_het_wb
dev.off()

tiff(filename="For paper/output/highlight2.tiff", width=550, height=480)
#highlight2
highlight
dev.off()

#tiff(filename="output/sarah/cell_het_age.tiff", height=960, width=480)
tiff(filename="For paper/output/cell_het_age_wide.tiff", width=1100, height=480)
cell_het_age
dev.off()



# now check stats for figure legend and this section of the paper.


library(nlme)

head(minusWB)

mod = lme(het ~ cell_type, random = ~ 1|pnums, data=subset(minusWB, !is.na(het)))
summary(mod)
sink(file="For paper/output/bulk_lmm.txt")
summary(mod)
sink()

# change order of factors so that naive is reference category
minusWB$maturity2 = factor(minusWB$maturity2, levels = c("Naive/Immature", "Memory/Mature"))
sink(file="For paper/output/bulk_lmm_mature_immature.txt")
mod4 = lme(het ~ age*maturity2, random = ~ 1|pnums, data=subset(minusWB, verbose_broad_cell_type=="CD4+ T"  & !is.na(het)))
summary(mod4)
mod8 = lme(het ~ age*maturity2, random = ~ 1|pnums, data=subset(minusWB, verbose_broad_cell_type=="CD8+ T"  & !is.na(het)))
summary(mod8)
modNK = lme(het ~ age*maturity2, random = ~ 1|pnums, data=subset(minusWB, verbose_broad_cell_type=="Natural Killer"  & !is.na(het)))
summary(modNK)
modB = lme(het ~ age*maturity2, random = ~ 1|pnums, data=subset(minusWB, verbose_broad_cell_type=="B Cells"  & !is.na(het)))
summary(modB)
sink()



### FIGURE 2 CN data ####

all$log_cn = log10(all$normalised_cn)

all$maturity2 = ifelse(all$maturity=="CD34+", NA, 
                             ifelse(all$maturity=="", NA, as.character(all$maturity)))
all$maturity2 = factor(all$maturity2, levels=c("Naive/Immature", "Memory/Mature"))

temp1 = unique(all$cell_type)
temp2 = c("Dendritic", "Dendritic",  "Monocytes","Monocytes", "CD34+", "CD4+ T", "CD4+ T", "Natural Killer", "Natural Killer", 
          "CD8+ T","CD8+ T","B Cells", "B Cells")
temp3 = data.frame(cell_type = temp1, verbose_broad_cell_type = temp2)

all = merge(all, temp3, all.x=T)

all$verbose_broad_cell_type = factor(all$verbose_broad_cell_type, levels=c("CD34+", "B Cells", "CD4+ T", "CD8+ T", "Natural Killer", "Monocytes", "Dendritic",     
                                                                                       "Granulocytes", "Platelets", "Whole Blood"))



logged_cn = 
  ggplot(all, aes(x=cell_type, y=log_cn)) + my_theme +
  geom_boxplot(outlier.shape=NA, aes(colour=type, fill = type), alpha=0.1)+
  geom_point(position=position_jitterdodge(jitter.height=0, jitter.width=0.1),alpha=0.8, size =1.8, aes(colour=type, group=type)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab("Log10 mtDNA copy number") +
  xlab("") +
#  scale_colour_manual(values=c("grey21", "steelblue4")) +
#  scale_fill_manual(values=c("gray86", "lightblue")) +
  scale_color_npg(name = "", na.translate = F) +
  scale_fill_npg(name = "", na.translate = F) +
  theme(legend.title = element_blank())  

# png("output/logged_cn.png", width=800, height=550)
#logged_cn
# dev.off()

## response to reviewer's comments. Investigate CN with age / 3243 level

ggplot(all, aes(x=cell_type, y=log_cn)) + my_theme +
  geom_boxplot(outlier.shape=NA, aes(fill = type), alpha=0.1)+
  geom_point(position=position_jitterdodge(jitter.height=0, jitter.width=0.1),alpha=0.8, size =1.8, aes(colour=het, group=type)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab("Log10 mtDNA copy number") +
  xlab("") +
  theme(legend.title = element_blank())  





wilcox_results = data.frame(cell_type=unique(all$cell_type), W = NA, p = NA)

for(cell in unique(all$cell_type)) {
  check = subset(all, cell_type==cell)
  w = wilcox.test(log_cn ~ type, data=check, na.rm = TRUE)
  wilcox_results[which(wilcox_results$cell_type==cell), "p"] = w$p.value
  wilcox_results[which(wilcox_results$cell_type==cell), "W"] = w$statistic
}

wilcox_results

wilcox_results$fdr = round(p.adjust(wilcox_results$p, method = "fdr", n = length(wilcox_results$p)), 4)
wilcox_results$p = round(wilcox_results$p, 4)

write.csv(wilcox_results, file="For paper/output/cn_log_wilcox_results.csv", quote=F, row.names=F)

  
  ggplot(all, aes(x=age, y=log_cn)) + my_theme +
    geom_point(aes(colour=type)) +
    geom_smooth(method="lm", aes(colour=type, fill=type), alpha=0.2) +
    scale_y_continuous("Normalised log10 mtDNA CN") +
    scale_x_continuous("Age") +
    facet_wrap(~~cell_type, ncol=4) +
    stat_cor(size=4, aes(colour = type, label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                                                   ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
                                                   sep = "*', '*")), label.x = c(50, 50), label.y = c(3.5, 3), key_glyph = draw_key_rect) +
    scale_color_npg(name = "", na.translate = F) +
    scale_fill_npg(name = "", na.translate = F) 
    #scale_colour_manual(values=c("grey21", "steelblue4")) +
    #scale_fill_manual(values=c("gray86", "lightblue")) 
    #scale_colour_manual(values=c("#6F99ADFF", "#E18727FF")) +
    #scale_fill_manual(values=c("#6F99ADFF", "#E18727FF")) 
    
# now not split by cell type
  
  
  ggplot(all, aes(x=age, y=log_cn)) + my_theme +
    geom_point(aes(colour=type)) +
    geom_smooth(method="lm", aes(colour=type, fill=type), alpha=0.2) +
    scale_y_continuous("Normalised log10 mtDNA CN") +
    scale_x_continuous("Age") +
    stat_cor(size=4, aes(colour = type, label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                                              ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
                                              sep = "*', '*")), label.x = c(50, 50), label.y = c(3.5, 3), key_glyph = draw_key_rect) +
    scale_color_npg(name = "", na.translate = F) +
    scale_fill_npg(name = "", na.translate = F) 

# need to run a linear mixed model
  
sink(file="For paper/output/cn_age_lme_results.txt")  
mod_p = lme(log_cn ~ age, random = ~ 1|pnums, data=subset(all,!is.na(log_cn) & type=="Patient"))
summary(mod_p)  
  
mod_c = lme(log_cn ~ age, random = ~ 1|patient_id, data=subset(all,!is.na(log_cn) & type=="Control"))
summary(mod_c) 
sink()
  
# extract these values and add regression lines to graph!

cn_age = ggplot(all, aes(x=age, y=log_cn)) + my_theme +
  geom_point(aes(colour=type)) +
  #geom_smooth(method="lm", aes(colour=type, fill=type), alpha=0.2) +
  scale_y_continuous("log10 mtDNA copy number") +
  scale_x_continuous("Age") +
  scale_color_npg(name = "", na.translate = F) +
  scale_fill_npg(name = "", na.translate = F) +
  geom_abline(intercept = mod_p$coefficients$fixed[1], slope = mod_p$coefficients$fixed[2], colour="#0072B5FF", size=1, alpha=0.8) +
  geom_abline(intercept = mod_c$coefficients$fixed[1], slope = mod_c$coefficients$fixed[2], colour="#BC3C29FF", size=1, alpha=0.8) +
  theme(legend.position = "right")
  
sink(file="For paper/output/cn_het_lme_results.txt")  
mod_p_het = lme(log_cn ~ het, random = ~ 1|pnums, data=subset(all,!is.na(log_cn) & type=="Patient"))
summary(mod_p_het)
sink()

cn_het = ggplot(all, aes(x=het, y=log_cn)) + my_theme +
  geom_point(aes(colour=age)) +
  scale_colour_gradient(high="red", low="royalblue", name="Age") +
  theme(legend.position = "right") +
  #geom_smooth(method="lm", aes(colour=type, fill=type), alpha=0.2) +
  scale_y_continuous("log10 mtDNA copy number") +
  scale_x_continuous("m.3243A>G level") +
  geom_abline(intercept = mod_p_het$coefficients$fixed[1], slope = mod_p_het$coefficients$fixed[2], size=1, alpha=0.8) 

tiff(filename="For paper/output/cn_age.tiff", width=520, height=400)
cn_age
dev.off()

tiff(filename="For paper/output/cn_het.tiff", width=520, height=400)
cn_het
dev.off()

tiff(filename="For paper/output/log_cn.tiff", width=1100, height=480)
logged_cn
dev.off()



#do the same with the cnvhet per cell type graph



# het vs cn

all_p = subset(all, type == "Patient")

cn_het_big = ggplot(subset(all_p, !is.na(normalised_cn)), aes(x=het, y=log_cn)) + theme_classic() +
  geom_point(aes(colour=maturity1), size=3) +
  geom_smooth(aes(colour=maturity1, fill=maturity1), method="lm", alpha=0.2) +
  ylab("Normalised mtDNA Copy Number") +
  xlab("m.3243A>G Level (%)") +
  facet_wrap(~broad_cell_type) + 
  stat_cor(size=4, size=12, aes(colour=maturity1, label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                                                        ifelse(round(..p.., 3)==0, "p < 0.001" , paste("p",round(..p.., 3),sep="*' = '*")),
                                                        sep = "*', '*")), label.y=c(230, 185), key_glyph = draw_key_rect) +
  scale_colour_manual(breaks = c("Naive/Immature", "Memory/Mature"), values=c("#377eb8", "#e41a1c","grey35","pink"))  +
  scale_fill_manual(breaks = c("Naive/Immature", "Memory/Mature"), values=c("#377eb8", "#e41a1c","grey35", "pink")) +
  theme(
    plot.title = element_text(size=30, face="bold.italic"),
    axis.title.x = element_text(size=30),
    axis.title.y = element_text(size=30),
    axis.text=element_text(size=30),
    legend.title = element_blank(),
    legend.text = element_text(size=30),
    strip.text.x = element_text(size = 30),
    legend.position="bottom"
  )

#png("output/cn_het_big.png", width=1200, height=1800)
#cn_het_big
#dev.off()

#make a logged one 

logcn_het_cell_type = 
  ggplot(subset(all, !is.na(normalised_cn)), aes(x=het, y=log_cn)) + my_theme +
  geom_point() +
  geom_smooth(method="lm") +
  #geom_point(aes(colour=maturity2)) +
  #geom_smooth(aes(colour=maturity2, fill=maturity2), method="lm", alpha=0.2) +
  #geom_point(data=subset(all, !is.na(normalised_cn) & is.na(maturity2)), aes(x=het, y=log_cn)) +
  #geom_smooth(data=subset(all, !is.na(normalised_cn) & is.na(maturity2)), alpha=0.2, method="lm", aes(x=het, y=log_cn)) +
  ylab("log10(Normalised mtDNA Copy Number)") +
  xlab("m.3243A>G Level (%)") +
  facet_wrap(~verbose_broad_cell_type, ncol=4) + 
  #stat_cor(aes(colour = maturity2, label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
  #                                                  ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
  #                                                  sep = "*', '*")),  label.y = c(3, 2.7), key_glyph = draw_key_rect)+
  #stat_cor(data=subset(all, is.na(maturity2)),aes(label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
  #                                                                      ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
  #                                                                       sep = "*', '*")), label.y = 2.7, key_glyph = draw_key_rect) +
  stat_cor(aes(label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                                                                        ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
                                                                         sep = "*', '*")), label.y = 2.7, key_glyph = draw_key_rect) 
  #scale_color_npg(name = "", na.translate = F) +
  #scale_fill_npg(name = "", na.translate = F) +
  #theme(legend.position = c(0.95, 0.1),legend.justification = c(0.95, 0.1))
 

logcn_age_cell_type = 
  ggplot(subset(all, !is.na(normalised_cn)), aes(x=age, y=log_cn)) + my_theme +
  geom_point(aes(colour=type)) +
  geom_smooth(aes(colour=type, fill=type), method="lm", alpha=0.2) +
  ylab("log10 mtDNA Copy Number") +
  xlab("Age") +
  facet_wrap(~verbose_broad_cell_type, ncol=4) + 
  stat_cor(aes(colour = type, label = paste(paste("italic(R)^2",round(..rr.., 2),sep="*' = '*"), 
                                                 ifelse(round(..p.., 4)==0, "p < 0.0001" , paste("p",round(..p.., 4),sep="*' = '*")),
                                                 sep = "*', '*")),  label.y = c(3, 2.7), key_glyph = draw_key_rect)+
  scale_color_npg(name = "", na.translate = F) +
  scale_fill_npg(name = "", na.translate = F) +
  theme(legend.position = c(0.95, 0.1),legend.justification = c(0.95, 0.1))

tiff(filename="output/sarah/logcn_het_cell_type.tiff", width=1100, height=600)
logcn_het_cell_type
dev.off()

tiff(filename="output/sarah/logcn_age_cell_type.tiff", width=1100, height=600)
logcn_age_cell_type
dev.off()

# stats for cell subtype CN correlations
sink(file="output/sarah/age_cn_cell_type_interaction.txt")
for(cell in unique(all$verbose_broad_cell_type)){
mod_age_cn = lm(log_cn ~ age*type, data=subset(all, !is.na(log_cn) & verbose_broad_cell_type==cell))
print(cell)
print(summary(mod_age_cn))
}
sink()

sink(file="output/sarah/age_cn_cell_type.txt")
for(cell in unique(all$verbose_broad_cell_type)){
  for(group in unique(all$type)){
    mod_age_cn = lm(log_cn ~ age, data=subset(all, !is.na(log_cn) & verbose_broad_cell_type==cell & type==group))
    print(c(cell, group))
    print(summary(mod_age_cn))
}
}
sink()


sink(file="output/sarah/het_cn_cell_type.txt")
for(cell in unique(all$verbose_broad_cell_type)){
  mod_het_cn = lm(log_cn ~ het, data=subset(all, !is.na(log_cn) & !is.na(het) & verbose_broad_cell_type==cell))
  print(cell)
  print(summary(mod_het_cn))
}
sink()


