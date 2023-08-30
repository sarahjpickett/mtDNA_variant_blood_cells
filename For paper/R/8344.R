library(ggplot2)
library(forcats)
library(dplyr)
#library(tidyverse)
library(patchwork)
library(nlme)

my_theme = theme_classic() + theme(text = element_text(size = 16))

# load in data ##

data = read.csv(file="For paper/input/All_8344.csv", stringsAsFactors = FALSE)
head(data)

# remove controls #

data_final = subset(data, Cell_type != "1567" & Cell_type != "502" & Cell_type != "2021" & Cell_type != "2022" & Cell_type != "Negative" & Cell_type != "CD4+ CM" & Cell_type != "CD8+ CM" & Cell_type != "CD8+ TEMRA")

data_final$Cell_type = ifelse(data_final$Cell_type=="CD34+ Progenitors", "CD34+",
                              ifelse(data_final$Cell_type=="CD56++ CD16- NK Cell","Immature NK",
                                     ifelse(data_final$Cell_type=="CD56+ CD16++ NK Cell","Mature NK",
                                            ifelse(data_final$Cell_type=="mDC","CD11c+ mDCs",
                                                   ifelse(data_final$Cell_type=="pDC","CD123+ pDCs",
                                                          ifelse(data_final$Cell_type=="CD16+ Monocyte","CD16+ Monocytes",
                                                                 ifelse(data_final$Cell_type=="Granulocyte","Granulocytes",
                                                                        ifelse(data_final$Cell_type=="Platelet","Platelets",
                                                                          ifelse(data_final$Cell_type=="NM B-cell","Naive B", data_final$Cell_type)))))))))
                                     
data_final$Cell_type = factor(data_final$Cell_type, levels=c("CD34+", "Naive B", "Memory B", "CD4+ Naive","CD4+ Memory", "CD8+ Naive", 
                                                             "CD8+ Memory", "Immature NK", "Mature NK", "CD14+ Monocytes", "CD16+ Monocytes",
                                                             "CD11c+ mDCs", "CD123+ pDCs","Granulocytes","Platelets", "PBMC", "WB"))


#remove DCs

minus_pdc = subset(data_final, Cell_type != "CD123+ pDCs")

#separate PBMC and WB measurements into a different dataset like in 3243 graphs. remove cell types we dont have for all patients

unique(minus_pdc$Cell_type)
m8344_clean = subset(minus_pdc, Cell_type != "PBMC" & Cell_type != "WB" & Cell_type != "CD11c+ mDCs")

forpaper =  ggplot(m8344_clean, aes(x=Cell_type, y=Het)) + my_theme +
  geom_point(aes( colour = Patient_ID)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylim(0,100) +
  geom_line(data = m8344_clean, aes(group = Patient_ID, colour = Patient_ID), size = 0.7) +
  ylab("") +
  xlab("")


final = forpaper +  scale_colour_discrete(name = "Patient", labels = c("23", "24", "25", "26"))


pdf("For paper/output/8344final.pdf", width=7, height=9, useDingbats = F)
final
dev.off()


box_8344 =  ggplot(m8344_clean, aes(x=Cell_type, y=Het)) + my_theme +
  geom_boxplot(outlier.shape=NA) +
  geom_point(aes(colour = Age)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylim(0,100) +
  ylab("m.8344A>G level") +
  xlab("") + 
  scale_colour_gradient(high="red", low="royalblue") 

box_8344

pdf("For paper/output/8344final_box.pdf", width=7, height=9, useDingbats = F)
box_8344
dev.off()

pdf("For paper/output/test.pdf", width=7, height=9, useDingbats = F)
final_het
dev.off()


#lets see if we have enough patients to see anything in the stats

#convert to factor and change levels

#m8344_clean$Cell_type = factor(m8344_clean$Cell_type, levels=c("CD34+ Progenitors", "CD4+ Naive","CD4+ Memory", "CD8+ Naive", 
#                                                               "CD8+ Memory", "NM B-cell", "Memory B", "CD56++ CD16- NK Cell",
#                                                               "CD56+ CD16++ NK Cell", "CD14+ Monocytes", "CD16+ Monocyte",
#                                                               "Platelet", "Granulocyte"))

mod = lme(Het ~ Cell_type, random = ~ 1|Patient_ID, data=subset(m8344_clean, !is.na(Het)))
summary(mod)

