library(ggplot2)
library(forcats)
library(dplyr)
library(tidyverse)
library(patchwork)
library(nlme)

# load in data ##

data = read.csv(file="For paper/input/All_8344.csv", stringsAsFactors = FALSE)
head(data)

# remove controls #

data_final = subset(data, Cell_type != "1567" & Cell_type != "502" & Cell_type != "2021" & Cell_type != "2022" & Cell_type != "Negative" & Cell_type != "CD4+ CM" & Cell_type != "CD8+ CM" & Cell_type != "CD8+ TEMRA")

#re order on graph

p <- data_final %>%
  mutate(Cell_type = fct_relevel(Cell_type, 
                                 "CD34+ Progenitors", "CD4+ Naive","CD4+ Memory", "CD8+ Naive", 
                                 "CD8+ Memory", "NM B-cell", "Memory B", "CD56++ CD16- NK Cell",
                                 "CD56+ CD16++ NK Cell", "CD14+ Monocytes", "CD16+ Monocyte", "mDC", "pDC",
                                 "Platelet", "Granulocyte")) %>%
  ggplot(aes(x=Cell_type, y=Het)) + theme_bw() +
  geom_point(aes( colour = Patient_ID)) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size = 15)) +
  ylim(0,100) +
  ylab("Heteroplasmy") +
  theme(axis.title = element_text(size = 15)) +
  ggtitle("m.8344A>G") +
  xlab("") 
#theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p 

#make pdf

png("For paper/output/8344bigtext.png", width=800, height=550)
p
dev.off()


#join with lines


p2 <- data_final %>%
  mutate(Cell_type = fct_relevel(Cell_type, 
                                 "CD34+ Progenitors", "CD4+ Naive","CD4+ Memory", "CD8+ Naive", 
                                 "CD8+ Memory", "NM B-cell", "Memory B", "CD56++ CD16- NK Cell",
                                 "CD56+ CD16++ NK Cell", "CD14+ Monocytes", "CD16+ Monocyte", "mDC", "pDC",
                                 "Platelet", "Granulocyte")) %>%
  ggplot(aes(x=Cell_type, y=Het)) + theme_classic() +
  geom_point(aes( colour = Patient_ID)) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size = 12)) +
  ylim(0,100) +
  geom_line(data = data_final, aes(group = Patient_ID, colour = Patient_ID), size = 1) +
  ylab("Heteroplasmy") +
  theme(axis.title = element_text(size = 15)) +
  ggtitle("m.8344A>G") +
  xlab("") 

p2


png("Output/8344lines.png", width=800, height=550)
p2
dev.off()

#remove pDC and redraw? just label as DC category?

minus_pdc = subset(data_final, Cell_type != "pDC")


#separate PBMC and WB measurements into a different dataset like in 3243 graphs. remove cell types we dont have for all patients

unique(minus_pdc$Cell_type)
m8344_clean = subset(minus_pdc, Cell_type != "PBMC" & Cell_type != "WB" & Cell_type != "DC")

forpaper <- m8344_clean %>%
  mutate(Cell_type = fct_relevel(Cell_type, 
                                 "CD34+ Progenitors", "CD4+ Naive","CD4+ Memory", "CD8+ Naive", 
                                 "CD8+ Memory", "NM B-cell", "Memory B", "CD56++ CD16- NK Cell",
                                 "CD56+ CD16++ NK Cell", "CD14+ Monocytes", "CD16+ Monocyte",
                                 "Platelet", "Granulocyte")) %>%
  ggplot(aes(x=Cell_type, y=Het)) + theme_classic() +
  geom_point(aes( colour = Patient_ID)) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  ylim(0,100) +
  geom_line(data = m8344_clean, aes(group = Patient_ID, colour = Patient_ID), size = 1) +
  ylab("") +
  theme(axis.title = element_text(size = 15)) +
  xlab("") 


final = forpaper +  scale_colour_discrete(name = "Patient", labels = c("23", "24", "25", "26"))


png("For paper/output/8344final.png", width=800, height=550)
p
dev.off()


#lets see if we have enough patients to see anything in the stats

#convert to factor and change levels

m8344_clean$Cell_type = factor(m8344_clean$Cell_type, levels=c("CD34+ Progenitors", "CD4+ Naive","CD4+ Memory", "CD8+ Naive", 
                                                               "CD8+ Memory", "NM B-cell", "Memory B", "CD56++ CD16- NK Cell",
                                                               "CD56+ CD16++ NK Cell", "CD14+ Monocytes", "CD16+ Monocyte",
                                                               "Platelet", "Granulocyte"))

mod = lme(Het ~ Cell_type, random = ~ 1|Patient_ID, data=subset(m8344_clean, !is.na(Het)))
summary(mod)
