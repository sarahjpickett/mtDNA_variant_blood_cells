## settings and load data ####

library(ggplot2)

# patients - het and cn data
het = read.csv(file="For paper/input/all_data_v6.csv")
head(het)
names(het)
het$type = "Patient"

# create a dataframe that shows wbcs numbers and P numbers for paper
WBC = (unique(het$WBC))[order(unique(het$WBC))]
pnums = sprintf("P%02d",1:length(WBC))
samples = data.frame(WBC, pnums)
het = merge(het, samples)

### control cn data
controls = read.csv(file="For paper/input/all_cn_10_11_21.csv")
head(controls)
names(controls)
controls_clean = controls[,c("patient", "absolute_cn","cell_type", "broad_cell_type", "control_cn", "cell_order","age")]
controls_clean = plyr::rename(controls_clean, c("patient" = "patient_id", "absolute_cn" = "cn"))
controls_clean$type = "Control"
head(controls_clean)
head(het)


unique(het$cell_type)

het$cell_type = plyr::revalue(het$cell_type, c(
  "CD34+" = "CD34+", 
  "Naive Mature B cell" = "Naive B", 
  "Memory B cell" = "Memory B", 
  "CD4+ Naive" = "CD4+ Naive T", 
  "CD8+ Naive" = "CD8+ Naive T",
  "CD4+ Memory" = "CD4+ Memory T", 
  "CD8+ Memory" = "CD8+ Memory T", 
  "CD14+ Monocyte" = "CD14+ Monocytes", 
  "CD16+ Monocyte" = "CD16+ Monocytes", 
  "CD56+ CD16++ NK Cell" = "Mature NK",
  "CD56++ CD16- NK Cell" = "Immature NK", 
  "11C+ Dendritic Cell" = "CD11c+ mDCs", 
  "123+ Dendritic Cell" = "CD123+ pDCs",
  "Neutrophils" = "Granulocytes"
))

het$maturity = plyr::revalue(het$maturity, c(
  "Mature" = "Memory/Mature",
  "Naive" = "Naive/Immature"))


controls_clean$cell_type = plyr::revalue(controls_clean$cell_type, c(
  "CD34+" = "CD34+", 
  "Naive Mature B cell" = "Naive B", 
  "Memory B cell" = "Memory B", 
  "CD4+ Naive" = "CD4+ Naive T", 
  "CD8+ Naive" = "CD8+ Naive T",
  "CD4+ Memory" = "CD4+ Memory T", 
  "CD8+ Memory" = "CD8+ Memory T", 
  "CD14+ Monocyte" = "CD14+ Monocytes", 
  "CD16+ Monocyte" = "CD16+ Monocytes", 
  "CD56+ CD16++ NK Cell" = "Mature NK",
  "CD56++ CD16- NK Cell" = "Immature NK", 
  "11C+ Dendritic Cell" = "CD11c+ mDCs", 
  "123+ Dendritic Cell" = "CD123+ pDCs"
))

dim(het)
dim(controls_clean)
all = merge(het, controls_clean, all=T)
dim(all)
head(all)
all= subset(all, cell_type!= "Neutrophils" & cell_type!= "Platelets" & cell_type!= "Whole Blood")

all$broad_cell_type2 = all$cell_type

all$broad_cell_type2 = factor(all$broad_cell_type2, levels = c("CD34+", "Naïve B","Memory B", "CD4+ Naïve T", "CD8+ Naïve T", "CD4+ Memory T", "CD8+ Memory T",
                                                               "Immature NK", "Mature NK", "CD14+ Monocytes", "CD16+ Monocytes","CD11c+ mDCs", "CD123+ pDCs"))

mean(unique(all$control_cn), na.rm=T)
subset(all, is.na(control_cn))
     
unique(controls$control_cn)
unique(het$control_cn)

## scale copy number according to plate controls ####
mean_control_cn = mean(unique(controls$control_cn))
controls$normalised_cn = (controls$absolute_cn / controls$control_cn)*mean_control_cn 

het_clean = het
## scale copy number according to plate controls ####

mean_ccn = mean(unique(het_clean$control_cn), na.rm = TRUE)

table(het_clean$control_cn)

het_clean$normalised_cn = (het_clean$cn / het_clean$control_cn)*mean_ccn


# Maturity

het_clean$maturity1 = ifelse(het_clean$maturity=="Naive/Immature", "Naive/Immature", ifelse(
  het_clean$maturity=="Memory/Mature", "Memory/Mature", ""
))

het_clean$maturity1 =  factor(het_clean$maturity1, levels=c("Naive/Immature", "Memory/Mature", ""))


#remove 4+ & 8+ CM from het cleam

het_clean = subset(het_clean, cell_type != "8+ CM" & cell_type != "4+ CM")

# and one for cn without platelets and neutrophils
cn_clean = subset(het_clean, cell_type!="Platelets" & cell_type!="Granulocytes" & cell_type!="Whole Blood")

# re order for nice graph

class(het_clean$cell_type)

het_clean$cell_type <- as.factor(het_clean$cell_type)

class(het_clean$cell_type)

levels(het_clean$cell_type)

het_clean$cell_type <- factor(het_clean$cell_type, levels=c("CD34+", "Naive B",
                                                            "Memory B", "CD4+ Naive T", "CD4+ Memory T",
                                                            "CD8+ Naive T", "CD8+ Memory T", "Immature NK",
                                                            "Mature NK", "CD14+ Monocytes", "CD16+ Monocytes", "CD11c+ mDCs",
                                                            "CD123+ pDCs", "Granulocytes", "Platelets", "Whole Blood"))


#add in column which depicts wb het

het_clean$WBHET =  ifelse(het_clean$patient_id=="R002", "24",
                          ifelse(het_clean$patient_id=="R004", "21",
                                 ifelse(het_clean$patient_id=="R006", "20",
                                        ifelse(het_clean$patient_id=="R007" , "18",
                                               ifelse(het_clean$patient_id=="R008", "6",
                                                      ifelse(het_clean$patient_id=="R009", "13", 
                                                             ifelse(het_clean$patient_id=="R010", "10", 
                                                                    ifelse(het_clean$patient_id=="R011", "32",
                                                                           ifelse(het_clean$patient_id=="R013", "59",
                                                                                  ifelse(het_clean$patient_id=="R016", "11",
                                                                                         ifelse(het_clean$patient_id=="R017", "32",
                                                                                                ifelse(het_clean$patient_id=="R023", "31",
                                                                                                       ifelse(het_clean$patient_id=="R024", "23",
                                                                                                              ifelse(het_clean$patient_id=="R025", "14",
                                                                                                                     ifelse(het_clean$patient_id=="R026", "23",
                                                                                                                            ifelse(het_clean$patient_id=="R027", "9",
                                                                                                                                   ifelse(het_clean$patient_id=="R028", "36",
                                                                                                                                          ifelse(het_clean$patient_id=="R044", "59",
                                                                                                                                                 ifelse(het_clean$patient_id=="R017", "32",
                                                                                                                                                        ifelse(het_clean$patient_id=="R067", "44",
                                                                                                                                                               ifelse(het_clean$patient_id=="R081", "12",
                                                                                                                                                                      ifelse(het_clean$patient_id=="R0149", "60",
                                                                                                                                                                             NA))))))))))))))))))))))


class(het_clean$WBHET)

het_clean$WBHET <- as.numeric(het_clean$WBHET)


#normalise all data by the mean control on all plate using controls

class(all$control_cn)

mean_control_cn = mean(unique(all$control_cn), na.rm = TRUE)

all$normalised_cn = (all$cn/ all$control_cn)*mean_control_cn

#remove cell types with no cn reading

all = subset(all, cell_type != "4+ CM" & cell_type != "8+ CM" & cell_type != "Granulocytes")

#reorder cell type levels

class(all$cell_type)

all$cell_type <- as.factor(all$cell_type)

class(all$cell_type)

levels(all$cell_type)

all$cell_type <- factor(all$cell_type, levels=c("CD34+", "Naive B",
                                                "Memory B", "CD4+ Naive T", "CD4+ Memory T",
                                                "CD8+ Naive T", "CD8+ Memory T", "Immature NK",
                                                "Mature NK", "CD14+ Monocytes", "CD16+ Monocytes", "CD11c+ mDCs",
                                                "CD123+ pDCs"))






