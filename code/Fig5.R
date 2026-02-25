### Figure 5 ###
### The codes are separated by Figures ###

### load packages ###
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(reshape2)
library(ggalluvial)
library(gtools)

### Figure 5A, 5G and 5M were created with Biorender ###

#### Figure 5B ####
crc_meta_to_ST = read_excel("./data_meta/ST_CRC_cohort_meta.xlsx")
crc_meta_to_ST = crc_meta_to_ST[crc_meta_to_ST$Sample_ID != 'TMA',]

crc_meta = read.csv("./data_meta/ST_CRC_cohort_meta2.csv")

ST_meta22 = read.csv("./data_meta/ST_CRC_cohort_meta3.csv")

ST_meta2 <- ST_meta22 %>%
  group_by(Patient) %>%
  summarize(
    APC.status = ifelse(any(APC.status == "APC mut."), "mutant", "WT"),
    KRAS.status = ifelse(any(KRAS.status == "KRAS mut."), "mutant", "WT"),
    TP53.status = ifelse(any(TP53.status == "TP53 mut."), "mutant", "WT"),
    BRAF.status = ifelse(any(BRAF.status == "BRAF mut."), "mutant", "WT"))

ST_meta2 = data.frame(ST_meta2)
ST_meta2 = ST_meta2[,c('Patient',"APC.status","KRAS.status","BRAF.status")]

crc_meta_combined_long = merge(crc_meta,crc_meta_to_ST,by.x='sample_key',by.y="Sample_ID")

ST_meta2_filtered = ST_meta2[ST_meta2$Patient %in% crc_meta_combined_long$patient_name,]
crc_meta_combined_long2 = merge(crc_meta_combined_long,ST_meta2_filtered,by.x='patient_name',by.y="Patient",all=T)

crc_meta = crc_meta[,c("sample_key",'patient_name','tumor_loc',"tumor_grade","tumor_type","CIN.Status","mets")]
crc_meta$mets = ifelse(crc_meta$mets == "",'Mets_No','Mets_Yes')
crc_meta_combined = merge(crc_meta,crc_meta_to_ST,by.x='sample_key',by.y="Sample_ID")
ST_meta2_filtered = ST_meta2[ST_meta2$Patient %in% crc_meta_combined$patient_name,]

crc_meta_combined2 = merge(crc_meta_combined,ST_meta2_filtered,by.x='patient_name','Patient',all=T)

tls_classification = readRDS("./data/TLS_classification.rds")
tls_classification$Sample_ID = sapply(strsplit(tls_classification$TLS_ID, "_"), function(x) paste(x[1]))

tls_classification2 = tls_classification[tls_classification$Sample_ID %in% crc_meta_combined2$ST_ID,]

tls_classification3 = tls_classification2 %>% group_by(Sample_ID,Cluster) %>% summarise(count = n()) %>% pivot_wider(names_from = Cluster, values_from = count, values_fill = 0)
tls_classification3 = data.frame(tls_classification3)
colnames(tls_classification3) = c("Sample_ID","PTLS","ETLS", "STLS")


crc_tls_count_with_meta = merge(crc_meta_combined2,tls_classification3,by.x='ST_ID',by.y='Sample_ID',all=T)

crc_tls_count_with_meta$tumor_side = ifelse(crc_tls_count_with_meta$tumor_loc %in% c('Ascending','Cecum','Hepatic Flexure','Transverse'),'Right_side',
                                            ifelse(crc_tls_count_with_meta$tumor_loc %in% c('Descending','Sigmoid'),'Left_side',
                                                   ''))

crc_tls_count_with_meta$ETLS = ifelse(is.na(crc_tls_count_with_meta$ETLS),0,crc_tls_count_with_meta$ETLS)
crc_tls_count_with_meta$PTLS = ifelse(is.na(crc_tls_count_with_meta$PTLS),0,crc_tls_count_with_meta$PTLS)
crc_tls_count_with_meta$STLS = ifelse(is.na(crc_tls_count_with_meta$STLS),0,crc_tls_count_with_meta$STLS)
crc_tls_count_with_meta$Total_counts = rowSums(crc_tls_count_with_meta[,c("PTLS","ETLS", "STLS")])
crc_tls_count_with_meta = crc_tls_count_with_meta %>% arrange(tumor_type,desc(Total_counts),patient_name)

crc_tls_count_with_meta = crc_tls_count_with_meta %>%
  group_by(patient_name) %>%
  mutate(Total_counts_sample = max(Total_counts)) %>%
  ungroup()

crc_tls_count_with_meta = crc_tls_count_with_meta %>% 
  arrange(tumor_side,desc(Total_counts_sample),patient_name,desc(Total_counts))


crc_tls_count_with_meta_long = crc_tls_count_with_meta %>% 
  arrange(desc(tumor_side),desc(Total_counts_sample),patient_name,desc(Total_counts),tumor_loc,tumor_grade) %>% 
  pivot_longer(cols = c(ETLS,PTLS,STLS),names_to = "TLS", values_to = "TLS_Count")
crc_tls_count_with_meta_long = data.frame(crc_tls_count_with_meta_long)

crc_tls_count_with_meta_long$ST_ID = factor(crc_tls_count_with_meta_long$ST_ID,levels = unique(crc_tls_count_with_meta_long$ST_ID))

samples <- c(
  "CO38","CO37","CO26","CO25","CO31","CO20","CO22","CO19","CO21",
  "CO52","CO51","CO39","CO44","CO53","CO54","CO33","CO40","CO32",
  "CO34","CO36","CO45","CO47","CO48","CO49","CO50","CO30","CO29",
  "CO27","CO28","CO41","CO35","CO42","CO43","CO23","CO24","CO46"
)

sum(as.character(crc_tls_count_with_meta_long$ST_ID) != rep(samples,each=3))

a= ggplot(crc_tls_count_with_meta_long, aes(x =ST_ID, y = TLS_Count, fill = TLS)) + 
  geom_bar(stat = "identity") + theme_classic() +
  scale_fill_manual(values = rev(c("#B3B1D8","#F16623","#71B28B"))) + 
  scale_y_continuous(breaks = seq(0, 10, 2),limits = c(0, 10), expand = c(0, 0))

pdf("./result/Fig5B_1.pdf",width = 15,height = 5)
print(a)
dev.off()

crc_tls_count_with_meta$mets = factor(crc_tls_count_with_meta$mets,levels = c('Mets_No','Mets_Yes'),
                                      labels=c('No Mets','Mets'))
crc_tls_count_with_meta_heatmap = crc_tls_count_with_meta
crc_tls_count_with_meta_heatmap = data.frame(crc_tls_count_with_meta_heatmap)
rownames(crc_tls_count_with_meta_heatmap) = crc_tls_count_with_meta_heatmap$ST_ID
crc_tls_count_with_meta_heatmap2 = crc_tls_count_with_meta_heatmap[,-which(colnames(crc_tls_count_with_meta_heatmap) %in% c('ST_ID','sample_key',"PTLS","ETLS","STLS","Total_counts","Total_counts_sample"))]
crc_tls_count_with_meta_heatmap2 = crc_tls_count_with_meta_heatmap2[,c('tumor_side','tumor_loc','tumor_grade','tumor_type','mets','KRAS.status','patient_name')]
crc_tls_count_with_meta_heatmap2 = crc_tls_count_with_meta_heatmap2[order(crc_tls_count_with_meta_heatmap2$tumor_side,decreasing = T),]
crc_tls_count_with_meta_heatmap2 = t(crc_tls_count_with_meta_heatmap2)

rownames(crc_tls_count_with_meta_heatmap2) <- c(
  "Tumor side",
  "Tumor location",
  "Tumor grade",
  "Tumor type",
  "Metastasis",
  "KRAS status",
  "Patient ID"
)
colors = c(
  "Right_side"="#79AF97FF",
  "Left_side"="#374E55FF",
  'Ascending' = "#B24745FF",
  'Cecum' = "#DF8F44FF",
  'Descending' = 'blue',
  'Hepatic Flexure' = "#80796BFF",
  'Sigmoid' = "#00A1D5FF",
  'Transverse' = "#6A6599FF",
  
  "NL" = "#D51317FF",
  "G1" = "#FFD9D9",
  "G2" = "#DE6464",
  "G3" = "#771213",
  
  "MSS" = "#F39200FF",
  "MSI-H" = "#EFD500FF",
  "SSL/HP" = "#95C11FFF",
  "TA/TVA" = "#007B3DFF",
  
  
  "No Mets" = "#fbb040",
  "Mets" = "#5e3f18",
  
  'mutant'='#CFDFF2', 
  'WT'='#17365E',
  
  "PAT54273" = "#d22753",
  "PAT73458" = "#913a1f",
  "HTA11_08622_A" = "#8e5725",
  "PAT71397" = "#977947",
  "SG00003" = "#ca6128",
  "HTA11_10711" = "#bc9f6f",
  "PAT59667" = "#ef7322",
  "SG00004" = "#f68526",
  "HTA11_01938" = "#ceba98",
  "HTA11_07663" = "#d79c66",
  "HTA11_08622_B" = "#e7c3a4",
  "PAT33430" = "#da6038",
  "PAT59460" = "#f57f20",
  "PAT59600" = "#fdc695",
  "SG00001" = "#e28362",
  "SG00002" = "#faece6",
  "PAT73899" = "#693514",
  "PAT71662" = "#f89a48",
  "HTA11_06134" = "#bc7332",
  "HTA11_07862" = "#e7decd",
  "PAT30884" = "#faaf6e",
  "PAT74143" = "#ffefe3",
  "SR00001" = "#d47e2e"
)

hmap2 <- Heatmap(crc_tls_count_with_meta_heatmap2, 
                 col = colors,
                 na_col = 'white', 
                 rect_gp = grid::gpar(col = NA),
                 row_title = NULL, 
                 height = unit(6*0.75,'cm'), 
                 width = unit(20*0.75, 'cm'), 
                 row_names_side = 'left', 
                 column_title = NULL)

pdf("./result/Fig5B_2.pdf",width = 15,height = 10)
draw(hmap2)
dev.off()
rm(list = ls())

#### Figure 5C ####
CRC_anno_TLS_loc = read.csv("./data/ST_CRC_maturation_location.csv")

CRC_cont_table <- table(CRC_anno_TLS_loc$Cluster, CRC_anno_TLS_loc$Location)

CRC_prop_table <- prop.table(CRC_cont_table, margin = 1)
CRC_prop_table_melted <- melt(CRC_prop_table)

CRC_prop_table_melted$Var2 = factor(CRC_prop_table_melted$Var2,
                                    levels = c('Normal mucosa','Smooth muscle','Adenoma','Carcinoma'))
CRC_prop_table_melted$Var1 = factor(CRC_prop_table_melted$Var1,
                                    levels = c('Secondary Mature','Primary Mature','Immature'),labels=c("S-TLS","P-TLS","E-TLS"))

cols <- c(
  "Carcinoma" = "#dca1a1",         
  "Adenoma" = "#7c9885",           
  "Smooth muscle" = "#d6ecf3",     
  "Normal mucosa" = "#f8b636"
)


b = ggplot(CRC_prop_table_melted, aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(position="stack", stat="identity") + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))  + 
  scale_fill_manual(values = cols) + coord_flip()

pdf("./result/Fig5C.pdf")
print(b)
dev.off()
rm(list = ls())

#### Figure 5D ####
CRC_anno_TLS_loc = read.csv("./data/ST_CRC_maturation_location.csv")
CRC_cont_table <- table(CRC_anno_TLS_loc$Cluster, CRC_anno_TLS_loc$Location)

CRC_observed <- as.matrix(CRC_cont_table)
CRC_expected <- chisq.test(CRC_observed)$expected
CRC_ro_e <- CRC_observed / CRC_expected
rownames(CRC_ro_e) = c('E-TLS','P-TLS','S-TLS')
CRC_ro_e_melted <- melt(CRC_ro_e)
CRC_ro_e_melted$Var2 = factor(CRC_ro_e_melted$Var2,
                              levels = rev(c('Normal mucosa','Smooth muscle','Adenoma','Carcinoma')))

b = ggplot(CRC_ro_e_melted, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#559B80", mid = "white", high = "red4", midpoint = 1) +
  labs(x = "Location", y = "Cluster", fill = "Ro/e") +
  coord_flip() +
  theme_minimal()

pdf("./result/Fig5D.pdf")
print(b)
dev.off()
rm(list = ls())

#### Figure 5E ####
crc_meta_to_ST = read_excel("./data_meta/ST_CRC_cohort_meta.xlsx")
crc_meta_to_ST = crc_meta_to_ST[crc_meta_to_ST$Sample_ID != 'TMA',]

crc_meta = read.csv("./data_meta/ST_CRC_cohort_meta2.csv")

ST_meta22 = read.csv("./data_meta/ST_CRC_cohort_meta3.csv")

ST_meta2 <- ST_meta22 %>%
  group_by(Patient) %>%
  summarize(
    APC.status = ifelse(any(APC.status == "APC mut."), "mutant", "WT"),
    KRAS.status = ifelse(any(KRAS.status == "KRAS mut."), "mutant", "WT"),
    TP53.status = ifelse(any(TP53.status == "TP53 mut."), "mutant", "WT"),
    BRAF.status = ifelse(any(BRAF.status == "BRAF mut."), "mutant", "WT"))

ST_meta2 = data.frame(ST_meta2)
ST_meta2 = ST_meta2[,c('Patient',"APC.status","KRAS.status","BRAF.status")]

crc_meta_combined_long = merge(crc_meta,crc_meta_to_ST,by.x='sample_key',by.y="Sample_ID")

ST_meta2_filtered = ST_meta2[ST_meta2$Patient %in% crc_meta_combined_long$patient_name,]
crc_meta_combined_long2 = merge(crc_meta_combined_long,ST_meta2_filtered,by.x='patient_name',by.y="Patient",all=T)


crc_meta = crc_meta[,c("sample_key",'patient_name',"tumor_type","tumor_grade","CIN.Status","mets")]
crc_meta_combined = merge(crc_meta,crc_meta_to_ST,by.x='sample_key',by.y="Sample_ID")
ST_meta2_filtered = ST_meta2[ST_meta2$Patient %in% crc_meta_combined$patient_name,]

crc_meta_combined2 = merge(crc_meta_combined,ST_meta2_filtered,by.x='patient_name','Patient',all=T)

tls_classification = readRDS("./data/TLS_classification.rds")
tls_classification$Sample_ID = sapply(strsplit(tls_classification$TLS_ID, "_"), function(x) paste(x[1]))

tls_classification2 = tls_classification[tls_classification$Sample_ID %in% crc_meta_combined2$ST_ID,]

tls_classification3 = tls_classification2 %>% group_by(Sample_ID,Cluster) %>% summarise(count = n()) %>% pivot_wider(names_from = Cluster, values_from = count, values_fill = 0)
tls_classification3 = data.frame(tls_classification3)
colnames(tls_classification3) = c("Sample_ID","PTLS","ETLS", "STLS")


crc_tls_count_with_meta = merge(crc_meta_combined2,tls_classification3,by.x='ST_ID',by.y='Sample_ID',all=T)

crc_tls_count_with_meta$ETLS = ifelse(is.na(crc_tls_count_with_meta$ETLS),0,crc_tls_count_with_meta$ETLS)
crc_tls_count_with_meta$PTLS = ifelse(is.na(crc_tls_count_with_meta$PTLS),0,crc_tls_count_with_meta$PTLS)
crc_tls_count_with_meta$STLS = ifelse(is.na(crc_tls_count_with_meta$STLS),0,crc_tls_count_with_meta$STLS)
crc_tls_count_with_meta$Total_counts = rowSums(crc_tls_count_with_meta[,c("PTLS","ETLS", "STLS")])

crc_tls_count_with_meta = crc_tls_count_with_meta[crc_tls_count_with_meta$tumor_type != 'NL',]

crc_tls_count_with_meta = crc_tls_count_with_meta <- crc_tls_count_with_meta %>%
  group_by(patient_name) %>%
  summarise(
    ETLS = sum(ETLS),
    PTLS = sum(PTLS),
    STLS = sum(STLS),
    Total_counts = sum(Total_counts),
    tumor_grade = first(tumor_grade),
    .groups = 'drop'
  )

crc_tls_count_with_meta$TLS_class = ifelse(rowSums(crc_tls_count_with_meta[,c('ETLS','PTLS','STLS')]) == 0,"No_TLS",
                                           apply(crc_tls_count_with_meta[, c('ETLS','PTLS','STLS')], 1, function(x) {
                                             categories <- c('E-TLS','P-TLS','S-TLS')
                                             categories[which.max(x)]
                                           }))

crc_tls_count_with_meta$TLS_class = factor(crc_tls_count_with_meta$TLS_class,levels = c("No_TLS","E-TLS","P-TLS","S-TLS"))

crc_tls_count_with_meta$TLS_class2 <-  ifelse(rowSums(crc_tls_count_with_meta[,c('ETLS','PTLS','STLS')]) == 0,"No_TLS",
                                              apply(crc_tls_count_with_meta[, c('ETLS', 'PTLS', 'STLS')], 1, function(x) {
                                                categories <- c('ETLS', 'PTLS', 'STLS')
                                                max_value <- max(x)
                                                max_indices <- which(x == max_value)
                                                tied_categories <- categories[max_indices]
                                                
                                                if ("STLS" %in% tied_categories) {
                                                  return("STLS")
                                                } else if ("PTLS" %in% tied_categories) {
                                                  return("PTLS")
                                                } else {
                                                  return("ETLS")
                                                }
                                              }))

crc_tls_count_with_meta$TLS_class2 = ifelse(crc_tls_count_with_meta$TLS_class2 == 'ETLS','E-TLS',
                                            ifelse(crc_tls_count_with_meta$TLS_class2 == 'PTLS','P-TLS',
                                                   ifelse(crc_tls_count_with_meta$TLS_class2 == 'STLS','S-TLS','No_TLS')))

crc_tls_count_with_meta$TLS_class2 = factor(crc_tls_count_with_meta$TLS_class2,levels = c("No_TLS","E-TLS","P-TLS","S-TLS")) 

crc_tls_count_with_meta_wo_no_TLS = crc_tls_count_with_meta[crc_tls_count_with_meta$TLS_class2 != "No_TLS",]
tumor_grade_tls_prop = prop.table(table(crc_tls_count_with_meta_wo_no_TLS$tumor_grade,crc_tls_count_with_meta_wo_no_TLS$TLS_class2),margin = 1)

tumor_grade_tls_df <- as.data.frame(tumor_grade_tls_prop)

colnames(tumor_grade_tls_df) <- c("Tumor_Grade", "TLS_Class", "Proportion")

tumor_grade_tls_df$TLS_Class = as.character(tumor_grade_tls_df$TLS_Class)
tumor_grade_tls_df = tumor_grade_tls_df[tumor_grade_tls_df$TLS_Class != 'No_TLS',]
a=ggplot(tumor_grade_tls_df, aes(x = Tumor_Grade, y = Proportion, fill = TLS_Class)) +
  geom_bar(stat = "identity", position = "stack") +      
  scale_fill_manual(values = c("E-TLS" = "#71B28B", "P-TLS" ="#F06825","S-TLS"="#B3B1D8")) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))  + 
  theme_classic()

pdf("./result/Fig5E.pdf")
print(a)
dev.off()
rm(list = ls())

#### Figure 5F ####
crc_meta_to_ST = read_excel("./data_meta/ST_CRC_cohort_meta.xlsx")
crc_meta_to_ST = crc_meta_to_ST[crc_meta_to_ST$Sample_ID != 'TMA',]

crc_meta = read.csv("./data_meta/ST_CRC_cohort_meta2.csv")

ST_meta22 = read.csv("./data_meta/ST_CRC_cohort_meta3.csv")

ST_meta2 <- ST_meta22 %>%
  group_by(Patient) %>%
  summarize(
    APC.status = ifelse(any(APC.status == "APC mut."), "mutant", "WT"),
    KRAS.status = ifelse(any(KRAS.status == "KRAS mut."), "mutant", "WT"),
    TP53.status = ifelse(any(TP53.status == "TP53 mut."), "mutant", "WT"),
    BRAF.status = ifelse(any(BRAF.status == "BRAF mut."), "mutant", "WT"))

ST_meta2 = data.frame(ST_meta2)
ST_meta2 = ST_meta2[,c('Patient',"APC.status","KRAS.status","BRAF.status")]

crc_meta_combined_long = merge(crc_meta,crc_meta_to_ST,by.x='sample_key',by.y="Sample_ID")

ST_meta2_filtered = ST_meta2[ST_meta2$Patient %in% crc_meta_combined_long$patient_name,]
crc_meta_combined_long2 = merge(crc_meta_combined_long,ST_meta2_filtered,by.x='patient_name',by.y="Patient",all=T)


crc_meta = crc_meta[,c("sample_key",'patient_name',"tumor_type","tumor_grade","CIN.Status","mets")]
crc_meta_combined = merge(crc_meta,crc_meta_to_ST,by.x='sample_key',by.y="Sample_ID")
ST_meta2_filtered = ST_meta2[ST_meta2$Patient %in% crc_meta_combined$patient_name,]

crc_meta_combined2 = merge(crc_meta_combined,ST_meta2_filtered,by.x='patient_name','Patient',all=T)

tls_classification = readRDS("./data/TLS_classification.rds")
tls_classification$Sample_ID = sapply(strsplit(tls_classification$TLS_ID, "_"), function(x) paste(x[1]))

tls_classification2 = tls_classification[tls_classification$Sample_ID %in% crc_meta_combined2$ST_ID,]

tls_classification3 = tls_classification2 %>% group_by(Sample_ID,Cluster) %>% summarise(count = n()) %>% pivot_wider(names_from = Cluster, values_from = count, values_fill = 0)
tls_classification3 = data.frame(tls_classification3)
colnames(tls_classification3) = c("Sample_ID","PTLS","ETLS", "STLS")

crc_tls_count_with_meta = merge(crc_meta_combined2,tls_classification3,by.x='ST_ID',by.y='Sample_ID',all=T)

crc_tls_count_with_meta$ETLS = ifelse(is.na(crc_tls_count_with_meta$ETLS),0,crc_tls_count_with_meta$ETLS)
crc_tls_count_with_meta$PTLS = ifelse(is.na(crc_tls_count_with_meta$PTLS),0,crc_tls_count_with_meta$PTLS)
crc_tls_count_with_meta$STLS = ifelse(is.na(crc_tls_count_with_meta$STLS),0,crc_tls_count_with_meta$STLS)
crc_tls_count_with_meta$Total_counts = rowSums(crc_tls_count_with_meta[,c("PTLS","ETLS", "STLS")])

crc_tls_count_with_meta = crc_tls_count_with_meta[crc_tls_count_with_meta$tumor_type != 'NL',]

crc_tls_count_with_meta = crc_tls_count_with_meta <- crc_tls_count_with_meta %>%
  group_by(patient_name) %>%
  summarise(
    ETLS = sum(ETLS),
    PTLS = sum(PTLS),
    STLS = sum(STLS),
    Total_counts = sum(Total_counts),
    tumor_type = first(tumor_type),
    tumor_grade = first(tumor_grade),
    CIN.Status = first(CIN.Status),
    APC.status = first(APC.status),
    KRAS.status = first(KRAS.status),
    BRAF.status = first(BRAF.status),
    mets = first(mets),
    .groups = 'drop'
  )

crc_tls_count_with_meta$TLS_class = ifelse(rowSums(crc_tls_count_with_meta[,c('ETLS','PTLS','STLS')]) == 0,"No_TLS",
                                           apply(crc_tls_count_with_meta[, c('ETLS','PTLS','STLS')], 1, function(x) {
                                             categories <- c('E-TLS','P-TLS','S-TLS')
                                             categories[which.max(x)]
                                           }))

crc_tls_count_with_meta$TLS_class = factor(crc_tls_count_with_meta$TLS_class,levels = c("No_TLS","E-TLS","P-TLS","S-TLS"))

crc_tls_count_with_meta$TLS_class2 <-  ifelse(rowSums(crc_tls_count_with_meta[,c('ETLS','PTLS','STLS')]) == 0,"No_TLS",
                                              apply(crc_tls_count_with_meta[, c('ETLS', 'PTLS', 'STLS')], 1, function(x) {
                                                categories <- c('ETLS', 'PTLS', 'STLS')
                                                max_value <- max(x)
                                                max_indices <- which(x == max_value)
                                                tied_categories <- categories[max_indices]
                                                
                                                if ("STLS" %in% tied_categories) {
                                                  return("STLS")
                                                } else if ("PTLS" %in% tied_categories) {
                                                  return("PTLS")
                                                } else {
                                                  return("ETLS")
                                                }
                                              }))

crc_tls_count_with_meta$TLS_class2 = ifelse(crc_tls_count_with_meta$TLS_class2 == 'ETLS','E-TLS',
                                            ifelse(crc_tls_count_with_meta$TLS_class2 == 'PTLS','P-TLS',
                                                   ifelse(crc_tls_count_with_meta$TLS_class2 == 'STLS','S-TLS','No_TLS')))

crc_tls_count_with_meta$TLS_class2 = factor(crc_tls_count_with_meta$TLS_class2,levels = c("No_TLS","E-TLS","P-TLS","S-TLS")) 
crc_tls_count_with_meta$mets = ifelse(crc_tls_count_with_meta$mets == "",'Mets_No','Mets_Yes')
var_names2=c("KRAS.status","mets")

crc_tls_count_with_meta_KRAS_filtered = crc_tls_count_with_meta[!is.na(crc_tls_count_with_meta[["KRAS.status"]]),]

p <- ggplot(crc_tls_count_with_meta_KRAS_filtered,
            aes(axis1 = TLS_class2, axis2 = KRAS.status, y = 1)) +
  geom_alluvium(aes(fill = TLS_class2), width = 1/6) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/6) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("TLS_class2", "KRAS.status"), expand = c(0.15, 0.15)) +
  scale_fill_manual(values = c(
    "No_TLS" = "#f7e6ca",
    "E-TLS"  = "#71B28B",
    "P-TLS"  = "#F06825",
    "S-TLS"  = "#B3B1D8",
    "WT"     = "#17365e",
    "mutant"    = "#cfdff2"
  )) +
  labs(title = '') +
  theme_void()

pdf("./result/Fig5F_1.pdf")
print(p)
dev.off()

crc_tls_count_with_meta_mets_filtered = crc_tls_count_with_meta[!is.na(crc_tls_count_with_meta[["mets"]]),]
crc_tls_count_with_meta_mets_filtered$mets = factor(crc_tls_count_with_meta_mets_filtered$mets,levels = c('Mets_No','Mets_Yes'),
                                                    labels=c('No Mets','Mets'))
p <- ggplot(crc_tls_count_with_meta_mets_filtered,
            aes(axis1 = TLS_class2, axis2 = mets, y = 1)) +
  geom_alluvium(aes(fill = TLS_class2), width = 1/6) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/6) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("TLS_class2", "mets"), expand = c(0.15, 0.15)) +
  scale_fill_manual(values = c(
    "No_TLS" = "#f7e6ca",
    "E-TLS"  = "#71B28B",
    "P-TLS"  = "#F06825",
    "S-TLS"  = "#B3B1D8",
    "No Mets"     = "#fbb040",
    "Mets"    = "#5f401a"
  )) +
  labs(title = '') +
  theme_void()

pdf("./result/Fig5F_2.pdf")
print(p)
dev.off()
rm(list = ls())

#### Figure 5H ####
tls_classification = readRDS("./data/TLS_classification.rds")
tls_classification$Sample_ID = sapply(strsplit(tls_classification$TLS_ID, "_"), function(x) paste(x[1]))

lihc_meta = read_excel('./data_meta/ST_LIHC_cohort_meta.xlsx')
lihc_meta$Stage = ifelse(lihc_meta$Stage %in% c("Ia","Ib"),'I',lihc_meta$Stage)
lihc_meta = lihc_meta[,colnames(lihc_meta) %in% c('Patient_ID','Etiology_of_liver_disease','BCLC_stage')]

lihc_meta_with_ST_ID2 = read_excel("./data_meta/ST_LIHC_cohort_meta2.xlsx")
lihc_meta_with_ST_ID2 = lihc_meta_with_ST_ID2[,c('Sample_ID','Patient_ID','Sample Description','Sample_ID2')]
lihc_meta_with_ST_ID2 = data.frame(lihc_meta_with_ST_ID2)

tls_classification2 = tls_classification[tls_classification$Sample_ID %in% lihc_meta_with_ST_ID2$Sample_ID2,]

tls_classification3 = tls_classification2 %>% group_by(Sample_ID,Cluster) %>% summarise(count = n()) %>% pivot_wider(names_from = Cluster, values_from = count, values_fill = 0)
tls_classification3 = data.frame(tls_classification3)
colnames(tls_classification3) = c("Sample_ID","ETLS","PTLS","STLS")

lihc_meta_combined_with_TLS = merge(lihc_meta_with_ST_ID2,tls_classification3,by.x='Sample_ID2',by.y='Sample_ID',all=T)
colnames(lihc_meta_combined_with_TLS)[1] = 'ST_ID'

lihc_meta_combined_with_TLS = merge(lihc_meta_combined_with_TLS,lihc_meta,by='Patient_ID',all=T)

lihc_meta_combined_with_TLS$ETLS = ifelse(is.na(lihc_meta_combined_with_TLS$ETLS),0,lihc_meta_combined_with_TLS$ETLS)
lihc_meta_combined_with_TLS$PTLS = ifelse(is.na(lihc_meta_combined_with_TLS$PTLS),0,lihc_meta_combined_with_TLS$PTLS)
lihc_meta_combined_with_TLS$STLS = ifelse(is.na(lihc_meta_combined_with_TLS$STLS),0,lihc_meta_combined_with_TLS$STLS)
lihc_meta_combined_with_TLS$Total_counts = rowSums(lihc_meta_combined_with_TLS[,c("PTLS","ETLS","STLS")])
lihc_meta_combined_with_TLS$Sample_ID2 = paste0(lihc_meta_combined_with_TLS$Sample_ID,"_",lihc_meta_combined_with_TLS$Sample.Description)
lihc_meta_combined_with_TLS = lihc_meta_combined_with_TLS %>%
  group_by(Sample_ID2) %>%
  mutate(Total_counts_sample = sum(Total_counts)) %>%
  ungroup()


order_vec <- c("Leading-edge", "Tumor", "Normal", "Portal vein tumor thrombus")

lihc_meta_combined_with_TLS <-
  lihc_meta_combined_with_TLS %>%
  arrange(
    match(Sample.Description, order_vec),
    desc(Total_counts_sample),
    Sample_ID,
    desc(Total_counts)
  )

lihc_meta_combined_with_TLS_long = lihc_meta_combined_with_TLS %>% 
  pivot_longer(cols = c(ETLS,PTLS,STLS),names_to = "TLS", values_to = "TLS_Count")
lihc_meta_combined_with_TLS_long = data.frame(lihc_meta_combined_with_TLS_long)

lihc_meta_combined_with_TLS_long$ST_ID = factor(lihc_meta_combined_with_TLS_long$ST_ID,levels = unique(lihc_meta_combined_with_TLS_long$ST_ID))


a1 = ggplot(lihc_meta_combined_with_TLS_long, aes(x=ST_ID, y=TLS_Count, fill = TLS)) + 
  geom_bar(stat = "identity") + theme_classic() +
  scale_fill_manual(values = rev(c("#B3B1D8","#F16623","#71B28B"))) + 
  scale_y_continuous(limits = c(0, 15), expand = c(0, 0))

pdf("./result/Fig5H_1.pdf",width = 15,height = 5)
print(a1)
dev.off()

lihc_meta_combined_with_TLS_heatmap = lihc_meta_combined_with_TLS
lihc_meta_combined_with_TLS_heatmap = data.frame(lihc_meta_combined_with_TLS_heatmap)
lihc_meta_combined_with_TLS_heatmap$BCLC_stage = ifelse(lihc_meta_combined_with_TLS_heatmap$BCLC_stage == 'B',
                                                        'B/C',lihc_meta_combined_with_TLS_heatmap$BCLC_stage)
lihc_meta_combined_with_TLS_heatmap$BCLC_stage = ifelse(lihc_meta_combined_with_TLS_heatmap$BCLC_stage == 'C',
                                                        'B/C',lihc_meta_combined_with_TLS_heatmap$BCLC_stage)
lihc_meta_combined_with_TLS_heatmap$BCLC_stage = ifelse(lihc_meta_combined_with_TLS_heatmap$BCLC_stage == '0',
                                                        '0/A',lihc_meta_combined_with_TLS_heatmap$BCLC_stage)
lihc_meta_combined_with_TLS_heatmap$BCLC_stage = ifelse(lihc_meta_combined_with_TLS_heatmap$BCLC_stage == 'A',
                                                        '0/A',lihc_meta_combined_with_TLS_heatmap$BCLC_stage)
lihc_meta_combined_with_TLS_heatmap$Etiology_of_liver_disease = factor(lihc_meta_combined_with_TLS_heatmap$Etiology_of_liver_disease,
                                                                       levels=c('HBV','non-HBV'),labels=c('HBV+','HBV-'))
lihc_meta_combined_with_TLS_heatmap$Sample.Description = ifelse(lihc_meta_combined_with_TLS_heatmap$Sample.Description == "Portal vein tumor thrombus",'Portal vein',lihc_meta_combined_with_TLS_heatmap$Sample.Description)
rownames(lihc_meta_combined_with_TLS_heatmap) = lihc_meta_combined_with_TLS_heatmap$ST_ID
lihc_meta_combined_with_TLS_heatmap2 = lihc_meta_combined_with_TLS_heatmap[,-which(colnames(lihc_meta_combined_with_TLS_heatmap) %in% c('ST_ID','Sample.Name',"Histology_Type","PTLS","ETLS","STLS","Total_counts","Total_counts_sample"))]
lihc_meta_combined_with_TLS_heatmap2 = lihc_meta_combined_with_TLS_heatmap2[,c('Sample.Description','BCLC_stage','Etiology_of_liver_disease','Patient_ID')]
lihc_meta_combined_with_TLS_heatmap2 = t(lihc_meta_combined_with_TLS_heatmap2)

colors = c(    "HBV+"     = "#18375e",
               "HBV-"    = "#cfddef",
               "0/A"     = "#fad8d8",
               "B/C"    = "#771215",
               "Normal" = "#d4ead4",
               "Leading-edge" = "#68ba7f",
               "Tumor" = "#2e7040",
               "Portal vein" = "#cfdff2",
               "HCC-3"="#9268ad",
               "cHC-1"="#2278b5",
               "HCC-4"="#3a3c7a",
               "HCC-1"="#f57f20",
               "HCC-2"="#2fa148",
               "ICC-1"="#bcbe32",
               "HCC-5"="#d62a28")

hmap2 <- Heatmap(lihc_meta_combined_with_TLS_heatmap2, 
                 col = colors,
                 na_col = 'grey60', 
                 rect_gp = grid::gpar(col = NA),
                 row_title = NULL, 
                 height = unit(6*0.75,'cm'), 
                 width = unit(20*0.75, 'cm'), 
                 row_names_side = 'left', 
                 column_title = NULL)

pdf("./result/Fig5H_2.pdf",width = 15,height = 10)
draw(hmap2)
dev.off()
rm(list = ls())

#### Figure 5I ####
LIHC_anno_TLS_loc = read.csv("./data/ST_LIHC_maturation_location.csv")
LIHC_anno_TLS_loc = LIHC_anno_TLS_loc[LIHC_anno_TLS_loc$Location != 'Out_of_tissue',]

LIHC_cont_table <- table(LIHC_anno_TLS_loc$Cluster, LIHC_anno_TLS_loc$Location)

LIHC_prop_table <- prop.table(LIHC_cont_table, margin = 1)
colnames(LIHC_prop_table)[which(colnames(LIHC_prop_table) == 'Normal')] = 'Non-tumor'
LIHC_prop_table_melted <- melt(LIHC_prop_table)

LIHC_prop_table_melted$Var2 = factor(LIHC_prop_table_melted$Var2,
                                     levels = c('Non-tumor','Tumor_Boundary','Tumor_region'))
LIHC_prop_table_melted$Var1 = factor(LIHC_prop_table_melted$Var1,
                                     levels = c('Secondary Mature','Primary Mature','Immature'),labels=c("S-TLS","P-TLS","E-TLS"))

cols <- c(
  "Tumor_Boundary" = "#a97c50",
  "Non-tumor" = "#f8b636",
  "Tumor_region" = "#dca1a1"
)


a = ggplot(LIHC_prop_table_melted, aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(position="stack", stat="identity") + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))  + 
  scale_fill_manual(values = cols ) + coord_flip()

pdf("./result/Fig5I.pdf")
print(a)
dev.off()
rm(list = ls())

#### Figure 5J ####
LIHC_anno_TLS_loc = read.csv("./data/ST_LIHC_maturation_location.csv")
LIHC_anno_TLS_loc = LIHC_anno_TLS_loc[LIHC_anno_TLS_loc$Location != 'Out_of_tissue',]

LIHC_cont_table <- table(LIHC_anno_TLS_loc$Cluster, LIHC_anno_TLS_loc$Location)

LIHC_observed <- as.matrix(LIHC_cont_table)
LIHC_expected <- chisq.test(LIHC_observed)$expected
LIHC_ro_e <- LIHC_observed / LIHC_expected
rownames(LIHC_ro_e) = c('E-TLS','P-TLS','S-TLS')
colnames(LIHC_ro_e)[which(colnames(LIHC_ro_e) == 'Normal')] = 'Non-tumor'

LIHC_ro_e_melted <- melt(LIHC_ro_e)
LIHC_ro_e_melted$Var2 = factor(LIHC_ro_e_melted$Var2,
                               levels = rev(c('Non-tumor','Tumor_Boundary','Tumor_region')))


a = ggplot(LIHC_ro_e_melted, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#559B80", mid = "white", high = "red4", midpoint = 1) +
  labs(x = "Location", y = "Cluster", fill = "Ro/e") +
  coord_flip() +
  theme_minimal()

pdf("./result/Fig5J.pdf")
print(a)
dev.off()
rm(list = ls())

#### Figure 5K ####
tls_classification = readRDS("./data/TLS_classification.rds")
tls_classification$Sample_ID = sapply(strsplit(tls_classification$TLS_ID, "_"), function(x) paste(x[1]))

lihc_meta = read_excel('./data_meta/ST_LIHC_cohort_meta.xlsx')
lihc_meta$Stage = ifelse(lihc_meta$Stage %in% c("Ia","Ib"),'I',lihc_meta$Stage)
lihc_meta = lihc_meta[,colnames(lihc_meta) %in% c('Patient_ID','Etiology_of_liver_disease','BCLC_stage')]

lihc_meta_with_ST_ID2 = read_excel("./data_meta/ST_LIHC_cohort_meta2.xlsx")
lihc_meta_with_ST_ID2 = lihc_meta_with_ST_ID2[,c('Sample_ID','Patient_ID','Sample Description','Sample_ID2')]
lihc_meta_with_ST_ID2 = data.frame(lihc_meta_with_ST_ID2)

tls_classification2 = tls_classification[tls_classification$Sample_ID %in% lihc_meta_with_ST_ID2$Sample_ID2,]

tls_classification3 = tls_classification2 %>% group_by(Sample_ID,Cluster) %>% summarise(count = n()) %>% pivot_wider(names_from = Cluster, values_from = count, values_fill = 0)
tls_classification3 = data.frame(tls_classification3)
colnames(tls_classification3) = c("Sample_ID","ETLS","PTLS","STLS")

lihc_meta_combined_with_TLS = merge(lihc_meta_with_ST_ID2,tls_classification3,by.x='Sample_ID2',by.y='Sample_ID')
colnames(lihc_meta_combined_with_TLS)[1] = 'ST_ID'

lihc_meta_combined_with_TLS = merge(lihc_meta_combined_with_TLS,lihc_meta,by='Patient_ID',all=T)

lihc_meta_combined_with_TLS$ETLS = ifelse(is.na(lihc_meta_combined_with_TLS$ETLS),0,lihc_meta_combined_with_TLS$ETLS)
lihc_meta_combined_with_TLS$PTLS = ifelse(is.na(lihc_meta_combined_with_TLS$PTLS),0,lihc_meta_combined_with_TLS$PTLS)
lihc_meta_combined_with_TLS$STLS = ifelse(is.na(lihc_meta_combined_with_TLS$STLS),0,lihc_meta_combined_with_TLS$STLS)

lihc_meta_combined_with_TLS$Total_TLS = rowSums(lihc_meta_combined_with_TLS[,c('ETLS','PTLS','STLS')])
lihc_meta_combined_with_TLS$ETLS_prop = lihc_meta_combined_with_TLS$ETLS / lihc_meta_combined_with_TLS$Total_TLS
lihc_meta_combined_with_TLS$PTLS_prop = lihc_meta_combined_with_TLS$PTLS / lihc_meta_combined_with_TLS$Total_TLS
lihc_meta_combined_with_TLS$STLS_prop = lihc_meta_combined_with_TLS$STLS / lihc_meta_combined_with_TLS$Total_TLS

lihc_HBV = lihc_meta_combined_with_TLS %>% group_by(Etiology_of_liver_disease) %>% summarise(across(c(ETLS_prop,PTLS_prop,STLS_prop),mean))
lihc_HBV_long = lihc_HBV %>% pivot_longer(cols= c(ETLS_prop,PTLS_prop,STLS_prop),names_to = "TLS", values_to = "Prop")

lihc_HBV_long$Etiology_of_liver_disease = factor(lihc_HBV_long$Etiology_of_liver_disease,levels = c("non-HBV","HBV"),labels=c("HBV-","HBV+"))
lihc_HBV_long$TLS = factor(lihc_HBV_long$TLS,levels = c("ETLS_prop","PTLS_prop","STLS_prop"),labels=c("E-TLS","P-TLS","S-TLS"))

a1 = ggplot(lihc_HBV_long, aes(x =Etiology_of_liver_disease, y = Prop, fill = TLS)) + 
  geom_bar(stat = "identity") + theme_classic() + labs(x = "")+
  scale_fill_manual(values = rev(c("#B3B1D8","#F16623","#71B28B"))) + 
  scale_x_discrete(limits=rev) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))

pdf("./result/Fig5K.pdf")
print(a1)
dev.off()
rm(list = ls())

#### Figure 5L ####
tls_classification = readRDS("./data/TLS_classification.rds")
tls_classification$Sample_ID = sapply(strsplit(tls_classification$TLS_ID, "_"), function(x) paste(x[1]))

lihc_meta = read_excel('./data_meta/ST_LIHC_cohort_meta.xlsx')
lihc_meta$Stage = ifelse(lihc_meta$Stage %in% c("Ia","Ib"),'I',lihc_meta$Stage)
lihc_meta = lihc_meta[,colnames(lihc_meta) %in% c('Patient_ID','Etiology_of_liver_disease','BCLC_stage')]

lihc_meta_with_ST_ID2 = read_excel("./data_meta/ST_LIHC_cohort_meta2.xlsx")
lihc_meta_with_ST_ID2 = lihc_meta_with_ST_ID2[,c('Sample_ID','Patient_ID','Sample Description','Sample_ID2')]
lihc_meta_with_ST_ID2 = data.frame(lihc_meta_with_ST_ID2)

tls_classification2 = tls_classification[tls_classification$Sample_ID %in% lihc_meta_with_ST_ID2$Sample_ID2,]

tls_classification3 = tls_classification2 %>% group_by(Sample_ID,Cluster) %>% summarise(count = n()) %>% pivot_wider(names_from = Cluster, values_from = count, values_fill = 0)
tls_classification3 = data.frame(tls_classification3)
colnames(tls_classification3) = c("Sample_ID","ETLS","PTLS","STLS")

lihc_meta_combined_with_TLS = merge(lihc_meta_with_ST_ID2,tls_classification3,by.x='Sample_ID2',by.y='Sample_ID')
colnames(lihc_meta_combined_with_TLS)[1] = 'ST_ID'

lihc_meta_combined_with_TLS = merge(lihc_meta_combined_with_TLS,lihc_meta,by='Patient_ID',all=T)

lihc_meta_combined_with_TLS$BCLC_stage = ifelse(lihc_meta_combined_with_TLS$BCLC_stage == 'B',
                                                'B/C',lihc_meta_combined_with_TLS$BCLC_stage)
lihc_meta_combined_with_TLS$BCLC_stage = ifelse(lihc_meta_combined_with_TLS$BCLC_stage == 'C',
                                                'B/C',lihc_meta_combined_with_TLS$BCLC_stage)
lihc_meta_combined_with_TLS$BCLC_stage = ifelse(lihc_meta_combined_with_TLS$BCLC_stage == '0',
                                                '0/A',lihc_meta_combined_with_TLS$BCLC_stage)
lihc_meta_combined_with_TLS$BCLC_stage = ifelse(lihc_meta_combined_with_TLS$BCLC_stage == 'A',
                                                '0/A',lihc_meta_combined_with_TLS$BCLC_stage)

lihc_meta_combined_with_TLS$ETLS = ifelse(is.na(lihc_meta_combined_with_TLS$ETLS),0,lihc_meta_combined_with_TLS$ETLS)
lihc_meta_combined_with_TLS$PTLS = ifelse(is.na(lihc_meta_combined_with_TLS$PTLS),0,lihc_meta_combined_with_TLS$PTLS)
lihc_meta_combined_with_TLS$STLS = ifelse(is.na(lihc_meta_combined_with_TLS$STLS),0,lihc_meta_combined_with_TLS$STLS)


lihc_meta_combined_with_TLS$TLS_class = ifelse(rowSums(lihc_meta_combined_with_TLS[,c('ETLS','PTLS','STLS')]) == 0,"No_TLS",
                                               apply(lihc_meta_combined_with_TLS[, c('ETLS', 'PTLS', 'STLS')], 1, function(x) {
                                                 categories <- c('ETLS', 'PTLS', 'STLS')
                                                 categories[which.max(x)]
                                               }))

lihc_meta_combined_with_TLS$TLS_class = ifelse(lihc_meta_combined_with_TLS$TLS_class == 'ETLS','E-TLS',
                                               ifelse(lihc_meta_combined_with_TLS$TLS_class == 'PTLS','P-TLS',
                                                      'S-TLS'))

lihc_meta_combined_with_TLS$Etiology_of_liver_disease = factor(lihc_meta_combined_with_TLS$Etiology_of_liver_disease,
                                                               levels=c('HBV','non-HBV'),labels=c('HBV+','HBV-'))

p <- ggplot(lihc_meta_combined_with_TLS,
            aes(axis1 = TLS_class, axis2 = Etiology_of_liver_disease, y = 1)) +
  geom_alluvium(aes(fill = TLS_class), width = 1/6) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/6) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("TLS_class2", "Etiology_of_liver_disease"), expand = c(0.15, 0.15)) +
  scale_fill_manual(values = c(
    "No_TLS" = "#f7e6ca",
    "E-TLS"  = "#71B28B",
    "P-TLS"  = "#F06825",
    "S-TLS"  = "#B3B1D8",
    "HBV+"     = "#18375e",
    "HBV-"    = "#cfddef"
  )) +
  labs(title = '') +
  theme_void()

pdf("./result/Fig5L_1.pdf")
print(p)
dev.off()

p <- ggplot(lihc_meta_combined_with_TLS,
            aes(axis1 = TLS_class, axis2 = BCLC_stage, y = 1)) +
  geom_alluvium(aes(fill = TLS_class), width = 1/6) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/6) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("TLS_class2", "BCLC_stage"), expand = c(0.15, 0.15)) +
  scale_fill_manual(values = c(
    "No_TLS" = "#f7e6ca",
    "E-TLS"  = "#71B28B",
    "P-TLS"  = "#F06825",
    "S-TLS"  = "#B3B1D8",
    "0/A"     = "#fad8d8",
    "B/C"    = "#771215"
  )) +
  labs(title = '') +
  theme_void()

pdf("./result/Fig5L_2.pdf")
print(p)
dev.off()
rm(list = ls())

#### Figure 5N ####
tls_classification = readRDS("./data/TLS_classification.rds")
tls_classification$Sample_ID = sapply(strsplit(tls_classification$TLS_ID, "_"), function(x) paste(x[1]))

bladder_sample_info = read_excel("./data_meta/ST_BLCA_cohort_meta.xlsx")
bladder_sample_info$`cStage(pre)` <- substr(bladder_sample_info$`cStage(pre)`, 1, 3)
bladder_sample_info$MTAP = ifelse(bladder_sample_info$MTAP == "P","Pos","Neg")

tls_classification2 = tls_classification[tls_classification$Sample_ID %in% bladder_sample_info$Sample_ID,]

tls_classification3 = tls_classification2 %>% group_by(Sample_ID,Cluster) %>% summarise(count = n()) %>% pivot_wider(names_from = Cluster, values_from = count, values_fill = 0)
tls_classification3 = data.frame(tls_classification3)

colnames(tls_classification3) = c("Sample_ID","ETLS","PTLS")

blca_meta_combined_with_TLS = merge(bladder_sample_info,tls_classification3,by.x = 'Sample_ID',by.y = 'Sample_ID',all=T)
colnames(blca_meta_combined_with_TLS)[1] = 'ST_ID'

blca_meta_combined_with_TLS$ETLS = ifelse(is.na(blca_meta_combined_with_TLS$ETLS),0,blca_meta_combined_with_TLS$ETLS)
blca_meta_combined_with_TLS$PTLS = ifelse(is.na(blca_meta_combined_with_TLS$PTLS),0,blca_meta_combined_with_TLS$PTLS)

blca_meta_combined_with_TLS$Total_counts = rowSums(blca_meta_combined_with_TLS[,c("ETLS", "PTLS")])
blca_meta_combined_with_TLS = blca_meta_combined_with_TLS %>% arrange(MTAP,desc(Total_counts),ST_ID)

blca_meta_combined_with_TLS$Responder = factor(blca_meta_combined_with_TLS$Responder,levels=c("N","Y"),labels=c("NR","R"))
blca_meta_combined_with_TLS$MTAP =  factor(blca_meta_combined_with_TLS$MTAP,levels = c("Neg","Pos"),labels=c("MTAP-","MTAP+"))

order_vec <- c("MTAP+", "MTAP-")

blca_meta_combined_with_TLS = blca_meta_combined_with_TLS %>% 
  arrange(match(MTAP, order_vec),desc(Total_counts),ST_ID)

blca_meta_combined_with_TLS_long = blca_meta_combined_with_TLS %>% 
  pivot_longer(cols = c(ETLS,PTLS),names_to = "TLS", values_to = "TLS_Count")

blca_meta_combined_with_TLS_long = data.frame(blca_meta_combined_with_TLS_long)

blca_meta_combined_with_TLS_long$ST_ID = factor(blca_meta_combined_with_TLS_long$ST_ID,levels = unique(blca_meta_combined_with_TLS_long$ST_ID))


a1 = ggplot(blca_meta_combined_with_TLS_long, aes(x=ST_ID, y=TLS_Count, fill = TLS)) + 
  geom_bar(stat = "identity") + theme_classic() +
  scale_fill_manual(values = rev(c("#B3B1D8","#F16623","#71B28B"))) + 
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0))

pdf("./result/Fig5N_1.pdf",width = 15,height = 5)
print(a1)
dev.off()

blca_meta_combined_with_TLS_heatmap = blca_meta_combined_with_TLS
blca_meta_combined_with_TLS_heatmap = data.frame(blca_meta_combined_with_TLS_heatmap)
rownames(blca_meta_combined_with_TLS_heatmap) = blca_meta_combined_with_TLS_heatmap$ST_ID
blca_meta_combined_with_TLS_heatmap2 = blca_meta_combined_with_TLS_heatmap[,-which(colnames(blca_meta_combined_with_TLS_heatmap) %in% c('ST_ID',"ETLS","PTLS","Total_counts"))]
blca_meta_combined_with_TLS_heatmap2 = blca_meta_combined_with_TLS_heatmap2[,c("MTAP","Responder","cStage.pre.")]
blca_meta_combined_with_TLS_heatmap2 = t(blca_meta_combined_with_TLS_heatmap2)
colors = c( "MTAP+"= "#17365e",
            "MTAP-"= "#cfdff2",
            "NR"   = "#e0f0e1",
            "R"    = "#0f8c44",
            "cT2"  = "#FFD9D9",
            "cT3"  = "#DE6464",
            "cT4"  = "#780606")

hmap2 <- Heatmap(blca_meta_combined_with_TLS_heatmap2, 
                 col = colors,
                 na_col = 'grey60', 
                 rect_gp = grid::gpar(col = NA),
                 row_title = NULL, 
                 height = unit(6*0.75,'cm'), 
                 width = unit(20*0.75, 'cm'), 
                 row_names_side = 'left', 
                 column_title = NULL)

pdf("./result/Fig5N_2.pdf",width = 15,height = 10)
draw(hmap2)
dev.off()
rm(list = ls())

#### Figure 5O ####
BLCA_anno_TLS_loc = read.csv("./data/ST_BLCA_maturation_location.csv")
BLCA_anno_TLS_loc$Location = ifelse(BLCA_anno_TLS_loc$Location == 'Normal_tissue','Normal',
                                    BLCA_anno_TLS_loc$Location)
BLCA_anno_TLS_loc$Location = ifelse(BLCA_anno_TLS_loc$Location == 'Normal stroma','Stroma',
                                    BLCA_anno_TLS_loc$Location)
BLCA_anno_TLS_loc = BLCA_anno_TLS_loc[BLCA_anno_TLS_loc$Location != 'Out_of_tissue',]

BLCA_cont_table <- table(BLCA_anno_TLS_loc$Cluster, BLCA_anno_TLS_loc$Location)

BLCA_prop_table <- prop.table(BLCA_cont_table, margin = 1)
BLCA_prop_table_melted <- melt(BLCA_prop_table)
BLCA_prop_table_melted$Var1 = factor(BLCA_prop_table_melted$Var1,
                                     levels = c('Primary Mature','Immature'),labels=c("P-TLS","E-TLS"))

cols <- c(
  "Stroma" = "#fdfbd4",
  "Tumor" = "#b22222"
)


c = ggplot(BLCA_prop_table_melted, aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(position="stack", stat="identity") + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))  + 
  scale_fill_manual(values = cols) + coord_flip()

pdf("./result/Fig5O.pdf")
print(c)
dev.off()

rm(list = ls())

#### Figure 5P ####
BLCA_anno_TLS_loc = read.csv("./data/ST_BLCA_maturation_location.csv")
BLCA_cont_table <- table(BLCA_anno_TLS_loc$Cluster, BLCA_anno_TLS_loc$Location)

BLCA_observed <- as.matrix(BLCA_cont_table)
BLCA_expected <- chisq.test(BLCA_observed)$expected
BLCA_ro_e <- BLCA_observed / BLCA_expected
rownames(BLCA_ro_e) = c('E-TLS','P-TLS')
BLCA_ro_e_melted <- melt(BLCA_ro_e)

BLCA_ro_e_melted$Var2 = factor(BLCA_ro_e_melted$Var2,
                               levels = rev(c('Stroma','Tumor')))


c = ggplot(BLCA_ro_e_melted, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#559B80", mid = "white", high = "red4", midpoint = 1) +
  labs(x = "Location", y = "Cluster", fill = "Ro/e") +
  coord_flip() + 
  theme_minimal()

pdf("./result/Fig5P.pdf")
print(c)
dev.off()
rm(list = ls())

#### Figure 5Q ####
tls_classification = readRDS("./data/TLS_classification.rds")
tls_classification$Sample_ID = sapply(strsplit(tls_classification$TLS_ID, "_"), function(x) paste(x[1]))

bladder_sample_info = read_excel("./data_meta/ST_BLCA_cohort_meta.xlsx")
bladder_sample_info$`cStage(pre)` <- substr(bladder_sample_info$`cStage(pre)`, 1, 3)
bladder_sample_info$MTAP = ifelse(bladder_sample_info$MTAP == "P","Pos","Neg")

tls_classification2 = tls_classification[tls_classification$Sample_ID %in% bladder_sample_info$Sample_ID,]

tls_classification3 = tls_classification2 %>% group_by(Sample_ID,Cluster) %>% summarise(count = n()) %>% pivot_wider(names_from = Cluster, values_from = count, values_fill = 0)
tls_classification3 = data.frame(tls_classification3)

colnames(tls_classification3) = c("Sample_ID","E-TLS","P-TLS")

blca_meta_combined_with_TLS = merge(bladder_sample_info,tls_classification3,by.x = 'Sample_ID',by.y = 'Sample_ID',all=T)
colnames(blca_meta_combined_with_TLS)[1] = 'ST_ID'

blca_meta_combined_with_TLS$'E-TLS' = ifelse(is.na(blca_meta_combined_with_TLS$'E-TLS'),0,blca_meta_combined_with_TLS$'E-TLS')
blca_meta_combined_with_TLS$'P-TLS' = ifelse(is.na(blca_meta_combined_with_TLS$'P-TLS'),0,blca_meta_combined_with_TLS$'P-TLS')

blca_meta_combined_with_TLS$TLS_class <- ifelse(rowSums(blca_meta_combined_with_TLS[,c('E-TLS','P-TLS')]) == 0,"No_TLS",
                                                apply(blca_meta_combined_with_TLS[, c('E-TLS', 'P-TLS')], 1, function(x) {
                                                  categories <- c('E-TLS', 'P-TLS', 'Secondary')
                                                  max_value <- max(x)
                                                  max_indices <- which(x == max_value)
                                                  tied_categories <- categories[max_indices]
                                                  if ("P-TLS" %in% tied_categories) {
                                                    return("P-TLS")
                                                  }  else {
                                                    return("E-TLS")
                                                  }
                                                }))

blca_meta_combined_with_TLS$TLS_class = factor(blca_meta_combined_with_TLS$TLS_class,levels = c("No_TLS","E-TLS","P-TLS"))

blca_meta_combined_with_TLS$Responder = factor(blca_meta_combined_with_TLS$Responder,levels=c("N","Y"),labels=c("NR","R"))
blca_meta_combined_with_TLS$MTAP =  factor(blca_meta_combined_with_TLS$MTAP,levels = c("Neg","Pos"),labels=c("MTAP-","MTAP+"))

p = ggplot(blca_meta_combined_with_TLS, aes_string(axis1 = 'TLS_class', axis2 = "MTAP", fill = 'TLS_class')) +
  geom_alluvium(aes(fill = TLS_class), width = 1/6) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/6) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("TLS_class", "MTAP"), expand = c(0.15, 0.15)) +
  scale_fill_manual(values = c(
    "No_TLS" = "#f7e6ca",
    "E-TLS"  = "#71B28B",
    "P-TLS"  = "#F06825",
    "MTAP+"     = "#17365e",
    "MTAP-"    = "#cfdff2"
  )) +
  labs(title = '') +
  theme_void()

pdf("./result/Fig5Q_1.pdf")
print(p)
dev.off()


p = ggplot(blca_meta_combined_with_TLS, aes_string(axis1 = 'TLS_class', axis2 = "Responder", fill = 'TLS_class')) +
  geom_alluvium(aes(fill = TLS_class), width = 1/6) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/6) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("TLS_class", "Responder"), expand = c(0.15, 0.15)) +
  scale_fill_manual(values = c(
    "No_TLS" = "#f7e6ca",
    "E-TLS"  = "#71B28B",
    "P-TLS"  = "#F06825",
    "NR"     = "#e0f0e1",
    "R"    = "#0f8c44"
  )) +
  labs(title = '') +
  theme_void()

pdf("./result/Fig5Q_2.pdf")
print(p)
dev.off()
rm(list = ls())




