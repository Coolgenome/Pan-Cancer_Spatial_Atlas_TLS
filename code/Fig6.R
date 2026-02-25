### Figure 6 ###
### The codes are separated by Figures ###

### load packages ###
library(readxl)
library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(ggsci)
library(adjustedCurves)
library(ggplot2)
library(ggpubr)
library(stringr)
library(grid)
library(ggalluvial)
library(ggridges)
library(purrr)
library(broom)
library(DirichletReg)
library(gridExtra)
library(mclust)
library(reshape2)

### Figure 6A, 6F and 6H were created with Biorender ###
### HE TLS detection and classification pipeline is described in Fig6.sh ###

#### Figure 6B ####
TLS_classification_TCGA <- readRDS("./data/TLS_classification_TCGA.rds")
TCGA_all_summary <- TLS_classification_TCGA %>%
  group_by(Sample_ID,Cancer_type) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

TCGA_all_summary <- data.frame(TCGA_all_summary)
TCGA_all_summary <- TCGA_all_summary[, -which(colnames(TCGA_all_summary) %in% "Slide_ID")]
TLS_negative_TCGA <- readRDS("./data/TLS_negative_TCGA.rds")
TCGA_all_summary = rbind(TCGA_all_summary,TLS_negative_TCGA)

tls_group_props <- TCGA_all_summary %>%
  mutate(
    Total = Immature + Primary + Secondary,
    Prop_Immature  = ifelse(Total > 0, Immature  / Total, 0),
    Prop_Primary   = ifelse(Total > 0, Primary   / Total, 0),
    Prop_Secondary = ifelse(Total > 0, Secondary / Total, 0),
    Prop_No_TLS    = ifelse(Total == 0, 1, 0)
  ) %>%
  select(Sample_ID, Cancer_type,
         Prop_Immature, Prop_Primary, Prop_Secondary, Prop_No_TLS) %>%
  pivot_longer(
    cols = starts_with("Prop_"),
    names_to = "Maturation",
    values_to = "Proportion"
  ) %>%
  mutate(
    Maturation = recode(Maturation,
                        Prop_Immature  = "E-TLS",
                        Prop_Primary   = "P-TLS",
                        Prop_Secondary = "S-TLS",
                        Prop_No_TLS    = "No_TLS"
    )
  ) %>%
  group_by(Cancer_type, Maturation) %>%
  summarise(mean_prop = mean(Proportion, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Maturation = factor(Maturation, levels = c("No_TLS", "E-TLS", "P-TLS", "S-TLS")))

tls_group_props = data.frame(tls_group_props)

tls_group_props$Cancer_type = factor(tls_group_props$Cancer_type,levels = c("LUSC","LUAD","STAD","BLCA","COAD","KIRC"))

a <- ggplot(tls_group_props, aes(x = Cancer_type, y = mean_prop, fill = Maturation)) +
  geom_bar(stat = "identity", position = "fill") + # fill → 100% stacked
  theme_classic(base_size = 14) +
  labs(x = "", y = "Proportion", fill = "Maturation") +
  scale_fill_manual(values = c(
    "E-TLS" = "#72B28B",
    "P-TLS" = "#F26522",
    "S-TLS" = "#B3B1D8",
    "No_TLS" = "#636363"
  )) +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    legend.text = element_text(size = 8, colour = "black"),
    legend.title = element_blank(),
    axis.text.y = element_text(margin = margin(r = 0)),
    axis.text.x = element_text(margin = margin(r = 0))
  ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.25), labels = scales::percent)

pdf("./result/Fig6B.pdf")
a
dev.off()
rm(list = ls())

#### Figure 6C ####
TLS_classification_TCGA <- readRDS("./data/TLS_classification_TCGA.rds")
TCGA_all_summary <- TLS_classification_TCGA %>%
  group_by(Sample_ID,Cancer_type) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

TCGA_all_summary <- data.frame(TCGA_all_summary)
TCGA_all_summary <- TCGA_all_summary[, -which(colnames(TCGA_all_summary) %in% "Slide_ID")]
TLS_negative_TCGA <- readRDS("./data/TLS_negative_TCGA.rds")
TCGA_all_summary = rbind(TCGA_all_summary,TLS_negative_TCGA)

TCGA_meta <- read_excel(paste0("./data_meta/HE_TCGA_meta2.xlsx"), sheet = "TCGA-CDR")[, -1]
TCGA_meta2 <- read_excel("./data_meta/HE_TCGA_meta4.xlsx")

TCGA_meta <- merge(TCGA_meta, TCGA_meta2, by.x = "bcr_patient_barcode", by.y = "TCGA barcode", all.x = TRUE)
tls_counts_with_meta <- merge(TCGA_meta, TCGA_all_summary, by.x = "bcr_patient_barcode", by.y = "Sample_ID")
tls_counts_with_meta$Tumor_Stage <- ifelse(tls_counts_with_meta$ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA", "Stage IB"), "Stage I",
                                           ifelse(tls_counts_with_meta$ajcc_pathologic_tumor_stage %in% c("Stage IIA", "Stage IIB", "Stage II"), "Stage II",
                                                  ifelse(tls_counts_with_meta$ajcc_pathologic_tumor_stage %in% c("Stage III", "Stage IIIC", "Stage IIIB", "Stage IIIA"), "Stage III",
                                                         ifelse(tls_counts_with_meta$ajcc_pathologic_tumor_stage %in% c("Stage IV", "Stage IVB", "Stage IVA"), "Stage IV", NA))))
tls_counts_with_meta <- tls_counts_with_meta[!(tls_counts_with_meta$Cancer_type == "KIRC" &tls_counts_with_meta$Tumor_Stage == "Stage IV"),]
TCGA_meta3 <- read_excel("./data_meta/HE_TCGA_meta3.xls")
TCGA_meta_COAD <- merge(TCGA_meta, TCGA_meta3, by.x = "bcr_patient_barcode", by.y = "patient", all.x = TRUE)

TCGA_meta_COAD <- TCGA_meta_COAD[TCGA_meta_COAD$bcr_patient_barcode %in% tls_counts_with_meta[tls_counts_with_meta$Cancer_type == 'COAD','bcr_patient_barcode'],c('bcr_patient_barcode','MSI_status')]

tls_counts_with_meta2 <- tls_counts_with_meta %>%
  left_join(
    TCGA_meta_COAD %>% select(bcr_patient_barcode, MSI_status),
    by = "bcr_patient_barcode"
  )

tls_counts_with_meta2 <- tls_counts_with_meta2 %>%
  mutate(
    Cancer_type = ifelse(
      Cancer_type == "COAD" & !is.na(MSI_status),
      paste0(Cancer_type, "_", MSI_status),
      Cancer_type
    )
  )

tls_group_props <- tls_counts_with_meta2 %>%
  transmute(
    Sample_ID = bcr_patient_barcode,
    Tumor_Stage,Cancer_type,
    Immature, Primary, Secondary) %>%
  mutate(
    Total = Immature + Primary + Secondary,
    Prop_Immature  = ifelse(Total > 0, Immature  / Total, 0),
    Prop_Primary   = ifelse(Total > 0, Primary   / Total, 0),
    Prop_Secondary = ifelse(Total > 0, Secondary / Total, 0),
    Prop_No_TLS    = ifelse(Total == 0, 1, 0)
  ) %>%
  select(Sample_ID, Cancer_type,Tumor_Stage,
         Prop_Immature, Prop_Primary, Prop_Secondary, Prop_No_TLS) %>%
  pivot_longer(
    cols = starts_with("Prop_"),
    names_to = "Maturation",
    values_to = "Proportion"
  ) %>%
  mutate(
    Maturation = recode(Maturation,
                        Prop_Immature  = "E-TLS",
                        Prop_Primary   = "P-TLS",
                        Prop_Secondary = "S-TLS",
                        Prop_No_TLS    = "No_TLS"
    )
  ) %>%
  group_by(Cancer_type,Tumor_Stage,Maturation) %>%
  summarise(mean_prop = mean(Proportion, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Maturation = factor(Maturation, levels = c("No_TLS", "E-TLS", "P-TLS", "S-TLS")))

tls_group_props = data.frame(tls_group_props)
tls_group_props = tls_group_props[!is.na(tls_group_props$Tumor_Stage),]
tls_group_props = tls_group_props[tls_group_props$Cancer_type %in% c("STAD",'COAD_MSI-L','COAD_MSS','KIRC'),]
tls_group_props$Cancer_type = factor(tls_group_props$Cancer_type,levels = c("STAD",'COAD_MSI-L','COAD_MSS','KIRC'))

plot_df <- droplevels(subset(tls_group_props, !(Cancer_type == "KIRC" & Tumor_Stage == "Stage IV")))

a <- ggplot(plot_df, aes(x = Tumor_Stage, y = mean_prop, fill = Maturation)) +
  geom_bar(stat = "identity", position = "fill") + # fill → 100% stacked
  theme_classic(base_size = 14) +
  facet_wrap(~Cancer_type, nrow = 1, scales = "free_x") +
  labs(x = "", y = "Proportion", fill = "Maturation") +
  scale_fill_manual(values = c(
    "E-TLS" = "#72B28B",
    "P-TLS" = "#F26522",
    "S-TLS" = "#B3B1D8",
    "No_TLS" = "#636363"
  )) +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    legend.text = element_text(size = 8, colour = "black"),
    legend.title = element_blank(),
    axis.text.y = element_text(margin = margin(r = 0)),
    axis.text.x = element_text(margin = margin(r = 0))
  ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.25), labels = scales::percent)
pdf("./result/Fig6C.pdf",width=10)
a
dev.off()
rm(list = ls())

