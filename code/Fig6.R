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
library(tibble)

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

maturation_levels <- c("Immature", "Primary", "Secondary", "No_TLS")

tls_counts_with_meta3 = tls_counts_with_meta2
tls_counts_with_meta3 = tls_counts_with_meta3[,c("bcr_patient_barcode", "Tumor_Stage","Cancer_type","Immature", "Primary", "Secondary")]
tls_counts_with_meta3 <- tls_counts_with_meta3 %>%
  pivot_longer(
    cols = c(Immature, Primary, Secondary),
    names_to = "maturation",
    values_to = "count"
  )

tls_props_fixed <- tls_counts_with_meta3 %>%
  group_by(bcr_patient_barcode) %>%
  mutate(n_rows = n()) %>%
  ungroup() %>%
  group_by(bcr_patient_barcode, Tumor_Stage, Cancer_type) %>%
  summarise(
    Immature  = sum(count[maturation == "Immature"],  na.rm = TRUE),
    Primary   = sum(count[maturation == "Primary"],   na.rm = TRUE),
    Secondary = sum(count[maturation == "Secondary"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    No_TLS = if_else(Immature + Primary + Secondary > 0, 0L, 1L)
  ) %>%
  pivot_longer(
    cols = all_of(maturation_levels),
    names_to = "maturation",
    values_to = "count"
  ) %>%
  group_by(bcr_patient_barcode) %>%
  mutate(
    TLS_total = sum(count),
    prop = if_else(TLS_total > 0, count / TLS_total, 0)
  ) %>%
  ungroup()

tls_props_fixed <- tls_props_fixed %>%
  select(bcr_patient_barcode,Tumor_Stage,TLS_total,maturation,count,prop,Cancer_type)

tls_props_long_all = tls_props_fixed

run_chi_for_one_cancer <- function(df_one) {
  chi_tbl <- df_one %>%
    group_by(Tumor_Stage, maturation) %>%
    summarise(n = sum(count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from  = maturation,
      values_from = n,
      values_fill = 0
    ) %>%
    arrange(Tumor_Stage)
  
  if (nrow(chi_tbl) < 2) {
    return(list(global = NULL, pairwise = NULL, chi_tbl = chi_tbl))
  }
  
  chi_mat <- chi_tbl %>% column_to_rownames("Tumor_Stage") %>% as.matrix()
  
  global_test <- tryCatch(chisq.test(chi_mat), error = function(e) NULL)
  
  global_df <- if (is.null(global_test)) {
    NULL
  } else {
    tibble(
      p_value = global_test$p.value,
      statistic = unname(global_test$statistic),
      df = unname(global_test$parameter),
      min_expected = min(global_test$expected)
    )
  }
  
  stages <- rownames(chi_mat)
  pairs <- combn(stages, 2, simplify = FALSE)
  
  pairwise_df <- map_dfr(pairs, function(x) {
    mat2 <- chi_mat[x, , drop = FALSE]
    
    tst <- tryCatch(chisq.test(mat2), error = function(e) NULL)
    if (is.null(tst)) {
      return(tibble(stage1 = x[1], stage2 = x[2],
                    p_value = NA_real_, statistic = NA_real_, df = NA_real_,
                    min_expected = NA_real_))
    }
    
    tibble(
      stage1 = x[1],
      stage2 = x[2],
      p_value = tst$p.value,
      statistic = unname(tst$statistic),
      df = unname(tst$parameter),
      min_expected = min(tst$expected)
    )
  }) %>%
    mutate(padj_fdr = p.adjust(p_value, method = "BH"))
  
  list(global = global_df, pairwise = pairwise_df, chi_tbl = chi_tbl)
}
tls_props_long_all = tls_props_long_all[!is.na(tls_props_long_all$Tumor_Stage),]
tls_props_long_all
res <- tls_props_long_all %>%
  group_by(Cancer_type) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.x$Cancer_type)[1])) %>%
  purrr::map(run_chi_for_one_cancer)

pairwise_by_cancer <- imap_dfr(res, ~ {
  if (is.null(.x$pairwise)) return(tibble())
  mutate(.x$pairwise, cancer_type = .y, .before = 1)
})


pairwise_by_cancer_filtered = pairwise_by_cancer[pairwise_by_cancer$padj_fdr < 0.05,]
pairwise_by_cancer_filtered = pairwise_by_cancer_filtered[!is.na(pairwise_by_cancer_filtered$cancer_type),]
pairwise_by_cancer_filtered$stage_all = paste0(pairwise_by_cancer_filtered$stage1,' vs ',pairwise_by_cancer_filtered$stage2)
STAD_pvalue = pairwise_by_cancer_filtered[pairwise_by_cancer_filtered$cancer_type=='STAD',]
STAD_pvalue = STAD_pvalue[STAD_pvalue$stage_all %in% c("Stage I vs Stage III","Stage II vs Stage III"),]
COAD_MSI_L_pvalue = pairwise_by_cancer_filtered[pairwise_by_cancer_filtered$cancer_type=='COAD_MSI-L',]
COAD_MSI_L_pvalue = COAD_MSI_L_pvalue[COAD_MSI_L_pvalue$stage_all %in% c("Stage I vs Stage III","Stage I vs Stage IV"),]
COAD_MSS_pvalue = pairwise_by_cancer_filtered[pairwise_by_cancer_filtered$cancer_type=='COAD_MSS',]
COAD_MSS_pvalue = COAD_MSS_pvalue[COAD_MSS_pvalue$stage_all %in% c("Stage I vs Stage III","Stage II vs Stage IV"),]
KIRC_pvalue = pairwise_by_cancer_filtered[pairwise_by_cancer_filtered$cancer_type=='KIRC',]

pvalues = rbind(STAD_pvalue, COAD_MSI_L_pvalue, COAD_MSS_pvalue, KIRC_pvalue)
pvalues
rm(list = ls())

#### Figure 6D ####
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

tls_counts_with_meta2 = tls_counts_with_meta2[!is.na(tls_counts_with_meta2$Tumor_Stage),]
tls_counts_with_meta2 = tls_counts_with_meta2[,c("bcr_patient_barcode","Tumor_Stage","Immature","Primary","Secondary","Cancer_type")]

maturation_levels <- c("Immature", "Primary", "Secondary", "No_TLS")

tls_counts_dedup <- tls_counts_with_meta2 %>%
  mutate(TLS_total = Immature + Primary + Secondary) %>%
  group_by(bcr_patient_barcode) %>%
  arrange(desc(TLS_total), .by_group = TRUE) %>%
  slice(1) %>%                      
  ungroup() %>%
  select(-TLS_total) 

tls_props_long_all = tls_counts_dedup %>%
  mutate(
    TLS_total = Immature + Primary + Secondary
  ) %>%
  pivot_longer(
    cols = c(Immature, Primary, Secondary),
    names_to = "maturation",
    values_to = "count"
  ) %>%
  mutate(
    prop = ifelse(TLS_total > 0, count / TLS_total, 0)
  ) %>%
  bind_rows(
    tls_counts_dedup %>%
      mutate(
        TLS_total = Immature + Primary + Secondary,
        maturation = "No_TLS",
        count = 0,
        prop = ifelse(TLS_total == 0, 1, 0)
      ) %>%
      select(bcr_patient_barcode, Tumor_Stage, TLS_total, maturation, count, prop, Cancer_type)
  ) %>%
  select(bcr_patient_barcode, Tumor_Stage, TLS_total, maturation, count, prop, Cancer_type)

tls_props_long_all$maturation <- factor(tls_props_long_all$maturation, levels = c("No_TLS", "Immature", "Primary", "Secondary"))
tls_props_long_sub <- tls_props_long_all[tls_props_long_all$Cancer_type %in% c("COAD_MSS", "COAD_MSI-L", "STAD", "KIRC"), ]

tls_props_long_all2 <- tls_props_long_all[tls_props_long_all$maturation %in% c("Immature", "Primary", "Secondary"), ]
tls_props_long_sub2 <- tls_props_long_all2[tls_props_long_all2$Cancer_type %in% c("COAD_MSS", "COAD_MSI-L", "STAD", "KIRC"), ]

a <- ggplot(tls_props_long_sub2,aes(x = prop, y = Tumor_Stage, fill = maturation)) +
  geom_density_ridges(alpha = 1, scale = 1, from = 0, to = 1) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_fill_manual(values = c(
    "Immature" = "#72B28B",
    "Primary" = "#F26522",
    "Secondary" = "#B3B1D8"
  )) +
  scale_y_discrete(limits = rev) +
  theme_classic(base_size = 14) +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    legend.text = element_text(size = 8, colour = "black"),
    legend.title = element_blank(),
    axis.text.y = element_text(margin = margin(r = 0)),
    axis.text.x = element_text(margin = margin(r = 0)),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )
pdf("./result/Fig6D.pdf", width = 5, height = 8)
a
dev.off()

library(DirichletReg)

tls_props_wide_all <- tls_props_long_all %>%
  select(bcr_patient_barcode, Tumor_Stage, Cancer_type, maturation, prop) %>%
  pivot_wider(
    id_cols = c(bcr_patient_barcode, Tumor_Stage, Cancer_type),
    names_from = maturation,
    values_from = prop,
    values_fill = 0
  )


Y <- DR_data(tls_props_wide_all[, c("Immature", "Primary", "Secondary", "No_TLS")])

fit_stage <- DirichReg(
  Y ~ Tumor_Stage,
  data = tls_props_wide_all
)

fit_null <- DirichReg(
  Y ~ 1,
  data = tls_props_wide_all
)

anova(fit_null, fit_stage)

fit_stage1 <- DirichReg(
  Y ~ Tumor_Stage2,
  data = tls_props_wide_all %>%
    mutate(Tumor_Stage2 = factor(ifelse(Tumor_Stage == "Stage I", "Stage I", "Other"))) %>%
    mutate(Tumor_Stage2 = factor(Tumor_Stage2))
)

a1 <- anova(fit_null, fit_stage1)

fit_stage2 <- DirichReg(
  Y ~ Tumor_Stage2,
  data = tls_props_wide_all %>%
    mutate(Tumor_Stage2 = factor(ifelse(Tumor_Stage == "Stage II", "Stage II", "Other"))) %>%
    mutate(Tumor_Stage2 = factor(Tumor_Stage2))
)

a2 <- anova(fit_null, fit_stage2)

fit_stage3 <- DirichReg(
  Y ~ Tumor_Stage2,
  data = tls_props_wide_all %>%
    mutate(Tumor_Stage2 = factor(ifelse(Tumor_Stage == "Stage III", "Stage III", "Other"))) %>%
    mutate(Tumor_Stage2 = factor(Tumor_Stage2))
)

a3 <- anova(fit_null, fit_stage3)

fit_stage4 <- DirichReg(
  Y ~ Tumor_Stage2,
  data = tls_props_wide_all %>%
    mutate(Tumor_Stage2 = factor(ifelse(Tumor_Stage == "Stage IV", "Stage IV", "Other"))) %>%
    mutate(Tumor_Stage2 = factor(Tumor_Stage2))
)

a4 <- anova(fit_null, fit_stage4)


pvals <- c(
  stage1 = a1$`Pr(>Chi)`[2],
  stage2 = a2$`Pr(>Chi)`[2],
  stage3 = a3$`Pr(>Chi)`[2],
  stage4 = a4$`Pr(>Chi)`[2]
)

p.adjust(pvals, method = "fdr")

p_i <- kruskal.test(Immature ~ Tumor_Stage, data = tls_props_wide_all)
p_p <- kruskal.test(Primary ~ Tumor_Stage, data = tls_props_wide_all)
p_s <- kruskal.test(Secondary ~ Tumor_Stage, data = tls_props_wide_all)

p.adjust(
  c(
    Immature  = p_i$p.value,
    Primary   = p_p$p.value,
    Secondary = p_s$p.value
  ),
  method = "fdr"
)

tls_props_wide_all2 <- tls_props_long_all2 %>%
  select(bcr_patient_barcode, Tumor_Stage, Cancer_type, maturation, prop) %>%
  pivot_wider(
    id_cols = c(bcr_patient_barcode, Tumor_Stage, Cancer_type),
    names_from = maturation,
    values_from = prop,
    values_fill = 0
  )


Y <- DR_data(tls_props_wide_all2[, c("Immature", "Primary", "Secondary")])

fit_stage <- DirichReg(
  Y ~ Tumor_Stage,
  data = tls_props_wide_all2
)

fit_null <- DirichReg(
  Y ~ 1,
  data = tls_props_wide_all2
)

anova(fit_null, fit_stage)

fit_stage1 <- DirichReg(
  Y ~ Tumor_Stage2,
  data = tls_props_wide_all2 %>%
    mutate(Tumor_Stage2 = factor(ifelse(Tumor_Stage == "Stage I", "Stage I", "Other"))) %>%
    mutate(Tumor_Stage2 = factor(Tumor_Stage2))
)

a1 <- anova(fit_null, fit_stage1)

fit_stage2 <- DirichReg(
  Y ~ Tumor_Stage2,
  data = tls_props_wide_all2 %>%
    mutate(Tumor_Stage2 = factor(ifelse(Tumor_Stage == "Stage II", "Stage II", "Other"))) %>%
    mutate(Tumor_Stage2 = factor(Tumor_Stage2))
)

a2 <- anova(fit_null, fit_stage2)

fit_stage3 <- DirichReg(
  Y ~ Tumor_Stage2,
  data = tls_props_wide_all2 %>%
    mutate(Tumor_Stage2 = factor(ifelse(Tumor_Stage == "Stage III", "Stage III", "Other"))) %>%
    mutate(Tumor_Stage2 = factor(Tumor_Stage2))
)

a3 <- anova(fit_null, fit_stage3)

fit_stage4 <- DirichReg(
  Y ~ Tumor_Stage2,
  data = tls_props_wide_all2 %>%
    mutate(Tumor_Stage2 = factor(ifelse(Tumor_Stage == "Stage IV", "Stage IV", "Other"))) %>%
    mutate(Tumor_Stage2 = factor(Tumor_Stage2))
)

a4 <- anova(fit_null, fit_stage4)


pvals <- c(
  stage1 = a1$`Pr(>Chi)`[2],
  stage2 = a2$`Pr(>Chi)`[2],
  stage3 = a3$`Pr(>Chi)`[2],
  stage4 = a4$`Pr(>Chi)`[2]
)

p.adjust(pvals, method = "fdr")

p_i <- kruskal.test(Immature ~ Tumor_Stage, data = tls_props_wide_all2)
p_p <- kruskal.test(Primary ~ Tumor_Stage, data = tls_props_wide_all2)
p_s <- kruskal.test(Secondary ~ Tumor_Stage, data = tls_props_wide_all2)

p.adjust(
  c(
    Immature  = p_i$p.value,
    Primary   = p_p$p.value,
    Secondary = p_s$p.value
  ),
  method = "fdr"
)


tls_props_wide_sub <- tls_props_long_sub %>%
  select(bcr_patient_barcode, Tumor_Stage, Cancer_type, maturation, prop) %>%
  pivot_wider(
    id_cols = c(bcr_patient_barcode, Tumor_Stage, Cancer_type),
    names_from = maturation,
    values_from = prop,
    values_fill = 0
  )


Y <- DR_data(tls_props_wide_sub[, c("Immature", "Primary", "Secondary", "No_TLS")])

fit_stage <- DirichReg(
  Y ~ Tumor_Stage,
  data = tls_props_wide_sub
)

fit_null <- DirichReg(
  Y ~ 1,
  data = tls_props_wide_sub
)

anova(fit_null, fit_stage)

fit_stage1 <- DirichReg(
  Y ~ Tumor_Stage2,
  data = tls_props_wide_sub %>%
    mutate(Tumor_Stage2 = factor(ifelse(Tumor_Stage == "Stage I", "Stage I", "Other"))) %>%
    mutate(Tumor_Stage2 = factor(Tumor_Stage2))
)

a1 <- anova(fit_null, fit_stage1)

fit_stage2 <- DirichReg(
  Y ~ Tumor_Stage2,
  data = tls_props_wide_sub %>%
    mutate(Tumor_Stage2 = factor(ifelse(Tumor_Stage == "Stage II", "Stage II", "Other"))) %>%
    mutate(Tumor_Stage2 = factor(Tumor_Stage2))
)

a2 <- anova(fit_null, fit_stage2)

fit_stage3 <- DirichReg(
  Y ~ Tumor_Stage2,
  data = tls_props_wide_sub %>%
    mutate(Tumor_Stage2 = factor(ifelse(Tumor_Stage == "Stage III", "Stage III", "Other"))) %>%
    mutate(Tumor_Stage2 = factor(Tumor_Stage2))
)

a3 <- anova(fit_null, fit_stage3)

fit_stage4 <- DirichReg(
  Y ~ Tumor_Stage2,
  data = tls_props_wide_sub %>%
    mutate(Tumor_Stage2 = factor(ifelse(Tumor_Stage == "Stage IV", "Stage IV", "Other"))) %>%
    mutate(Tumor_Stage2 = factor(Tumor_Stage2))
)

a4 <- anova(fit_null, fit_stage4)


pvals <- c(
  stage1 = a1$`Pr(>Chi)`[2],
  stage2 = a2$`Pr(>Chi)`[2],
  stage3 = a3$`Pr(>Chi)`[2],
  stage4 = a4$`Pr(>Chi)`[2]
)

p.adjust(pvals, method = "fdr")

rm(list = ls())

#### Figure 6E ####
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

tls_counts_BLCA = tls_counts_with_meta2[tls_counts_with_meta2$Cancer_type == 'BLCA',]
tls_counts_BLCA = tls_counts_BLCA[tls_counts_BLCA$Tumor_Stage != 'Stage I',]
tls_counts_BLCA = tls_counts_BLCA[!is.na(tls_counts_BLCA$OS.time),]
tls_counts_BLCA$OS.Time.months <- tls_counts_BLCA$OS.time / 30.44

tls_counts_BLCA$total_TLS <- with(tls_counts_BLCA, Immature + Primary + Secondary)
tls_counts_BLCA$ratio_Immature <- ifelse(tls_counts_BLCA$total_TLS == 0, NA, tls_counts_BLCA$Immature / tls_counts_BLCA$total_TLS)
tls_counts_BLCA$ratio_Primary <- ifelse(tls_counts_BLCA$total_TLS == 0, NA, tls_counts_BLCA$Primary / tls_counts_BLCA$total_TLS)
tls_counts_BLCA$ratio_Secondary <- ifelse(tls_counts_BLCA$total_TLS == 0, NA, tls_counts_BLCA$Secondary / tls_counts_BLCA$total_TLS)

tls_counts_BLCA$TLS_bin <- ifelse(tls_counts_BLCA$total_TLS == 0, 1, 0)
tls_counts_BLCA$TLS_status <- ifelse(tls_counts_BLCA$total_TLS == 0, "Without TLS", "With TLS")
tls_counts_BLCA$TLS_status <- factor(tls_counts_BLCA$TLS_status, levels = c("With TLS", "Without TLS"))

tls_counts_BLCA$TLS_most_common <- ifelse(
  rowSums(tls_counts_BLCA[, c("Immature", "Primary", "Secondary")] == 0) == 3,
  "No_TLS",
  apply(tls_counts_BLCA[, c("Immature", "Primary", "Secondary")], 1, function(x) {
    categories <- c("Immature", "Primary", "Secondary")
    categories[which.max(x)]
  })
)
tls_counts_BLCA$TLS_most_common <- factor(tls_counts_BLCA$TLS_most_common, levels = c("No_TLS", "Immature", "Primary", "Secondary"), labels = c("No TLS", "E-TLS dom", "P-TLS dom", "S-TLS dom"))

tls_counts_BLCA$ETLS_status <- ifelse(tls_counts_BLCA$total_TLS == 0, "No TLS",
                                      ifelse(tls_counts_BLCA$Immature > 0 & tls_counts_BLCA$Primary == 0 & tls_counts_BLCA$Secondary == 0, "E-TLS", "mTLS")
)
tls_counts_BLCA$ETLS_status <- factor(tls_counts_BLCA$ETLS_status, levels = c("No TLS", "E-TLS", "mTLS"))

fit_tls1 <- survfit(Surv(OS.Time.months, OS) ~ TLS_status, data = tls_counts_BLCA)

plot_tls1 <- ggsurvplot(fit_tls1,
                        data = tls_counts_BLCA, title = "BLCA",
                        pval = TRUE, risk.table = FALSE, xlab = "Time (months)", palette = c("#68a0d6", "#de6465"))

pdf("./result/Fig6E_1.pdf", width = 12, height = 8)
plot_tls1$plot
dev.off()

fit_tls2 <- survfit(Surv(OS.Time.months, OS) ~ TLS_most_common, data = tls_counts_BLCA)

plot_tls2 <- ggsurvplot(fit_tls2,
                        data = tls_counts_BLCA, title = "BLCA",
                        pval = TRUE, risk.table = FALSE, xlab = "Time (months)", palette = c("#646666", "#6eb28b", "#ed6825", "#7b799a"))

pdf("./result/Fig6E_2.pdf", width = 12, height = 8)
plot_tls2$plot
dev.off()

fit_etls <- survfit(Surv(OS.Time.months, OS) ~ ETLS_status, data = tls_counts_BLCA)

plot_etls <- ggsurvplot(fit_etls,
                        data = tls_counts_BLCA, title = "BLCA",
                        pval = TRUE, risk.table = FALSE, xlab = "Time (months)", palette = c("#646666", "#6eb28b", "#7092cb"))

pdf("./result/Fig6E_3.pdf", width = 12, height = 8)
plot_etls$plot
dev.off()

rm(list = ls())

#### Figure 6G ####
tls_df = read_excel("./data/Table7.xlsx")
colnames(tls_df)[1:2] = c("Sample_ID","Cancer_type")
TCGA_meta <- read_excel(paste0("./data_meta/HE_TCGA_meta2.xlsx"), sheet = "TCGA-CDR")[, -1]
TCGA_meta$OS.Time.months = TCGA_meta$OS.time / 30.44
TCGA_meta = TCGA_meta[,c('bcr_patient_barcode','OS','OS.Time.months')]
df_meta = merge(tls_df,TCGA_meta,by.x='Sample_ID',by.y='bcr_patient_barcode')

df_meta_BLCA = df_meta[df_meta$Cancer_type == 'BLCA',]
df_meta_LUSC = df_meta[df_meta$Cancer_type == 'LUSC',]
df_meta_LUAD = df_meta[df_meta$Cancer_type == 'LUAD',]
df_meta_STAD = df_meta[df_meta$Cancer_type == 'STAD',]
df_meta_COAD = df_meta[df_meta$Cancer_type == 'COAD',]
df_meta_KIRC = df_meta[df_meta$Cancer_type == 'KIRC',]

fit_BLCA <- survfit(Surv(OS.Time.months, OS) ~ TLS_Group, data = df_meta_BLCA)
fit_LUSC <- survfit(Surv(OS.Time.months, OS) ~ TLS_Group, data = df_meta_LUSC)
fit_LUAD <- survfit(Surv(OS.Time.months, OS) ~ TLS_Group, data = df_meta_LUAD)
fit_STAD <- survfit(Surv(OS.Time.months, OS) ~ TLS_Group, data = df_meta_STAD)
fit_COAD <- survfit(Surv(OS.Time.months, OS) ~ TLS_Group, data = df_meta_COAD)
fit_KIRC <- survfit(Surv(OS.Time.months, OS) ~ TLS_Group, data = df_meta_KIRC)

plot_BLCA <- ggsurvplot(fit_BLCA, data = df_meta_BLCA, title = "BLCA", risk.table = TRUE, pval = TRUE, palette = c("#7B799A", "#70B28B"))
plot_LUSC <- ggsurvplot(fit_LUSC, data = df_meta_LUSC, title = "LUSC", risk.table = TRUE, pval = TRUE, palette = c("#7B799A", "#70B28B"))
plot_LUAD <- ggsurvplot(fit_LUAD, data = df_meta_LUAD, title = "LUAD", risk.table = TRUE, pval = TRUE, palette = c("#7B799A", "#70B28B"))
plot_STAD <- ggsurvplot(fit_STAD, data = df_meta_STAD, title = "STAD", risk.table = TRUE, pval = TRUE, palette = (c("#7B799A", "#70B28B")))
plot_COAD <- ggsurvplot(fit_COAD, data = df_meta_COAD, title = "COAD", risk.table = TRUE, pval = TRUE, palette = c("#7B799A", "#70B28B"))
plot_KIRC <- ggsurvplot(fit_KIRC, data = df_meta_KIRC, title = "KIRC", risk.table = TRUE, pval = TRUE, palette = c("#7B799A", "#70B28B"))

p_BLCA <- plot_BLCA$plot
p_LUSC <- plot_LUSC$plot
p_LUAD <- plot_LUAD$plot + scale_x_continuous(breaks = c(0, 60, 120, 180, 240))
p_STAD <- plot_STAD$plot
p_COAD <- plot_COAD$plot
p_KIRC <- plot_KIRC$plot

t_BLCA <- plot_BLCA$table
t_LUSC <- plot_LUSC$table
t_LUAD <- plot_LUAD$table + scale_x_continuous(breaks = c(0, 60, 120, 180, 240))
t_STAD <- plot_STAD$table
t_COAD <- plot_COAD$table
t_KIRC <- plot_KIRC$table

pdf("./result/Fig6G_1.pdf", width = 12, height = 8)
print(p_BLCA)
dev.off()

pdf("./result/Fig6G_2.pdf", width = 12, height = 8)
print(p_LUSC)
dev.off()

pdf("./result/Fig6G_3.pdf", width = 12, height = 8)
print(p_LUAD)
dev.off()

pdf("./result/Fig6G_4.pdf", width = 12, height = 8)
print(p_STAD)
dev.off()

pdf("./result/Fig6G_5.pdf", width = 12, height = 8)
print(p_COAD)
dev.off()

pdf("./result/Fig6G_6.pdf", width = 12, height = 8)
print(p_KIRC)
dev.off()

pdf("./result/Fig6G_1_legend.pdf", width = 12, height = 8)
print(t_BLCA)
dev.off()

pdf("./result/Fig6G_2_legend.pdf", width = 12, height = 8)
print(t_LUSC)
dev.off()

pdf("./result/Fig6G_3_legend.pdf", width = 12, height = 8)
print(t_LUAD)
dev.off()

pdf("./result/Fig6G_4_legend.pdf", width = 12, height = 8)
print(t_STAD)
dev.off()

pdf("./result/Fig6G_5_legend.pdf", width = 12, height = 8)
print(t_COAD)
dev.off()

pdf("./result/Fig6G_6_legend.pdf", width = 12, height = 8)
print(t_KIRC)
dev.off()

rm(list = ls())

#### Figure 6I ####
occc_df = readRDS("./data/TLS_classification_OCCC.rds")
colnames(occc_df)[which(colnames(occc_df)=='Immature')] = 'E-TLS'
colnames(occc_df)[which(colnames(occc_df)=='Primary')] = 'P-TLS'
colnames(occc_df)[which(colnames(occc_df)=='Secondary')] = 'S-TLS'

occc_survival <- read_excel("./data_meta/HE_OCCC_meta.xlsx")
occc_survival <- occc_survival[occc_survival$`Patient ID` %in% occc_df$Patient_ID, ]
occc_survival$`Patient ID` <- as.character(occc_survival$`Patient ID`)
occc_df_survival <- merge(occc_df, occc_survival, by.x = "Patient_ID", by.y = "Patient ID", all.y = TRUE)

fit <- survfit(Surv(`OS time (month)`, `OS stauts`) ~ 1, data = occc_df_survival)
t_med <- unname(summary(fit)$table["median"])
t_75 <- unname(quantile(fit, probs = c(0.75))$quantile)

occc_df_survival <- occc_df_survival %>%
  mutate(
    km_75 = case_when(
      `OS stauts` == 1 & `OS time (month)` < t_75 ~ "Short",
      `OS time (month)` >= t_med ~ "Long",
      TRUE ~ NA_character_
    )
  )

occc_df <- merge(occc_df, occc_df_survival)

occc_df <- occc_df %>%
  rowwise() %>%
  mutate(total_TLS = sum(c(`E-TLS`, `P-TLS`, `S-TLS`))) %>%
  mutate(
    prop_E = ifelse(total_TLS > 0, `E-TLS` / total_TLS, 0),
    prop_P = ifelse(total_TLS > 0, `P-TLS` / total_TLS, 0),
    prop_S = ifelse(total_TLS > 0, `S-TLS` / total_TLS, 0),
    No_TLS = ifelse(total_TLS == 0, 1, 0)
  ) %>%
  ungroup()
occc_df <- occc_df %>% rename(E_TLS = `E-TLS`, P_TLS = `P-TLS`, S_TLS = `S-TLS`)
occc_df_pre <- occc_df[occc_df$time_point == "pre", ]

pal_all <- c(
  "Short" = "#bc83a1",
  "Long" = "#a1bc83",
  "WT" = "#bcba83",
  "Mut" = "#8385bc"
)
occc_df_pre$`PPP2R1A mutation` <- ifelse(occc_df_pre$`PPP2R1A mutation` == 0, "WT", "Mut")
for (i in c("PPP2R1A mutation", "km_75")) {
  for (j in c("E_TLS", "P_TLS", "S_TLS")) {
    p1 <- ggviolin(occc_df_pre,
                   x = i, y = j, color = i, fill = i, palette = pal_all,
                   trim = TRUE, width = 1,
                   outlier.shape = NA
    ) +
      stat_summary(fun = mean, geom = "point", size = 1, color = "black") +
      stat_compare_means(method = "wilcox.test", size = 2, label.y = max(occc_df_pre[j]) - 0.5) +
      labs(title = j, y = "Counts") +
      theme(
        text = element_text(size = 8, colour = "black"),
        axis.text = element_text(size = 8, colour = "black"),
        axis.title = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 8, colour = "black"),
        plot.subtitle = element_text(size = 8, colour = "black"),
        legend.position = "none"
      )
    
    g <- ggplotGrob(p1)
    panel_idx <- g$layout[g$layout$name == "panel", ]
    g$widths[panel_idx$l] <- unit(100, "pt") 
    g$heights[panel_idx$t] <- unit(120, "pt") 
    
    w_in <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)
    h_in <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)
    title1 <- ifelse(i == "PPP2R1A mutation", "Fig6I_1", "Fig6I_2")
    title2 <- ifelse(j == "E_TLS", "_1",
                     ifelse(j == "P_TLS", "_2", "_3")
    )
    
    pdf(paste0("./result/",title1, title2, ".pdf"),
        width = w_in, height = h_in, useDingbats = FALSE
    )
    grid.newpage()
    grid.draw(g)
    dev.off()
  }
}

occc_df_long <- occc_df %>%
  pivot_longer(
    cols = c(prop_E, prop_P, prop_S, No_TLS),
    names_to = "maturation",
    values_to = "proportion"
  ) %>%
  mutate(
    maturation = factor(
      recode(
        maturation,
        prop_E = "Immature",
        prop_P = "Primary",
        prop_S = "Secondary",
        No_TLS = "No_TLS"
      ),
      levels = c("No_TLS", "Immature", "Primary", "Secondary")
    )
  )

for (i in c("PPP2R1A mutation", "km_75")) {
  occc_df_long2 <- occc_df_long %>%
    group_by(.data[[i]], time_point, maturation) %>%
    summarise(proportion = mean(proportion, na.rm = TRUE), .groups = "drop")
  
  occc_df_long2$var <- paste0(occc_df_long2[[i]], "_", occc_df_long2$time_point)
  if (i == "PPP2R1A mutation") {
    occc_df_long2$var <- factor(occc_df_long2$var, levels = c("0_pre", "0_post", "1_pre", "1_post"), labels = c("WT pre", "WT on", "Mut pre", "Mut on"))
  }
  if (i == "km_75") {
    occc_df_long2$var <- factor(occc_df_long2$var, levels = c("Short_pre", "Short_post", "Long_pre", "Long_post"), labels = c("Short pre", "Short on", "Long pre", "Long on"))
  }
  
  a <- ggplot(occc_df_long2, aes(x = var, y = proportion, fill = maturation)) +
    geom_bar(stat = "identity", position = "fill", width = 0.6) + # fill → 100% stacked
    theme_classic(base_size = 6) +
    labs(x = "", y = "Proportion", fill = "Maturation") +
    scale_fill_manual(values = c(
      "Immature" = "#72B28B",
      "Primary" = "#F26522",
      "Secondary" = "#B3B1D8",
      "No_TLS" = "#636363"
    )) +
    theme(
      text = element_text(size = 8, colour = "black"),
      axis.text = element_text(size = 8, colour = "black"),
      axis.title = element_text(size = 8, colour = "black"),
      plot.title = element_text(size = 8, colour = "black"),
      plot.subtitle = element_text(size = 8, colour = "black"),
      legend.position = "none"
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
  
  g <- ggplotGrob(a)
  panel_idx <- g$layout[g$layout$name == "panel", ]
  g$widths[panel_idx$l] <- unit(200, "pt") 
  g$heights[panel_idx$t] <- unit(400, "pt")
  w_in <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)
  h_in <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  
  title <- ifelse(i == "PPP2R1A mutation", "Fig6I_1_4", "Fig6I_2_4")
  pdf(paste0("./result/",title, ".pdf"),
      width = w_in, height = h_in, useDingbats = FALSE
  )
  grid.newpage()
  grid.draw(g)
  dev.off()
}


for (i in c("PPP2R1A mutation", "km_75")) {
  df <- occc_df
  df$var <- paste0(df[[i]], "_", df$time_point)
  df2 <- df[, c("var", "No_TLS", "E_TLS", "P_TLS", "S_TLS")]
  df_long <- df2 %>%
    group_by(var) %>%
    summarise(
      No_TLS = sum(No_TLS, na.rm = TRUE),
      E_TLS = sum(E_TLS, na.rm = TRUE),
      P_TLS = sum(P_TLS, na.rm = TRUE),
      S_TLS = sum(S_TLS, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_long <- df_long %>%
    pivot_longer(
      cols = c(No_TLS, E_TLS, P_TLS, S_TLS),
      names_to = "maturation",
      values_to = "counts"
    )
  
  tab <- df_long %>%
    tidyr::pivot_wider(
      names_from = maturation,
      values_from = counts,
      values_fill = 0
    ) %>%
    tibble::column_to_rownames("var") %>%
    as.matrix()
  tab_pre <- tab[grepl("pre$", rownames(tab)), , drop = FALSE]
  tab_post <- tab[grepl("post$", rownames(tab)), , drop = FALSE]
  tab_pre <- tab_pre[, colSums(tab_pre) > 0, drop = FALSE]
  tab_post <- tab_post[, colSums(tab_post) > 0, drop = FALSE]
  
  print(chisq.test(tab_pre))
  print(chisq.test(tab_post))
}

rm(list = ls())

### Figure 6J,K ####
tls_counts = readRDS("./data/TLS_classification_HGSOC.rds")
response <- read_excel("./data_meta/HE_HGSOC_meta.xlsx")

response$sample_ID <- gsub(".svs", "", response$`Image No.`)

tls_counts_with_meta <- merge(tls_counts, response, by.x = "sample_id", by.y = "sample_ID")
tls_counts_with_meta$response <- tls_counts_with_meta$`Treatment effect`
tls_long <- tls_counts_with_meta %>%
  pivot_longer(
    cols = c(Immature, Primary, Secondary),
    names_to = "maturation",
    values_to = "count"
  ) %>%
  mutate(
    maturation = factor(maturation,
                        levels = c("Immature", "Primary", "Secondary")
    )
  )

tls_props <- tls_counts_with_meta %>%
  mutate(
    TLS_total = Immature + Primary + Secondary
  ) %>%
  filter(TLS_total > 0)


tls_props_long <- tls_props %>%
  pivot_longer(
    cols = c(Immature, Primary, Secondary),
    names_to = "maturation",
    values_to = "count"
  ) %>%
  mutate(
    maturation = factor(maturation,
                        levels = c("Immature", "Primary", "Secondary")
    ),
    prop = count / TLS_total
  )

tls_props_long <- tls_props %>%
  pivot_longer(
    cols = c(Immature, Primary, Secondary),
    names_to = "maturation",
    values_to = "count"
  ) %>%
  mutate(
    maturation = factor(maturation,
                        levels = c("Immature", "Primary", "Secondary")
    ),
    prop = count / TLS_total
  )

tls_group_props <- tls_props_long %>%
  group_by(response, maturation) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop")

tls_group_props$response <- factor(tls_group_props$response, levels = c("invalid", "effective"), labels = c("NR", "R"))

a <- ggplot(tls_group_props, aes(x = response, y = mean_prop, fill = maturation)) +
  geom_bar(stat = "identity", position = "fill") + 
  theme_classic(base_size = 6) +
  labs(x = "effective group", y = "Proportion", fill = "Maturation") +
  scale_fill_manual(values = c(
    "Immature" = "#72B28B",
    "Primary" = "#F26522",
    "Secondary" = "#B3B1D8"
  )) +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    legend.text = element_text(size = 8, colour = "black"),
    legend.title = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_flip()


g <- ggplotGrob(a)
panel_idx <- g$layout[g$layout$name == "panel", ]
g$widths[panel_idx$l] <- unit(80, "pt")
g$heights[panel_idx$t] <- unit(72, "pt") 

w_in <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)
h_in <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)

pdf("./result/Fig6K_1.pdf",
    width = w_in, height = h_in, useDingbats = FALSE
)
grid.newpage()
grid.draw(g)
dev.off()

patient_table=read.csv("./data_meta/HE_HGSOC_meta2.csv")
tls_counts_with_meta2 <- merge(tls_counts_with_meta, patient_table, by = "sample_id")

c1_c2_obs <- as.matrix(table(tls_counts_with_meta2$TLS_Group, tls_counts_with_meta2$response))

c1_c2_chi_res <- chisq.test(c1_c2_obs)
c1_c2_expected <- c1_c2_chi_res$expected

c1_c2_ro_e <- c1_c2_obs / c1_c2_expected

rownames(c1_c2_ro_e) <- c("C1", "C2")

c1_c2_ro_e_melted <- melt(c1_c2_ro_e)
colnames(c1_c2_ro_e_melted) <- c("Group", "Treatment", "RoE")
c1_c2_ro_e_melted$Group <- factor(c1_c2_ro_e_melted$Group, levels = rev(c("C2", "C1")))
c1_c2_ro_e_melted$Treatment <- factor(c1_c2_ro_e_melted$Treatment, levels = c("invalid", "effective"), labels = c("NR", "R"))

p <- ggplot(c1_c2_ro_e_melted, aes(x = Treatment, y = Group, fill = RoE)) +
  geom_tile() +
  scale_fill_gradient2(low = "#6BAED6", mid = "white", high = "#FB6A6A", midpoint = 1) +
  labs(x = "Treatment", y = "TLS Group", fill = "Ro/e") +
  theme_minimal() +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 9, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    legend.position = "none"
  ) +
  coord_flip()

ggsave("./result/Fig6K_2.pdf", p, width = 5, height = 4)

vmin <- floor(min(c1_c2_ro_e_melted$RoE, na.rm = TRUE) * 10) / 10
vmax <- ceiling(max(c1_c2_ro_e_melted$RoE, na.rm = TRUE) * 10) / 10

p_leg <- ggplot(data.frame(x = 1, y = 1, z = 0), aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#6BAED6",
    mid = "white",
    high = "#FB6A6A",
    midpoint = 1,
    limits = c(vmin, vmax),
    breaks = seq(vmin, vmax, by = 0.1),
    oob = scales::squish,
    name = "Ro/e"
  ) +
  guides(colour = "none") + 
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )

leg <- cowplot::get_legend(p_leg) 

rot_leg <- grid::editGrob(
  leg,
  vp = grid::viewport(angle = 270)
)

grid::grid.newpage()
grid::grid.draw(rot_leg)

pdf("./result/Fig6K_2_legend.pdf", width = 2, height = 3)
grid::grid.draw(rot_leg)
dev.off()

tls_counts_with_meta$response2 <- tls_counts_with_meta$response
tls_counts_with_meta$response2 <- ifelse(tls_counts_with_meta$response2 == "invalid", "NR", "R")
tls_counts_with_meta$response2 <- factor(tls_counts_with_meta$response2, levels = c("NR", "R"))

p1 <- ggboxplot(tls_counts_with_meta,
                x = "response2", y = "Immature",
                color = "response2",
                palette = c(
                  "R"  = "#a1bc83",
                  "NR" = "#bc83a1"
                ),
                outlier.shape = NA
) + 
  geom_jitter(aes(color = response2),
              width = 0.2, size = 0.5, alpha = 0.5
  ) +
  stat_compare_means(method = "wilcox.test", label.y = 20, size = 2) +
  labs(title = "E-TLS", y = "Counts") +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 25),
    breaks = seq(0, 25, 5),
    oob = scales::squish
  )

g <- ggplotGrob(p1)
panel_idx <- g$layout[g$layout$name == "panel", ]
g$widths[panel_idx$l] <- unit(100, "pt") 
g$heights[panel_idx$t] <- unit(120, "pt") 

w_in <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)
h_in <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)

pdf("./result/Fig6J_1.pdf",
    width = w_in, height = h_in, useDingbats = FALSE
)
grid.newpage()
grid.draw(g)
dev.off()

p2 <- ggboxplot(tls_counts_with_meta,
                x = "response2", y = "Primary",
                color = "response2",
                palette = c(
                  "R"  = "#a1bc83",
                  "NR" = "#bc83a1"
                ),
                outlier.shape = NA
) + 
  geom_jitter(aes(color = response2),
              width = 0.2, size = 0.5, alpha = 0.5
  ) +
  stat_compare_means(method = "wilcox.test", label.y = 13, size = 2) +
  labs(title = "P-TLS", y = "Counts") +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 15),
    breaks = seq(0, 15, 5),
    oob = scales::squish
  )

g <- ggplotGrob(p2)
panel_idx <- g$layout[g$layout$name == "panel", ]
g$widths[panel_idx$l] <- unit(100, "pt") 
g$heights[panel_idx$t] <- unit(120, "pt") 

w_in <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)
h_in <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)

pdf("./result/Fig6J_2.pdf",
    width = w_in, height = h_in, useDingbats = FALSE
)
grid.newpage()
grid.draw(g)
dev.off()

p3 <- ggboxplot(tls_counts_with_meta,
                x = "response2", y = "Secondary",
                color = "response2",
                palette = c(
                  "R"  = "#a1bc83",
                  "NR" = "#bc83a1"
                ),
                outlier.shape = NA
) + 
  geom_jitter(aes(color = response2),
              width = 0.2, size = 0.5, alpha = 0.5
  ) +
  stat_compare_means(method = "wilcox.test", label.y = 3, size = 2) +
  labs(title = "S-TLS", y = "Counts") +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 3),
    breaks = seq(0, 3, 1),
    oob = scales::squish
  )


g <- ggplotGrob(p3)
panel_idx <- g$layout[g$layout$name == "panel", ]
g$widths[panel_idx$l] <- unit(100, "pt") 
g$heights[panel_idx$t] <- unit(120, "pt") 

w_in <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)
h_in <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)

pdf("./result/Fig6J_3.pdf",
    width = w_in, height = h_in, useDingbats = FALSE
)
grid.newpage()
grid.draw(g)
dev.off()

df_long <- tls_counts_with_meta %>%
  group_by(response) %>%
  summarise(
    E_TLS = sum(Immature, na.rm = TRUE),
    P_TLS = sum(Primary, na.rm = TRUE),
    S_TLS = sum(Secondary, na.rm = TRUE),
    .groups = "drop"
  )

df_long <- df_long %>%
  pivot_longer(
    cols = c(E_TLS, P_TLS, S_TLS),
    names_to = "maturation",
    values_to = "counts"
  )

tab <- df_long %>%
  tidyr::pivot_wider(
    names_from = maturation,
    values_from = counts,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("response") %>%
  as.matrix()

chisq.test(tab)

rm(list = ls())

#### Figure 6L,M ####
tls_counts = readRDS("./data/TLS_classification_BRCA_Sammut.rds")
breast_neoadjuvant_meta <- read.csv("./data_meta/HE_BRCA_Sammut_meta.csv")
breast_neoadjuvant_meta <- breast_neoadjuvant_meta[breast_neoadjuvant_meta$Slide.ID %in% tls_counts$sample_id, ]
breast_neoadjuvant_meta <- breast_neoadjuvant_meta[, -which(colnames(breast_neoadjuvant_meta) %in% c("Filename", "Cohort", "MD5", "ER.Allred", "RCB.score", "RCB.category"))]
tls_counts_with_meta <- merge(tls_counts, breast_neoadjuvant_meta, by.x = "sample_id", by.y = "Slide.ID")

tls_counts_with_meta <- tls_counts_with_meta[!is.na(tls_counts_with_meta$pCR.RD), ]

tls_long <- tls_counts_with_meta %>%
  pivot_longer(
    cols = c(Immature, Primary, Secondary),
    names_to = "maturation",
    values_to = "count"
  ) %>%
  mutate(
    maturation = factor(maturation,
                        levels = c("Immature", "Primary", "Secondary")
    )
  )

tls_long$Grade.pre.NAT <- factor(tls_long$Grade.pre.NAT)

tls_props <- tls_counts_with_meta %>%
  mutate(
    TLS_total = Immature + Primary + Secondary
  ) %>%
  filter(TLS_total > 0)

tls_props_long <- tls_props %>%
  pivot_longer(
    cols = c(Immature, Primary, Secondary),
    names_to = "maturation",
    values_to = "count"
  ) %>%
  mutate(
    maturation = factor(maturation,
                        levels = c("Immature", "Primary", "Secondary")
    ),
    prop = count / TLS_total
  )

tls_props_long <- tls_props %>%
  pivot_longer(
    cols = c(Immature, Primary, Secondary),
    names_to = "maturation",
    values_to = "count"
  ) %>%
  mutate(
    maturation = factor(maturation,
                        levels = c("Immature", "Primary", "Secondary")
    ),
    prop = count / TLS_total
  )

tls_group_props <- tls_props_long %>%
  group_by(pCR.RD, maturation) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop")
tls_group_props$pCR.RD <- factor(tls_group_props$pCR.RD, levels = c("RD", "pCR"), labels = c("NR", "R"))

a <- ggplot(tls_group_props, aes(x = pCR.RD, y = mean_prop, fill = maturation)) +
  geom_bar(stat = "identity", position = "fill") + 
  theme_classic(base_size = 6) +
  labs(x = "Response group", y = "Proportion", fill = "Maturation") +
  scale_fill_manual(values = c(
    "Immature" = "#72B28B",
    "Primary" = "#F26522",
    "Secondary" = "#B3B1D8"
  )) +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    legend.text = element_text(size = 8, colour = "black"),
    legend.title = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_flip()

g <- ggplotGrob(a)
panel_idx <- g$layout[g$layout$name == "panel", ]
g$widths[panel_idx$l] <- unit(80, "pt") 
g$heights[panel_idx$t] <- unit(72, "pt") 

w_in <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)
h_in <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)

pdf("./result/Fig6M_1.pdf",
    width = w_in, height = h_in, useDingbats = FALSE
)
grid.newpage()
grid.draw(g)
dev.off()

patient_table=read.csv("./data_meta/HE_BRCA_Sammut_meta2.csv")
tls_counts_with_meta2 <- merge(tls_counts_with_meta, patient_table, by = "sample_id")

c1_c2_obs <- as.matrix(table(tls_counts_with_meta2$TLS_Group, tls_counts_with_meta2$pCR.RD))

c1_c2_chi_res <- chisq.test(c1_c2_obs)
c1_c2_expected <- c1_c2_chi_res$expected

c1_c2_ro_e <- c1_c2_obs / c1_c2_expected

rownames(c1_c2_ro_e) <- c("C1", "C2")

c1_c2_ro_e_melted <- melt(c1_c2_ro_e)
colnames(c1_c2_ro_e_melted) <- c("Group", "Treatment", "RoE")
c1_c2_ro_e_melted$Group <- factor(c1_c2_ro_e_melted$Group, levels = rev(c("C2", "C1")))
c1_c2_ro_e_melted$Treatment <- factor(c1_c2_ro_e_melted$Treatment, levels = c("RD", "pCR"), labels = c("NR", "R"))

p <- ggplot(c1_c2_ro_e_melted, aes(x = Treatment, y = Group, fill = RoE)) +
  geom_tile() +
  scale_fill_gradient2(low = "#6BAED6", mid = "white", high = "#FB6A6A", midpoint = 1) +
  labs(x = "Treatment", y = "TLS Group", fill = "Ro/e") +
  theme_minimal() +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 9, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    legend.position = "none"
  ) +
  coord_flip()

ggsave("./result/Fig6M_2.pdf", p, width = 5, height = 4)

vmin <- floor(min(c1_c2_ro_e_melted$RoE, na.rm = TRUE) * 10) / 10
vmax <- ceiling(max(c1_c2_ro_e_melted$RoE, na.rm = TRUE) * 10) / 10

p_leg <- ggplot(data.frame(x = 1, y = 1, z = 0), aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#6BAED6",
    mid = "white",
    high = "#FB6A6A",
    midpoint = 1,
    limits = c(vmin, vmax),
    breaks = seq(vmin, vmax, by = 0.1),
    oob = scales::squish,
    name = "Ro/e"
  ) +
  guides(colour = "none") + 
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )

leg <- cowplot::get_legend(p_leg) 

rot_leg <- grid::editGrob(
  leg,
  vp = grid::viewport(angle = 270)
)

grid::grid.newpage()
grid::grid.draw(rot_leg)

pdf("./result/Fig6M_2_legend.pdf", width = 2, height = 3)
grid::grid.draw(rot_leg)
dev.off()

tls_counts_with_meta3 <- tls_counts_with_meta
tls_counts_with_meta3 <- tls_counts_with_meta3[!is.na(tls_counts_with_meta3$pCR.RD), ]
tls_counts_with_meta3$response2 <- tls_counts_with_meta3$pCR.RD
tls_counts_with_meta3$response2 <- ifelse(tls_counts_with_meta3$response2 == "RD", "NR", "R")
tls_counts_with_meta3$response2 <- factor(tls_counts_with_meta3$response2, levels = c("NR", "R"))

p1 <- ggboxplot(tls_counts_with_meta3,
                x = "response2", y = "Immature",
                color = "response2",
                palette = c(
                  "R"  = "#a1bc83",
                  "NR" = "#bc83a1"
                ),
                outlier.shape = NA
) + 
  geom_jitter(aes(color = response2),
              width = 0.2, size = 0.5, alpha = 0.5
  ) +
  stat_compare_means(method = "wilcox.test", label.y = 13, size = 2) +
  labs(title = "E-TLS", y = "Counts") +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 15),
    breaks = seq(0, 15, 5),
    oob = scales::squish
  )

g <- ggplotGrob(p1)
panel_idx <- g$layout[g$layout$name == "panel", ]
g$widths[panel_idx$l] <- unit(100, "pt") 
g$heights[panel_idx$t] <- unit(120, "pt") 

w_in <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)
h_in <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)

pdf("./result/Fig6L_1.pdf",
    width = w_in, height = h_in, useDingbats = FALSE
)
grid.newpage()
grid.draw(g)
dev.off()

p2 <- ggboxplot(tls_counts_with_meta3,
                x = "response2", y = "Primary",
                color = "response2",
                palette = c(
                  "R"  = "#a1bc83",
                  "NR" = "#bc83a1"
                ),
                outlier.shape = NA
) + 
  stat_compare_means(method = "wilcox.test", label.y = 10, size = 2) +
  labs(title = "P-TLS", y = "Counts") +
  geom_jitter(aes(color = response2),
              width = 0.2, size = 0.5, alpha = 0.5
  ) +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 11),
    breaks = seq(0, 11, 5),
    oob = scales::squish
  )

g <- ggplotGrob(p2)
panel_idx <- g$layout[g$layout$name == "panel", ]
g$widths[panel_idx$l] <- unit(100, "pt") 
g$heights[panel_idx$t] <- unit(120, "pt") 

w_in <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)
h_in <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)

pdf("./result/Fig6L_2.pdf",
    width = w_in, height = h_in, useDingbats = FALSE
)
grid.newpage()
grid.draw(g)
dev.off()

p3 <- ggboxplot(tls_counts_with_meta3,
                x = "response2", y = "Secondary",
                color = "response2",
                palette = c(
                  "R"  = "#a1bc83",
                  "NR" = "#bc83a1"
                ),
                outlier.shape = NA
) + 
  stat_compare_means(method = "wilcox.test", label.y = 3, size = 2) +
  labs(title = "S-TLS", y = "Counts") +
  geom_jitter(aes(color = response2),
              width = 0.2, size = 0.5, alpha = 0.5
  ) +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 3),
    breaks = seq(0, 3, 1),
    oob = scales::squish
  )

g <- ggplotGrob(p3)
panel_idx <- g$layout[g$layout$name == "panel", ]
g$widths[panel_idx$l] <- unit(100, "pt") 
g$heights[panel_idx$t] <- unit(120, "pt") 

w_in <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)
h_in <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)

pdf("./result/Fig6L_3.pdf",
    width = w_in, height = h_in, useDingbats = FALSE
)
grid.newpage()
grid.draw(g)
dev.off()

df_long <- tls_counts_with_meta %>%
  group_by(pCR.RD) %>%
  summarise(
    E_TLS = sum(Immature, na.rm = TRUE),
    P_TLS = sum(Primary, na.rm = TRUE),
    S_TLS = sum(Secondary, na.rm = TRUE),
    .groups = "drop"
  )

df_long <- df_long %>%
  pivot_longer(
    cols = c(E_TLS, P_TLS, S_TLS),
    names_to = "maturation",
    values_to = "counts"
  )

tab <- df_long %>%
  tidyr::pivot_wider(
    names_from = maturation,
    values_from = counts,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("pCR.RD") %>%
  as.matrix()

chisq.test(tab)

rm(list = ls())

#### Figure 6N,O ####
tls_counts = readRDS("./data/TLS_classification_BRCA_Farahmand.rds")
breast_trastuzumab_meta <- read_excel("./data_meta/HE_BRCA_Farahmand_meta.xlsx")

tls_counts_with_meta <- merge(tls_counts, breast_trastuzumab_meta, by.x = "sample_id", by.y = "Patient")

tls_long <- tls_counts_with_meta %>%
  pivot_longer(
    cols = c(Immature, Primary, Secondary),
    names_to = "maturation",
    values_to = "count"
  ) %>%
  mutate(
    maturation = factor(maturation,
                        levels = c("Immature", "Primary", "Secondary")
    )
  )

tls_props <- tls_counts_with_meta %>%
  mutate(
    TLS_total = Immature + Primary + Secondary
  ) %>%
  filter(TLS_total > 0)


tls_props_long <- tls_props %>%
  pivot_longer(
    cols = c(Immature, Primary, Secondary),
    names_to = "maturation",
    values_to = "count"
  ) %>%
  mutate(
    maturation = factor(maturation,
                        levels = c("Immature", "Primary", "Secondary")
    ),
    prop = count / TLS_total
  )


tls_props_long <- tls_props %>%
  pivot_longer(
    cols = c(Immature, Primary, Secondary),
    names_to = "maturation",
    values_to = "count"
  ) %>%
  mutate(
    maturation = factor(maturation,
                        levels = c("Immature", "Primary", "Secondary")
    ),
    prop = count / TLS_total
  )



tls_group_props <- tls_props_long %>%
  group_by(Responder, maturation) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop")
tls_group_props$Responder <- factor(tls_group_props$Responder, levels = c("nonresponder", "responder"), labels = c("NR", "R"))

a <- ggplot(tls_group_props, aes(x = Responder, y = mean_prop, fill = maturation)) +
  geom_bar(stat = "identity", position = "fill") + 
  theme_classic(base_size = 6) +
  labs(x = "Response group", y = "Proportion", fill = "Maturation") +
  scale_fill_manual(values = c(
    "Immature" = "#72B28B",
    "Primary" = "#F26522",
    "Secondary" = "#B3B1D8"
  )) +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    legend.text = element_text(size = 8, colour = "black"),
    legend.title = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_flip()

g <- ggplotGrob(a)
panel_idx <- g$layout[g$layout$name == "panel", ]
g$widths[panel_idx$l] <- unit(80, "pt")
g$heights[panel_idx$t] <- unit(72, "pt")

w_in <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)
h_in <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)

pdf("./result/Fig6O_1.pdf",
    width = w_in, height = h_in, useDingbats = FALSE
)
grid.newpage()
grid.draw(g)
dev.off()

patient_table=read.csv("./data_meta/HE_BRCA_Farahmand_meta2.csv")
tls_counts_with_meta2 <- merge(tls_counts_with_meta, patient_table, by = "sample_id")

c1_c2_obs <- as.matrix(table(tls_counts_with_meta2$TLS_Group, tls_counts_with_meta2$Responder))

c1_c2_chi_res <- chisq.test(c1_c2_obs)
c1_c2_expected <- c1_c2_chi_res$expected

c1_c2_ro_e <- c1_c2_obs / c1_c2_expected

rownames(c1_c2_ro_e) <- c("C1", "C2")

c1_c2_ro_e_melted <- melt(c1_c2_ro_e)
colnames(c1_c2_ro_e_melted) <- c("Group", "Treatment", "RoE")
c1_c2_ro_e_melted$Group <- factor(c1_c2_ro_e_melted$Group, levels = rev(c("C2", "C1")))
c1_c2_ro_e_melted$Treatment <- factor(c1_c2_ro_e_melted$Treatment, levels = c("nonresponder", "responder"), labels = c("NR", "R"))

p <- ggplot(c1_c2_ro_e_melted, aes(x = Treatment, y = Group, fill = RoE)) +
  geom_tile() +
  scale_fill_gradient2(low = "#6BAED6", mid = "white", high = "#FB6A6A", midpoint = 1) +
  labs(x = "Treatment", y = "TLS Group", fill = "Ro/e") +
  theme_minimal() +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 9, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    legend.position = "none"
  ) +
  coord_flip()


ggsave("./result/Fig6O_2.pdf", p, width = 5, height = 4)

vmin <- floor(min(c1_c2_ro_e_melted$RoE, na.rm = TRUE) * 10) / 10
vmax <- ceiling(max(c1_c2_ro_e_melted$RoE, na.rm = TRUE) * 10) / 10

p_leg <- ggplot(data.frame(x = 1, y = 1, z = 0), aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#6BAED6",
    mid = "white",
    high = "#FB6A6A",
    midpoint = 1,
    limits = c(vmin, vmax),
    breaks = seq(vmin, vmax, by = 0.1),
    oob = scales::squish,
    name = "Ro/e"
  ) +
  guides(colour = "none") + 
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )

leg <- cowplot::get_legend(p_leg) 

rot_leg <- grid::editGrob(
  leg,
  vp = grid::viewport(angle = 270)
)

grid::grid.newpage()
grid::grid.draw(rot_leg)

pdf("./result/Fig6O_2_legend.pdf", width = 2, height = 3)
grid::grid.draw(rot_leg)
dev.off()


tls_counts_with_meta3 <- tls_counts_with_meta
tls_counts_with_meta3$response2 <- tls_counts_with_meta3$Responder
tls_counts_with_meta3$response2 <- ifelse(tls_counts_with_meta3$response2 == "nonresponder", "NR", "R")
tls_counts_with_meta3$response2 <- factor(tls_counts_with_meta3$response2, levels = c("NR", "R"))

p1 <- ggboxplot(tls_counts_with_meta3,
                x = "response2", y = "Immature",
                color = "response2",
                palette = c(
                  "R"  = "#a1bc83",
                  "NR" = "#bc83a1"
                ),
                outlier.shape = NA
) + 
  geom_jitter(aes(color = response2),
              width = 0.2, size = 0.5, alpha = 0.5,
  ) +
  stat_compare_means(method = "wilcox.test", size = 2) +
  labs(title = "E-TLS", y = "Counts") +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 10),
    breaks = seq(0, 10, 5),
    oob = scales::squish
  )

g <- ggplotGrob(p1)
panel_idx <- g$layout[g$layout$name == "panel", ]
g$widths[panel_idx$l] <- unit(100, "pt") 
g$heights[panel_idx$t] <- unit(120, "pt") 

w_in <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)
h_in <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)

pdf("./result/Fig6N_1.pdf",
    width = w_in, height = h_in, useDingbats = FALSE
)
grid.newpage()
grid.draw(g)
dev.off()




p2 <- ggboxplot(tls_counts_with_meta3,
                x = "response2", y = "Primary",
                color = "response2",
                palette = c(
                  "R"  = "#a1bc83",
                  "NR" = "#bc83a1"
                ),
                outlier.shape = NA
) + 
  geom_jitter(aes(color = response2),
              width = 0.2, size = 0.5, alpha = 0.5
  ) +
  stat_compare_means(method = "wilcox.test", size = 2) +
  labs(title = "P-TLS", y = "Counts") +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 8),
    breaks = seq(0, 8, 4),
    oob = scales::squish
  )

g <- ggplotGrob(p2)
panel_idx <- g$layout[g$layout$name == "panel", ]
g$widths[panel_idx$l] <- unit(100, "pt") 
g$heights[panel_idx$t] <- unit(120, "pt") 

w_in <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)
h_in <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)

pdf("./result/Fig6N_2.pdf",
    width = w_in, height = h_in, useDingbats = FALSE
)
grid.newpage()
grid.draw(g)
dev.off()



p3 <- ggboxplot(tls_counts_with_meta3,
                x = "response2", y = "Secondary",
                color = "response2",
                palette = c(
                  "R"  = "#a1bc83",
                  "NR" = "#bc83a1"
                ),
                outlier.shape = NA
) + 
  geom_jitter(aes(color = response2),
              width = 0.2, size = 0.5, alpha = 0.5
  ) +
  stat_compare_means(method = "wilcox.test", size = 2) +
  labs(title = "S-TLS", y = "Counts") +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, colour = "black"),
    plot.subtitle = element_text(size = 8, colour = "black"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 2),
    breaks = seq(0, 2, 1),
    oob = scales::squish
  )


g <- ggplotGrob(p3)
panel_idx <- g$layout[g$layout$name == "panel", ]
g$widths[panel_idx$l] <- unit(100, "pt") 
g$heights[panel_idx$t] <- unit(120, "pt") 

w_in <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)
h_in <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)

pdf("./result/Fig6N_3.pdf",
    width = w_in, height = h_in, useDingbats = FALSE
)
grid.newpage()
grid.draw(g)
dev.off()

df_long <- tls_counts_with_meta %>%
  group_by(Responder) %>%
  summarise(
    E_TLS = sum(Immature, na.rm = TRUE),
    P_TLS = sum(Primary, na.rm = TRUE),
    S_TLS = sum(Secondary, na.rm = TRUE),
    .groups = "drop"
  )

df_long <- df_long %>%
  pivot_longer(
    cols = c(E_TLS, P_TLS, S_TLS),
    names_to = "maturation",
    values_to = "counts"
  )

tab <- df_long %>%
  tidyr::pivot_wider(
    names_from = maturation,
    values_from = counts,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("Responder") %>%
  as.matrix()

chisq.test(tab)

rm(list = ls())
