### Figure 2 ###
### The codes are separated by Figures ###

### load packages ###
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(ComplexHeatmap)
library(rstatix)
library(tidyr)
library(stats)
library(viridis)
library(png)
library(gtools)
library(stringr)

#### Figure 2A ####
tls_classification = readRDS("./data/TLS_classification.rds")
tls_classification_table = table(tls_classification$Cluster)
tls_classification_table_df = data.frame(tls_classification_table)

tls_classification_table_df <- tls_classification_table_df %>%
  mutate(Var1 = recode(
    Var1,
    "Immature" = "E-TLS",
    "Primary Mature" = "P-TLS",
    "Secondary Mature" = "S-TLS"
  ))

tls_classification_table_df <- tls_classification_table_df %>%
  mutate(
    percentage = round(100 * Freq / sum(Freq), 1),
    label = paste0(Freq, " (", percentage, "%)")
  )

tls_classification_table_df$Var1 <- factor(
  tls_classification_table_df$Var1,
  levels = c("E-TLS","P-TLS","S-TLS")
)

p <- ggplot(tls_classification_table_df,
            aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 2
  ) +
  theme_void() +
  scale_fill_manual(values = c(
    "E-TLS" = "#72b28b",
    "P-TLS" = "#f26522",
    "S-TLS" = "#b3b1d8"
  ))

pdf("./result/Fig2A.pdf", width = 4, height = 4)
print(p)
dev.off()
rm(list = ls())

#### Figure 2C ####
tls_score = readRDS("./data/TLS_score.rds")
tls_classification = readRDS("./data/TLS_classification.rds")

TLS_score_with_classification  = merge(tls_score,tls_classification,by='TLS_ID')

my_comparisons <- list( c("Immature", "Primary Mature"), c("Immature", "Secondary Mature"),
                        c("Primary Mature", "Secondary Mature"))

a5=ggplot(TLS_score_with_classification, aes(x = Cluster, y = tls,fill=Cluster)) +
  geom_boxplot(na.rm = TRUE,alpha=1) + theme_bw() + scale_fill_manual(values = c("#72B28B", "#F26522", "#B3B1D8")) +
  scale_x_discrete(limits=c("Immature","Primary Mature", "Secondary Mature"),labels=c("Immature"="E-TLS", "Primary Mature"="P-TLS", "Secondary Mature"="S-TLS")) +
  stat_compare_means(na.rm=TRUE,comparisons = my_comparisons) + theme_classic()

pdf("./result/Fig2C_1.pdf",width=5,height=5)
a5
dev.off()
rm(list=ls())

tls_CD3D = readRDS("./data/TLS_CD3D.rds")
tls_classification = readRDS("./data/TLS_classification.rds")

tls_CD3D_with_classification  = merge(tls_CD3D,tls_classification,by='TLS_ID')

my_comparisons <- list( c("Immature", "Primary Mature"), c("Immature", "Secondary Mature"),
                        c("Primary Mature", "Secondary Mature"))


a2=ggplot(tls_CD3D_with_classification, aes(x = Cluster, y = CD3D,fill=Cluster)) +
  geom_boxplot(na.rm = TRUE,alpha=1) + theme_bw() + scale_fill_manual(values = c("#72B28B", "#F26522", "#B3B1D8")) +
  scale_x_discrete(limits=c("Immature", "Primary Mature", "Secondary Mature"),labels=c("Immature"="E-TLS", "Primary Mature"="P-TLS", "Secondary Mature"="S-TLS")) +
  stat_compare_means(na.rm=TRUE,comparisons = my_comparisons) + theme_classic()

pdf("./result/Fig2C_2.pdf",width=5,height=5)
a2
dev.off()
rm(list = ls())

tls_MS4A1 = readRDS("./data/TLS_MS4A1.rds")
tls_classification = readRDS("./data/TLS_classification.rds")

tls_MS4A1_with_classification  = merge(tls_MS4A1,tls_classification,by='TLS_ID')

my_comparisons <- list( c("Immature", "Primary Mature"), c("Immature", "Secondary Mature"),
                        c("Primary Mature", "Secondary Mature"))


a3=ggplot(tls_MS4A1_with_classification, aes(x = Cluster, y = MS4A1,fill=Cluster)) +
  geom_boxplot(na.rm = TRUE,alpha=1) + theme_bw() + scale_fill_manual(values = c("#72B28B", "#F26522", "#B3B1D8")) +
  scale_x_discrete(limits=c("Immature", "Primary Mature", "Secondary Mature"),labels=c("Immature"="E-TLS", "Primary Mature"="P-TLS", "Secondary Mature"="S-TLS")) +
  stat_compare_means(na.rm=TRUE,comparisons = my_comparisons) + theme_classic()

pdf("./result/Fig2C_3.pdf",width=5,height=5)
a3
dev.off()
rm(list = ls())

tls_CR2 = readRDS("./data/TLS_CR2.rds")
tls_classification = readRDS("./data/TLS_classification.rds")

tls_CR2_with_classification  = merge(tls_CR2,tls_classification,by='TLS_ID')

my_comparisons <- list( c("Immature", "Primary Mature"), c("Immature", "Secondary Mature"),
                        c("Primary Mature", "Secondary Mature"))


a4=ggplot(tls_CR2_with_classification, aes(x = Cluster, y = CR2,fill=Cluster)) +
  geom_boxplot(na.rm = TRUE,alpha=1) + theme_bw() + scale_fill_manual(values = c("#72B28B", "#F26522", "#B3B1D8")) +
  scale_x_discrete(limits=c("Immature", "Primary Mature", "Secondary Mature"),labels=c("Immature"="E-TLS", "Primary Mature"="P-TLS", "Secondary Mature"="S-TLS")) +
  stat_compare_means(na.rm=TRUE,comparisons = my_comparisons) + theme_classic()

pdf("./result/Fig2C_4.pdf",width=5,height=5)
a4
dev.off()
rm(list = ls())

tls_FCER2 = readRDS("./data/TLS_FCER2.rds")
tls_classification = readRDS("./data/TLS_classification.rds")

tls_FCER2_with_classification  = merge(tls_FCER2,tls_classification,by='TLS_ID')

my_comparisons <- list( c("Immature", "Primary Mature"), c("Immature", "Secondary Mature"),
                        c("Primary Mature", "Secondary Mature"))


a5=ggplot(tls_FCER2_with_classification, aes(x = Cluster, y = FCER2,fill=Cluster)) +
  geom_boxplot(na.rm = TRUE,alpha=1) + theme_bw() + scale_fill_manual(values = c("#72B28B", "#F26522", "#B3B1D8")) +
  scale_x_discrete(limits=c("Immature", "Primary Mature", "Secondary Mature"),labels=c("Immature"="E-TLS", "Primary Mature"="P-TLS", "Secondary Mature"="S-TLS")) +
  stat_compare_means(na.rm=TRUE,comparisons = my_comparisons) + theme_classic()

pdf("./result/Fig2C_5.pdf",width=5,height=5)
a5
dev.off()
rm(list = ls())

tls_MKI67 = readRDS("./data/TLS_MKI67.rds")
tls_classification = readRDS("./data/TLS_classification.rds")

tls_MKI67_with_classification  = merge(tls_MKI67,tls_classification,by='TLS_ID')

my_comparisons <- list( c("Immature", "Primary Mature"), c("Immature", "Secondary Mature"),
                        c("Primary Mature", "Secondary Mature"))


a6=ggplot(tls_MKI67_with_classification, aes(x = Cluster, y = MKI67,fill=Cluster)) +
  geom_boxplot(na.rm = TRUE,alpha=1) + theme_bw() + scale_fill_manual(values = c("#72B28B", "#F26522", "#B3B1D8")) +
  scale_x_discrete(limits=c("Immature", "Primary Mature", "Secondary Mature"),labels=c("Immature"="E-TLS", "Primary Mature"="P-TLS", "Secondary Mature"="S-TLS")) +
  stat_compare_means(na.rm=TRUE,comparisons = my_comparisons) + theme_classic()

pdf("./result/Fig2C_6.pdf",width=5,height=5)
a6
dev.off()
rm(list = ls())

tls_BCL6 = readRDS("./data/TLS_BCL6.rds")
tls_classification = readRDS("./data/TLS_classification.rds")

tls_BCL6_with_classification  = merge(tls_BCL6,tls_classification,by='TLS_ID')

my_comparisons <- list( c("Immature", "Primary Mature"), c("Immature", "Secondary Mature"),
                        c("Primary Mature", "Secondary Mature"))


a7=ggplot(tls_BCL6_with_classification, aes(x = Cluster, y = BCL6,fill=Cluster)) +
  geom_boxplot(na.rm = TRUE,alpha=1) + theme_bw() + scale_fill_manual(values = c("#72B28B", "#F26522", "#B3B1D8")) +
  scale_x_discrete(limits=c("Immature", "Primary Mature", "Secondary Mature"),labels=c("Immature"="E-TLS", "Primary Mature"="P-TLS", "Secondary Mature"="S-TLS")) +
  stat_compare_means(na.rm=TRUE,comparisons = my_comparisons) + theme_classic()

pdf("./result/Fig2C_7.pdf",width=5,height=5)
a7
dev.off()
rm(list = ls())


#### Figure 2E ####
tls_classification = readRDS("./data/TLS_classification.rds")
tls_location <- readRDS("./data/TLS_location.rds")
tls_location = tls_location[,c("TLS_ID","Location")]
tls_classification2 = tls_classification
tls_classification2 = tls_classification2[tls_classification2$TLS_ID %in% tls_location$TLS_ID,]

tls_classification2 = merge(tls_classification2,tls_location,by='TLS_ID')

maturation_counts2 <- tls_classification2 %>%
  dplyr::group_by(Cluster,Location) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  ungroup() %>%
  group_by(Location) %>%
  mutate(total_count = sum(count),
         proportion = count / total_count) %>%
  ungroup() %>%
  arrange(Cluster, count)


maturation_counts2 = data.frame(maturation_counts2)
maturation_counts2$Location = factor(maturation_counts2$Location,levels=rev(c('IT',"PT","DT")))

plot_df = maturation_counts2

maturation_counts2_long = maturation_counts2[,c('Cluster','Location','count')]

maturation_counts2_long = maturation_counts2_long %>% pivot_wider(names_from = Cluster, values_from = count, values_fill = list(count = 0))
maturation_counts2_long = data.frame(maturation_counts2_long)
rownames(maturation_counts2_long) = maturation_counts2_long$Location
maturation_counts2_long = maturation_counts2_long[,-c(1)]

chisq_test_result <- chisq.test(maturation_counts2_long)
print(chisq_test_result)

residuals <- chisq_test_result$residuals

adjusted_residuals <- residuals / sqrt(1 - chisq_test_result$observed / sum(chisq_test_result$observed))

p_values_adjusted_residuals <- 2 * pnorm(-abs(adjusted_residuals))

categorize_enrichment_pvalue <- function(p, residual) {
  if (residual > 0) {
    if (p < 0.001) {
      return("***")
    } else if (p < 0.01) {
      return("**")
    } else if (p < 0.05) {
      return("*")
    } else {
      return("")
    }
  } else {
    return("")
  }
}
p_values_enrichment_categories <- mapply(categorize_enrichment_pvalue, p_values_adjusted_residuals, adjusted_residuals)

p_values_enrichment_categories_matrix <- matrix(p_values_enrichment_categories,
                                                nrow = nrow(p_values_adjusted_residuals),
                                                dimnames = dimnames(p_values_adjusted_residuals))

print(p_values_enrichment_categories_matrix)

ann_df <- as.data.frame(as.table(p_values_enrichment_categories_matrix))
colnames(ann_df) <- c('Location','Cluster', "sig")
ann_df <- ann_df %>% filter(sig != "")


plot_df$Cluster     <- factor(as.character(plot_df$Cluster), levels = c("Immature","Primary Mature","Secondary Mature"))
ann_df$Cluster      <- factor(as.character(ann_df$Cluster), levels = c("Immature","Primary.Mature","Secondary.Mature"))

recode_map <- c(
  "Immature" = "E-TLS",
  "Primary Mature" = "P-TLS",
  "Primary.Mature" = "P-TLS",
  "Secondary Mature" = "S-TLS",
  "Secondary.Mature" = "S-TLS"
)

plot_df$Cluster <- factor(
  recode_map[as.character(plot_df$Cluster)],
  levels = c("E-TLS","P-TLS","S-TLS")
)

ann_df$Cluster <- factor(
  recode_map[as.character(ann_df$Cluster)],
  levels = c("E-TLS","P-TLS","S-TLS")
)
plot_df2 <- plot_df %>%
  complete( Location, Cluster,fill = list(count = 0, proportion = 0)) %>%
  group_by(Location) %>%
  mutate(total_count = sum(count),
         proportion  = ifelse(total_count > 0, count / total_count, 0)) %>%
  ungroup() %>%
  left_join(ann_df, by = c("Cluster","Location"))

y_var= "proportion"
p <- ggplot(plot_df2, aes(x = Location, y = proportion, fill = Cluster)) +
  geom_col() +
  geom_text(
    aes(label = sig),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 4,
    na.rm = TRUE
  ) +
  theme_classic() +
  scale_fill_manual(values = c("#72B28B", "#F26522", "#B3B1D8")) +
  scale_y_continuous(limits = if (y_var == "proportion") c(0, 1) else waiver(),
                     expand = c(0, 0))

pdf("./result/Fig2E.pdf", width = 6.5, height = 4.5)
print(p)
dev.off()

rm(list = ls())

#### Figure 2F ####
tls_density <- readRDS("./data/TLS_density.rds")
tls_density = tls_density[,c("TLS_ID","density")]
tls_classification = readRDS("./data/TLS_classification.rds")
tls_classification2 = tls_classification
tls_classification2 = tls_classification2[tls_classification2$TLS_ID %in% tls_density$TLS_ID,]

tls_classification2 = merge(tls_classification2,tls_density,by='TLS_ID')

tls_classification2$Cluster = ifelse(tls_classification2$Cluster == "Immature","E-TLS",
                                     ifelse(tls_classification2$Cluster == "Primary Mature","P-TLS","S-TLS"))

my_comparisons <- list( c("E-TLS", "P-TLS"), c("E-TLS", "S-TLS"),
                        c("P-TLS", "S-TLS"))

p_tls <- ggplot(tls_classification2, aes(x=Cluster, y=density, fill=Cluster)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  stat_compare_means(comparisons = my_comparisons, label = "p.format",label.y = c(360000, 410000, 460000),tip.length = 0.00) +
  scale_y_continuous(breaks=seq(0, 500000, by=100000)) +
  coord_cartesian(ylim = c(0, 600000)) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  scale_fill_manual(values = c("S-TLS"="#B3B1D8","P-TLS"="#F26522","E-TLS"="#72B28B"))


pdf("./result/Fig2F.pdf",width=4,height=5)
p_tls
dev.off()
rm(list = ls())

#### Figure 2G ####
tls_classification = readRDS("./data/TLS_classification.rds")

cancer_type_color <- c(
  "LUAD" = "#A6CEE3",
  "STAD" = "#1F78B4",
  "KIRC" ="#B2DF8A",
  "BRCA"= "#33A02C",
  "CSCC"= "#FB9A99",
  "OSCC" = "#E31A1C",
  "BLCA" ="#FDBF6F",
  "PAAD" = "#FF7F00",
  "CRC" = "#CAB2D6",
  "OVCA"= "#6A3D9A",
  "LIHC" = "#B15928"
)


x_immature = tls_classification[tls_classification$Cluster %in% c("Immature"),]
cluster_immature_counts <- x_immature %>%
  dplyr::group_by(Cluster, Cancer_type) %>%
  dplyr::summarise(count =  dplyr::n(),.groups = "drop") %>%
  ungroup()

cluster_immature_counts = data.frame(cluster_immature_counts)
cluster_immature_counts = rbind(cluster_immature_counts,c("Immature","OSCC",0))
cluster_immature_counts$count = as.numeric(cluster_immature_counts$count)

cluster_immature_counts = cluster_immature_counts[order(cluster_immature_counts$count),]
cluster_immature_counts$Cancer_type <- factor(cluster_immature_counts$Cancer_type,levels=cluster_immature_counts$Cancer_type)
cluster_immature_counts$Cancer_type2 <- factor(cluster_immature_counts$Cancer_type,levels=rev(cluster_immature_counts$Cancer_type))

b = ggplot(cluster_immature_counts,aes(x=Cancer_type,y=count))+geom_bar(stat='identity',fill="#71b28b") +
  theme_classic()  +
  coord_flip() + labs(fill = "") + ylab("") + xlab("")

pdf("./result/Fig2G_1-1.pdf")
print(b)
dev.off()


cluster_immature_counts$prop = cluster_immature_counts$count / sum(cluster_immature_counts$count)
cluster_immature_counts$prop = cluster_immature_counts$prop *100

cluster_immature_counts <- cluster_immature_counts %>%
  mutate(lab.ypos = cumsum(round(prop)) - 0.5*round(prop))

p_tls <- ggplot(cluster_immature_counts, aes(x='', y=prop,fill = Cancer_type2)) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0)  + theme_void() +
  scale_fill_manual(values = cancer_type_color)

pdf("./result/Fig2G_1-2.pdf")
p_tls
dev.off()

x_primary = tls_classification[tls_classification$Cluster %in% c("Primary Mature"),]

cluster_primary_counts <- x_primary %>%
  dplyr::group_by(Cluster, Cancer_type) %>%
  dplyr::summarise(count =  dplyr::n(),.groups = "drop") %>%
  ungroup()

cluster_primary_counts = data.frame(cluster_primary_counts)

cluster_primary_counts = rbind(cluster_primary_counts,c("Primary Mature","OSCC",0))
cluster_primary_counts = rbind(cluster_primary_counts,c("Primary Mature","OVCA",0))

cluster_primary_counts$count = as.numeric(cluster_primary_counts$count)

cluster_primary_counts$Cancer_type <- factor(cluster_primary_counts$Cancer_type,levels=cluster_immature_counts$Cancer_type)
cluster_primary_counts$Cancer_type2 <- factor(cluster_primary_counts$Cancer_type,levels=rev(cluster_immature_counts$Cancer_type))

b2 = ggplot(cluster_primary_counts,aes(x=Cancer_type,y=count))+geom_bar(stat='identity',fill="#f16623") +
  theme_classic()  +
  coord_flip() + labs(fill = "") + ylab("") + xlab("")

pdf("./result/Fig2G_2-1.pdf")
print(b2)
dev.off()

cluster_primary_counts$prop = cluster_primary_counts$count / sum(cluster_primary_counts$count)
cluster_primary_counts$prop = cluster_primary_counts$prop *100

cluster_primary_counts <- cluster_primary_counts %>%
  mutate(lab.ypos = cumsum(round(prop)) - 0.5*round(prop))

p_tls2 <- ggplot(cluster_primary_counts, aes(x='', y=prop,fill = Cancer_type2)) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0)  + theme_void() +
  scale_fill_manual(values = cancer_type_color)

pdf("./result/Fig2G_2-2.pdf")
p_tls2
dev.off()


x_secondary = tls_classification[tls_classification$Cluster %in% c("Secondary Mature"),]


cluster_secondary_counts <- x_secondary %>%
  dplyr::group_by(Cluster, Cancer_type) %>%
  dplyr::summarise(count =  dplyr::n(),.groups = "drop") %>%
  ungroup()

cluster_secondary_counts = data.frame(cluster_secondary_counts)

cluster_secondary_counts = rbind(cluster_secondary_counts,c("Secondary Mature","BLCA",0))
cluster_secondary_counts = rbind(cluster_secondary_counts,c("Secondary Mature","PAAD",0))

cluster_secondary_counts$count = as.numeric(cluster_secondary_counts$count)

cluster_secondary_counts$Cancer_type <- factor(cluster_secondary_counts$Cancer_type,levels=cluster_immature_counts$Cancer_type)
cluster_secondary_counts$Cancer_type2 <- factor(cluster_secondary_counts$Cancer_type,levels=rev(cluster_immature_counts$Cancer_type))

b3 = ggplot(cluster_secondary_counts,aes(x=Cancer_type,y=count))+geom_bar(stat='identity',fill="#b3b1d8") +
  theme_classic()  +
  coord_flip() + labs(fill = "") + ylab("") + xlab("")

pdf("./result/Fig2G_3-1.pdf")
print(b3)
dev.off()


cluster_secondary_counts$prop = cluster_secondary_counts$count / sum(cluster_secondary_counts$count)
cluster_secondary_counts$prop = cluster_secondary_counts$prop *100

cluster_secondary_counts <- cluster_secondary_counts %>%
  mutate(lab.ypos = cumsum(round(prop)) - 0.5*round(prop))

p_tls3 <- ggplot(cluster_secondary_counts, aes(x='', y=prop,fill = Cancer_type2)) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0)  + theme_void() +
  scale_fill_manual(values = cancer_type_color)

pdf("./result/Fig2G_3-2.pdf")
p_tls3
dev.off()

rm(list = ls())

#### Figure 2H ####
tls_classification = readRDS("./data/TLS_classification.rds")

cnc_type <- c("LUAD","STAD","LIHC","BRCA","CRC","CSCC","BLCA","KIRC","OSCC","PAAD","OVCA","GBM")

maturation_counts <- tls_classification %>%
  dplyr::group_by(Cluster, Cancer_type) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  group_by(Cancer_type) %>%
  mutate(total_count = sum(count),
         proportion  = count / total_count) %>%
  ungroup()

maturation_counts$Cancer_type <- factor(maturation_counts$Cancer_type, levels = cnc_type[1:11])

plot_df <- maturation_counts %>% filter(Cancer_type != "OSCC")

cluster_levels <- sort(unique(as.character(plot_df$Cluster)))
plot_df$Cluster <- factor(plot_df$Cluster, levels = cluster_levels)

count_wide <- plot_df %>%
  select(Cancer_type, Cluster, count) %>%
  tidyr::pivot_wider(
    names_from  = Cluster,
    values_from = count,
    values_fill = list(count = 0)
  ) %>%
  as.data.frame()

rownames(count_wide) <- count_wide$Cancer_type
count_mat <- as.matrix(count_wide[, setdiff(colnames(count_wide), "Cancer_type"), drop = FALSE])

chisq_test_result <- chisq.test(count_mat)
message("Chi-square test:")
print(chisq_test_result)

residuals <- chisq_test_result$residuals
adjusted_residuals <- residuals / sqrt(1 - chisq_test_result$observed / sum(chisq_test_result$observed))
p_values_adjusted_residuals <- 2 * pnorm(-abs(adjusted_residuals))

categorize_enrichment_pvalue <- function(p, residual) {
  if (residual > 0) {
    if (p < 0.001) "***"
    else if (p < 0.01) "**"
    else if (p < 0.05) "*"
    else ""
  } else {
    ""
  }
}

p_values_enrichment_categories <- mapply(
  categorize_enrichment_pvalue,
  p_values_adjusted_residuals,
  adjusted_residuals
)

p_values_enrichment_categories_matrix <- matrix(
  p_values_enrichment_categories,
  nrow = nrow(p_values_adjusted_residuals),
  dimnames = dimnames(p_values_adjusted_residuals)
)

message("Enrichment categories matrix (rows=cancer, cols=cluster):")
print(p_values_enrichment_categories_matrix)

ann_df <- as.data.frame(as.table(p_values_enrichment_categories_matrix))
colnames(ann_df) <- c("Cancer_type", "Cluster", "sig")

ann_df <- ann_df %>% filter(sig != "Neg")

ann_df$Cancer_type <- factor(as.character(ann_df$Cancer_type), levels = levels(plot_df$Cancer_type))
ann_df$Cluster     <- factor(as.character(ann_df$Cluster), levels = levels(plot_df$Cluster))

plot_df$Cluster <- factor(as.character(plot_df$Cluster), levels = (c("Immature","Primary Mature","Secondary Mature")))

plot_df2 <- plot_df %>%
  arrange(Cancer_type, Cluster) %>%
  group_by(Cancer_type) %>%
  mutate(y_mid = cumsum(proportion) - proportion / 2) %>%
  ungroup() %>%
  left_join(ann_df, by = c("Cancer_type", "Cluster"))


p <- ggplot(plot_df2, aes(fill = Cluster, y = proportion, x = Cancer_type)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_text(
    aes(label = sig),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 4,
    na.rm = TRUE
  ) +
  theme_classic() +
  scale_fill_manual(values = (c("#72B28B", "#F26522", "#B3B1D8"))) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))

pdf("./result/Fig2H.pdf", width = 6.5, height = 4.5)
print(p)
dev.off()

rm(list = ls())

#### Figure 2I ####
tls_classification = readRDS("./data/TLS_classification.rds")
tls_location <- readRDS("./data/TLS_location.rds")
tls_location = tls_location[,c("TLS_ID","Location")]
tls_classification2 = tls_classification
tls_classification2 = tls_classification2[tls_classification2$TLS_ID %in% tls_location$TLS_ID,]

tls_classification2 = merge(tls_classification2,tls_location,by='TLS_ID')

cnc_type <- c("LUAD","STAD","LIHC","BRCA","CRC","CSCC","BLCA","KIRC","OSCC","PAAD","OVCA","GBM")

tls_classification2$mat_loc = paste0(tls_classification2$Cluster,"_",tls_classification2$Location)

maturation_counts3 <- tls_classification2 %>%
  dplyr::group_by(mat_loc, Cancer_type) %>%
  dplyr::summarise(count = n(),.groups='drop') %>%
  ungroup() %>%
  group_by(Cancer_type) %>%
  mutate(total_count = sum(count),
         proportion = count / total_count) %>%
  ungroup() %>%
  arrange(Cancer_type, count)

maturation_counts3 = data.frame(maturation_counts3)

maturation_counts3$Cancer_type <- factor(maturation_counts3$Cancer_type, levels = rev(cnc_type[1:11]))

plot_df <- maturation_counts3 %>% filter(Cancer_type != "OSCC")

plot_df <- plot_df %>%
  group_by(Cancer_type) %>%
  mutate(total_count = sum(count),
         proportion  = ifelse(total_count > 0, count / total_count, 0)) %>%
  ungroup()

count_wide <- plot_df %>%
  select(Cancer_type, mat_loc, count) %>%
  pivot_wider(names_from = mat_loc, values_from = count, values_fill = list(count = 0)) %>%
  as.data.frame()

rownames(count_wide) <- count_wide$Cancer_type
count_mat <- as.matrix(count_wide[, setdiff(colnames(count_wide), "Cancer_type"), drop = FALSE])

chisq_test_result <- chisq.test(count_mat)
print(chisq_test_result)

residuals <- chisq_test_result$residuals
adjusted_residuals <- residuals / sqrt(1 - chisq_test_result$observed / sum(chisq_test_result$observed))
pvals <- 2 * pnorm(-abs(adjusted_residuals))

categorize_enrichment_pvalue <- function(p, residual) {
  if (residual > 0) {
    if (p < 0.001) "***"
    else if (p < 0.01) "**"
    else if (p < 0.05) "*"
    else ""
  } else {
    ""
  }
}

p_cat <- mapply(categorize_enrichment_pvalue, pvals, adjusted_residuals)
p_values_enrichment_categories_matrix <- matrix(
  p_cat,
  nrow = nrow(pvals),
  dimnames = dimnames(pvals)
)
print(p_values_enrichment_categories_matrix)

ann_df <- as.data.frame(as.table(p_values_enrichment_categories_matrix))
colnames(ann_df) <- c("Cancer_type", "mat_loc", "sig")
ann_df <- ann_df %>% filter(sig != "")

mcnv_levels <- colnames(count_mat)

plot_df$Cancer_type     <- factor(as.character(plot_df$Cancer_type), levels = levels(maturation_counts3$Cancer_type))
plot_df$mat_loc  <- factor(as.character(plot_df$mat_loc), levels = mcnv_levels)

ann_df$Cancer_type      <- factor(as.character(ann_df$Cancer_type), levels = levels(plot_df$Cancer_type))
ann_df$mat_loc   <- factor(as.character(ann_df$mat_loc), levels = mcnv_levels)

plot_df2 <- plot_df %>%
  complete(Cancer_type, mat_loc, fill = list(count = 0, proportion = 0)) %>%
  group_by(Cancer_type) %>%
  mutate(total_count = sum(count),
         proportion  = ifelse(total_count > 0, count / total_count, 0)) %>%
  ungroup() %>%
  left_join(ann_df, by = c("Cancer_type","mat_loc"))
plot_df2 = plot_df2[plot_df2$Cancer_type != 'OSCC',]
plot_df2$Cancer_type <- droplevels(plot_df2$Cancer_type)
plot_df2$Cancer_type <- factor(plot_df2$Cancer_type,levels = rev(levels(plot_df2$Cancer_type)))
plot_df2$mat_loc = factor(plot_df2$mat_loc,
                          levels = c(
                            "Immature_DT",
                            "Immature_PT",
                            "Immature_IT",
                            "Primary Mature_DT",
                            "Primary Mature_PT",
                            "Primary Mature_IT",
                            "Secondary Mature_DT",
                            "Secondary Mature_PT",
                            "Secondary Mature_IT"
                          ),
                          labels = c(
                            "E-TLS_DT",
                            "E-TLS_PT",
                            "E-TLS_IT",
                            "P-TLS_DT",
                            "P-TLS_PT",
                            "P-TLS_IT",
                            "S-TLS_DT",
                            "S-TLS_PT",
                            "S-TLS_IT"
                          ))
y_var='proportion'

p <- ggplot(plot_df2, aes(x = Cancer_type, y = proportion, fill = mat_loc)) +
  geom_col() +
  geom_text(
    aes(label = sig),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 4,
    na.rm = TRUE
  ) +
  theme_classic() +
  scale_fill_manual(values = c("#365f3c","#03a087","#c1e2c2","#bf4327","#ea5a24","#f39b80","#7c287d","#6a52a2","#bda0cc")) +
  scale_y_continuous(limits = if (y_var == "proportion") c(0, 1) else waiver(),
                     expand = c(0, 0))

pdf("./result/Fig2I.pdf", width = 6.5, height = 4.5)
print(p)
dev.off()

rm(list = ls())
