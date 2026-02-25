### Figure 1 ###
### The codes are separated by Figures ###

### load packages ###
library(gtools)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(png)
library(dplyr)
library(tidyr)

### Figure 1A and 1B were created with Biorender ###
### ST TLS detection pipeline is described in Fig1.sh ###

#### Figure 1C ####
tls_counts = readRDS("./data/TLS_counts.rds")
tls_counts_table = tls_counts[!duplicated(tls_counts$Sample_ID),'Cancer_type']
tls_counts_table = table(tls_counts_table)

cancer_types <- sapply(strsplit(names(tls_counts_table), "_"), function(x) x[1])

cancer_counts <- setNames(as.numeric(tls_counts_table), cancer_types)

aggregated_counts <- tapply(cancer_counts, INDEX = names(cancer_counts), FUN = mean)

summary_data <- data.frame(cancer_type = names(aggregated_counts), total_counts = aggregated_counts)
summary_data[12,] = c("GBM",0)
summary_data$cancer_type = c("BLCA","BRCA","KIRC","CRC","CSCC","STAD","LIHC","LUAD","OSCC","OVCA","PAAD","GBM")
cnc_type= c("LUAD","STAD","KIRC","BRCA","CSCC","OSCC","BLCA","PAAD","CRC","OVCA","LIHC","GBM")
summary_data = summary_data[match(cnc_type,summary_data$cancer_type),]

cnt= c(10,40,26,61,16,12,14,26,54,20,25,36)
summary_data = cbind(summary_data,cnt)

summary_data$prop = as.numeric(summary_data$total_counts) / as.numeric(summary_data$cnt)
TLS_samples = summary_data[,c('cancer_type','prop')]
TLS_samples$type = 'With TLS'
TLS_samples = TLS_samples[order(TLS_samples$prop,decreasing=T),]
prop = 1 - TLS_samples$prop
type = 'Without TLS'
TLS_samples2 = data.frame(cbind(TLS_samples$cancer_type,prop,type))
colnames(TLS_samples2) = c('cancer_type','prop','type')
TLS_samples = rbind(TLS_samples,TLS_samples2)

TLS_samples$type = factor(TLS_samples$type, levels = c("Without TLS","With TLS"))
TLS_samples$cancer_type <- factor(TLS_samples$cancer_type, levels = rev(TLS_samples$cancer_type[1:12]))
TLS_samples$prop = as.numeric(TLS_samples$prop)

a = ggplot(TLS_samples,aes(x=cancer_type,y=prop,fill=type))+geom_bar(stat='identity')+
  scale_fill_manual(values=c( "yellow3","green4"))+theme_classic() +
  scale_y_continuous(limits = c(0, 1.01), expand = c(0,0)) +
  coord_flip()

pdf("./result/Fig1C_1.pdf",width=5,height=5)
print(a)
dev.off()

summary_data$total_counts = as.numeric(summary_data$total_counts)
summary_data2 = summary_data
summary_data2$cancer_type = factor(summary_data2$cancer_type,levels = rev(as.character(TLS_samples$cancer_type)[1:12]))
b = ggplot(summary_data2,aes(x=cancer_type,y=total_counts))+geom_bar(stat='identity',fill="green4")+
  theme_classic() + scale_y_continuous(breaks=c(1,5,10,15,20,25,30,35,40)) +
  coord_flip() + labs(fill = "") + ylab("") + xlab("")

pdf("./result/Fig1C_2.pdf")
print(b)
dev.off()
rm(list = ls())

#### Figure 1D ####
tls_density <- readRDS("./data/TLS_density.rds")
cnc_type= c("LUAD","STAD","LIHC","BRCA",
            "CRC","CSCC","BLCA","KIRC",
            "OSCC","PAAD","OVCA","GBM")
tls_density$Cancer_type <- factor(tls_density$Cancer_type, levels = rev(cnc_type[1:11]))

tls_density = tls_density[!tls_density$Cancer_type %in% c("OSCC","PAAD","OVCA"),]

p_tls <- ggplot(tls_density, aes(x=Cancer_type, y=density, fill=Cancer_type)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  scale_y_continuous(limits=c(0, 500000), breaks=seq(0, 500000, by=100000)) +
  scale_fill_brewer(palette="Paired") +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  coord_flip()

pdf("./result/Fig1D_1.pdf",width=4,height=5)
p_tls
dev.off()
rm(list = ls())

tls_density <- readRDS("./data/TLS_density.rds")
tls_density$Sample_ID = sapply(strsplit(tls_density$TLS_ID, "_"), function(x) paste(x[1]))
tls_sample_count <- tls_density %>%
  count(Cancer_type, Sample_ID, name = "Freq")

cnc_type= c("LUAD","STAD","LIHC","BRCA",
            "CRC","CSCC","BLCA","KIRC",
            "OSCC","PAAD","OVCA","GBM")

tls_sample_count$Cancer_type <- factor(tls_sample_count$Cancer_type, levels = rev(cnc_type[1:11]))
tls_sample_count = tls_sample_count[!tls_sample_count$Cancer_type %in% c("OSCC","PAAD","OVCA"),]

a <- ggplot(tls_sample_count, aes(x=Cancer_type, y=Freq,fill = Cancer_type)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette="Paired")+
  theme(legend.position = "none",axis.title.x = element_blank())  + coord_flip()

pdf("./result/Fig1D_2.pdf",width=3,height=5)
a
dev.off()

rm(list = ls())

#### Figure 1E ####
tls_location <- readRDS("./data/TLS_location.rds")

tls_location_table = table(tls_location$Location)
tls_location_df = data.frame(tls_location_table)
tls_location_df$percentage = 100* (tls_location_df$Freq / sum(tls_location_df$Freq))
tls_location_df$percentage = round(tls_location_df$percentage,1)

a <- ggplot(tls_location_df, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(
    aes(label = paste0(percentage, "%")),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 5
  ) +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(
    values = c(
      DT = "#3d5589",
      PT = "#8491b4",
      IT = "#c0d5ee"
    )
  )

pdf("./result/Fig1E_1.pdf")
print(a)
dev.off()

rm(list = ls())

tls_location <- readRDS("./data/TLS_location.rds")

tls_location <- tls_location[tls_location$Cancer_type != "OSCC", ]

cluster_type_counts <- tls_location %>%
  dplyr::group_by(Location, Cancer_type) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Cancer_type) %>%
  dplyr::mutate(
    total_count = sum(count),
    proportion  = count / total_count
  ) %>%
  dplyr::ungroup()

cluster_type_counts$Cancer_type <- factor(
  cluster_type_counts$Cancer_type,
  levels = rev(cnc_type[1:11])
)

cluster_type_counts$Location <- factor(cluster_type_counts$Location, levels = c("DT","PT","IT"))

cluster_wide <- cluster_type_counts %>%
  dplyr::select(Location, Cancer_type, count) %>%
  tidyr::pivot_wider(names_from = Location, values_from = count, values_fill = 0) %>%
  as.data.frame()

rownames(cluster_wide) <- cluster_wide$Cancer_type
cluster_wide <- cluster_wide[, setdiff(colnames(cluster_wide), "Cancer_type"), drop = FALSE]

chisq_test_result <- chisq.test(cluster_wide)
print(chisq_test_result)

residuals <- chisq_test_result$residuals
adjusted_residuals <- residuals / sqrt(1 - chisq_test_result$observed / sum(chisq_test_result$observed))
p_values_adjusted_residuals <- 2 * pnorm(-abs(adjusted_residuals))

categorize_enrichment_pvalue <- function(p, residual) {
  if (is.na(p) || is.na(residual)) return("")
  if (residual > 0) {
    if (p < 0.001) return("***")
    if (p < 0.01)  return("**")
    if (p < 0.05)  return("*")
    return("")
  } else {
    return("")
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

desired_order <- c(
  "LUAD",
  "STAD",
  "LIHC",
  "BRCA",
  "CRC",
  "CSCC",
  "BLCA",
  "KIRC",
  "PAAD"
)

keep_rows <- intersect(desired_order, rownames(adjusted_residuals))

adjusted_residuals <- adjusted_residuals[keep_rows, , drop = FALSE]
p_values_adjusted_residuals <- p_values_adjusted_residuals[keep_rows, , drop = FALSE]
p_values_enrichment_categories_matrix <- p_values_enrichment_categories_matrix[keep_rows, , drop = FALSE]

pval_long <- as.data.frame(as.table(p_values_enrichment_categories_matrix)) %>%
  dplyr::rename(Cancer_type = Var1, Location = Var2, sig = Freq) %>%
  mutate(
    Cancer_type = as.character(Cancer_type),
    Location     = as.character(Location),
    sig         = as.character(sig)
  )

plot_df <- cluster_type_counts %>%
  left_join(pval_long, by = c("Cancer_type", "Location")) %>%
  mutate(sig = ifelse(is.na(sig), "", sig))

plot_df <- plot_df %>%
  dplyr::group_by(Cancer_type) %>%
  dplyr::arrange(Location) %>%
  dplyr::mutate(ypos = cumsum(proportion) - 0.5 * proportion) %>%
  dplyr::ungroup()

fill_cols <- c(
  DT = "#3d5589",
  PT = "#8491b4",
  IT = "#c0d5ee"
)

plot_df <- plot_df %>%
  dplyr::filter(Cancer_type %in% desired_order) %>%
  dplyr::mutate(Cancer_type = factor(Cancer_type, levels = rev(desired_order)),
                Location = factor(Location, levels = c("DT","IT","PT")))


p <- ggplot(plot_df, aes(x = Cancer_type, y = proportion, fill = Location)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = ifelse(sig == "", "", sig)),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 4
  ) +
  coord_flip() +
  theme_classic() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(values = c(DT="#3d5589", IT="#c0d5ee", PT="#8491b4"))

pdf("./result/Fig1E_2.pdf", width = 5, height = 5)
print(p)
dev.off()

