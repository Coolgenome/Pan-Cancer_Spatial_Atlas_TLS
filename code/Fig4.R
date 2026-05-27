### Figure 4 ###
### The codes are separated by Figures ###

### load packages ###
library(monocle)
library(pheatmap)
library(reshape2)
library(dplyr)
library(grid)
library(ggplot2)
library(RColorBrewer)

### Figure 4A was created with Biorender ###

#### Figure 4B ####
cor_matrix = list()

for(j in c(1,3:8,10)){
  tls_object = "IT"
  cancer_index = j
  load(paste0("./data/","ST_",cancer_index,"_","combined_sig_",tls_object,"_objects.RData"))
  sig_variable_name <- paste0("combined_sig_", tls_object)
  meta_data_variable_name <- paste0("combined_meta_data_", tls_object)
  
  expr_mat <- get(sig_variable_name)
  cell_metadata_df <- get(meta_data_variable_name)
  
  gene_metadata <- data.frame(gene_short_name = rownames(expr_mat), row.names = rownames(expr_mat))
  cell_metadata <- cell_metadata_df
  pd <- new("AnnotatedDataFrame", data = cell_metadata)
  fd <- new("AnnotatedDataFrame", data = gene_metadata)
  cds <- newCellDataSet(expr_mat,
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0,
                        expressionFamily = uninormal())
  
  cds$Pseudotime = cds@phenoData@data[,grep('r_dist_',colnames(cds@phenoData@data))]
  model_formula_str <- "~ Pseudotime"
  
  diff_test_res <- differentialGeneTest(cds,
                                        fullModelFormulaStr = model_formula_str,
                                        cores = 1)
  
  sig_genes <- subset(diff_test_res, qval < 0.05)
  significant_genes_sorted <- sig_genes[order(sig_genes$qval), ]
  
  sig_genes_sorted = row.names(significant_genes_sorted)
  cds_subset = cds[sig_genes_sorted,]
  
  trend_formula = "~sm.ns(Pseudotime, df=3)"
  pseudocount <- 1
  norm_method = c("log")
  scale_max = 3
  scale_min = -3
  cores=1
  pseudotime_grid_df <- data.frame(Pseudotime = seq(min(pData(cds)$Pseudotime),
                                         max(pData(cds_subset)$Pseudotime), length.out = 100))
  m <- genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula,
                       relative_expr = T, new_data = pseudotime_grid_df)
  m = m[!apply(m, 1, sum) == 0, ]
  
  if(norm_method == "log") {
    m = log10(m + pseudocount)
  }
  m = m[!apply(m, 1, sd) == 0, ]
  m = Matrix::t(scale(Matrix::t(m), center = TRUE))
  m = m[is.na(row.names(m)) == FALSE, ]
  m[is.nan(m)] = 0
  m[m > scale_max] = scale_max
  m[m < scale_min] = scale_min
  
  cor_matrix[[j]] <- vector("numeric", nrow(m))

  for(i in 1:nrow(m)){
    cor_test_result <- cor.test(m[i,], pseudotime_grid_df$Pseudotime,method = "pearson")
    cor_matrix[[j]][i] = cor_test_result$estimate
  }
  names(cor_matrix[[j]]) = rownames(m)
}


high_low_pathways = c("HALLMARK_COMPLEMENT","HALLMARK_IL2_STAT5_SIGNALING",
                      "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE",
                      "MP18.Interferon.MHC.II..II.","HALLMARK_APOPTOSIS",
                      "HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_INTERFERON_ALPHA_RESPONSE")
low_high_pathways = c("HALLMARK_E2F_TARGETS","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                      "HALLMARK_NOTCH_SIGNALING","HALLMARK_WNT_BETA_CATENIN_SIGNALING",
                      "MP1..Cell.Cycle...G2.M",  "MP20.MYC",
                      "MP30.PDAC.classical","HALLMARK_MYOGENESIS","HALLMARK_XENOBIOTIC_METABOLISM","Mesenchymal")
all_pathways = c(high_low_pathways,low_high_pathways)

pathway_cor_by_cancer = list()
for(i in c(1,3:8,10)){
  pathway_cor_by_cancer[[i]] = cor_matrix[[i]][which(names(cor_matrix[[i]]) %in% c(high_low_pathways,low_high_pathways))]
}

cancer_names <- c("Bladder", "Brain", "Breast", "ccRCC", "Colorectal", "Gastric", "Liver", "Lung", "OSCC", "Ovarian", "Pancreatic")

names(pathway_cor_by_cancer) = cancer_names[c(1,3:11)]

pathway_cor_mat <- matrix(NA, nrow = length(all_pathways), ncol = length(pathway_cor_by_cancer))
rownames(pathway_cor_mat) <- all_pathways
colnames(pathway_cor_mat) <- cancer_names[c(1,3:11)]

for(i in c(1,3:8,10)){
  cancer_name <- cancer_names[i]
  cor_vector <- pathway_cor_by_cancer[[i]]
  
  cor_vector <- sapply(cor_vector, function(x) if(is.null(x)) NA else x)
  cor_vector <- cor_vector[match(all_pathways, names(cor_vector))]
  cor_vector = unname(cor_vector)
  col_index <- match(cancer_name, colnames(pathway_cor_mat))
  
  pathway_cor_mat[, col_index] <- unlist(unname(cor_vector))
}
pathway_cor_mat = pathway_cor_mat[,!colnames(pathway_cor_mat) %in% c("OSCC","Pancreatic")]
pathway_cor_mat = data.frame(pathway_cor_mat)
pathway_cor_mat$pathway = rownames(pathway_cor_mat)
pathway_cor_mat$category = c(rep("High_Low",length(high_low_pathways)),rep("Low_High",length(low_high_pathways)))

IT_high_low_pathways = read.csv('./data/IT_sig_High_Low.csv')
IT_low_high_pathways = read.csv('./data/IT_sig_Low_High.csv')

IT_high_low_pathways$X = gsub("HALLMARK_","",IT_high_low_pathways$X)
IT_low_high_pathways$X = gsub("HALLMARK_","",IT_low_high_pathways$X)
IT_high_low_pathways$X = gsub("MP18.","",IT_high_low_pathways$X)
IT_low_high_pathways$X = gsub("MP18.","",IT_low_high_pathways$X)

pathway_cor_long_df = data.frame()
for(i in 1:length(high_low_pathways)){
  cancer_high_low = unlist(strsplit(IT_high_low_pathways[i,]$overlap_cancer_types,', '))
  long_format <- melt(pathway_cor_mat[i,colnames(pathway_cor_mat) %in% c(cancer_high_low,'pathway')]
                      , id.vars = "pathway", variable.name = "cancer_type", value.name = "value")
  pathway_cor_long_df = rbind(pathway_cor_long_df,long_format)
}

for(i in 9:nrow(pathway_cor_mat)){
  k = i - 8
  cancer_low_high = unlist(strsplit(IT_low_high_pathways[k,]$overlap_cancer_types,', '))
  long_format <- melt(pathway_cor_mat[i,colnames(pathway_cor_mat) %in% c(cancer_low_high,'pathway')]
                      , id.vars = "pathway", variable.name = "cancer_type", value.name = "value")
  pathway_cor_long_df = rbind(pathway_cor_long_df,long_format)
}

pathway_cor_long_df$pathway = factor(pathway_cor_long_df$pathway,levels = (c(all_pathways)))

pathway_cor_summary_df <- pathway_cor_long_df %>%
  group_by(pathway) %>%
  summarise(average_value = mean(value, na.rm = TRUE))

pathway_cor_summary_df$pathway = gsub('HALLMARK_','',pathway_cor_summary_df$pathway)

pathway_cor_summary_df$category = c(rep("High_Low",length(high_low_pathways)),rep("Low_High",length(low_high_pathways)))


pathway_cor_summary_df = pathway_cor_summary_df[pathway_cor_summary_df$pathway %in% c('INTERFERON_ALPHA_RESPONSE','INTERFERON_GAMMA_RESPONSE',"MP18.Interferon.MHC.II..II.",'INFLAMMATORY_RESPONSE',"APOPTOSIS",
                                                                                      'EPITHELIAL_MESENCHYMAL_TRANSITION','Mesenchymal','MP30.PDAC.classical','MP1..Cell.Cycle...G2.M','MP20.MYC'),]

pathway_cor_summary_df = pathway_cor_summary_df[order(pathway_cor_summary_df$average_value),]
pathway_cor_summary_df$pathway[which(pathway_cor_summary_df=="INTERFERON_ALPHA_RESPONSE")] = "Interferon-alpha response"
pathway_cor_summary_df$pathway[which(pathway_cor_summary_df=="INTERFERON_GAMMA_RESPONSE")] = "Interferon-gamma response"
pathway_cor_summary_df$pathway[which(pathway_cor_summary_df=="MP18.Interferon.MHC.II..II.")] = "MHC-II"
pathway_cor_summary_df$pathway[which(pathway_cor_summary_df=="INFLAMMATORY_RESPONSE")] = "Inflammatory"
pathway_cor_summary_df$pathway[which(pathway_cor_summary_df=="APOPTOSIS")] = "Apoptosis"

pathway_cor_summary_df$pathway[which(pathway_cor_summary_df=="EPITHELIAL_MESENCHYMAL_TRANSITION")] = "EMT"
pathway_cor_summary_df$pathway[which(pathway_cor_summary_df=="Mesenchymal")] = "Mesenchymal"
pathway_cor_summary_df$pathway[which(pathway_cor_summary_df=="MP20.MYC")] = "MYC"
pathway_cor_summary_df$pathway[which(pathway_cor_summary_df=="MP1..Cell.Cycle...G2.M")] = "G2M Cell Cycle"
pathway_cor_summary_df$pathway[which(pathway_cor_summary_df=="MP30.PDAC.classical")] = "KRAS signaling"

pathway_cor_summary_df$pathway = factor(pathway_cor_summary_df$pathway,levels = (pathway_cor_summary_df$pathway))

a= ggplot(pathway_cor_summary_df, aes(x = pathway, y = average_value,fill=category)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_classic() +
  scale_y_reverse() +
  scale_fill_manual(values = c("#F47F72",'#8491B4')) +
  labs(title = "",x="",y = "rho") +
  theme(text = element_text(size = 8))

pdf("./result/Fig4B.pdf",width=20)
print(a)
dev.off()

cancer_type_map <- c(
  "Lung" = "LUAD",
  "Gastric" = "STAD",
  "CCRCC" = "KIRC",
  "Breast" = "BRCA",
  "CSCC" = "CSCC",
  "OSCC" = "OSCC",
  "Bladder" = "BLCA",
  "Pancreatic" = "PAAD",
  "Colorectal" = "CRC",
  "Ovarian" = "OVCA",
  "Liver" = "LIHC"
)


target_n <- 8

cor_mat_cancer_prop <- pathway_cor_long_df %>%
  mutate(
    cancer_type = as.character(cancer_type),
    cancer_type = ifelse(cancer_type == "ccRCC", "CCRCC", cancer_type),
    cancer_type = unname(cancer_type_map[cancer_type]),
    cancer_type = ifelse(is.na(cancer_type), "", cancer_type),
    value = 1/target_n
  ) %>%
  group_by(pathway) %>%
  group_modify(~{
    n_real <- nrow(.x)
    n_dummy <- max(0, target_n - n_real)
    if (n_dummy == 0) return(.x)
    bind_rows(.x, tibble(pathway = unique(.x$pathway), cancer_type = "", value = 1/target_n)[rep(1, n_dummy), ])
  }) %>%
  ungroup()

colors <- c('white',"#6B3F98",brewer.pal(8, 'Paired'))
cancer_types <- c("","OVCA","KIRC","BLCA","empty","CRC","BRCA","LIHC","STAD","LUAD")
color_mapping <- setNames(colors[1:length(cancer_types)], cancer_types)

cor_mat_cancer_prop$cancer_type = factor(cor_mat_cancer_prop$cancer_type,levels = rev(cancer_types))

ten_pathways = c("HALLMARK_INTERFERON_ALPHA_RESPONSE","MP18.Interferon.MHC.II..II.",
                 "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE",
                 "HALLMARK_APOPTOSIS","MP20.MYC",
                 "MP1..Cell.Cycle...G2.M","MP30.PDAC.classical",
                 "Mesenchymal","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

cor_mat_cancer_prop = cor_mat_cancer_prop[cor_mat_cancer_prop$pathway %in% ten_pathways,]

cor_mat_cancer_prop$pathway <- as.character(cor_mat_cancer_prop$pathway)

cor_mat_cancer_prop$pathway[which(cor_mat_cancer_prop=="HALLMARK_INTERFERON_ALPHA_RESPONSE")] = "Interferon-alpha response"
cor_mat_cancer_prop$pathway[which(cor_mat_cancer_prop=="HALLMARK_INTERFERON_GAMMA_RESPONSE")] = "Interferon-gamma response"
cor_mat_cancer_prop$pathway[which(cor_mat_cancer_prop=="MP18.Interferon.MHC.II..II.")] = "MHC-II"
cor_mat_cancer_prop$pathway[which(cor_mat_cancer_prop=="HALLMARK_INFLAMMATORY_RESPONSE")] = "Inflammatory"
cor_mat_cancer_prop$pathway[which(cor_mat_cancer_prop=="HALLMARK_APOPTOSIS")] = "Apoptosis"

cor_mat_cancer_prop$pathway[which(cor_mat_cancer_prop=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")] = "EMT"
cor_mat_cancer_prop$pathway[which(cor_mat_cancer_prop=="Mesenchymal")] = "Mesenchymal"
cor_mat_cancer_prop$pathway[which(cor_mat_cancer_prop=="MP20.MYC")] = "MYC"
cor_mat_cancer_prop$pathway[which(cor_mat_cancer_prop=="MP1..Cell.Cycle...G2.M")] = "G2M Cell Cycle"
cor_mat_cancer_prop$pathway[which(cor_mat_cancer_prop=="MP30.PDAC.classical")] = "KRAS signaling"

for(l in 1:length(unique(cor_mat_cancer_prop$pathway))){
  cor_mat_cancer_prop_sub = cor_mat_cancer_prop[cor_mat_cancer_prop$pathway == unique(cor_mat_cancer_prop$pathway)[l],]
  cor_mat_cancer_prop_sub = cor_mat_cancer_prop_sub %>% mutate(lab.ypos = cumsum(value) - 0.5*value)
  
  b = ggplot(cor_mat_cancer_prop_sub, aes(x = "", y = value, fill = cancer_type)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0,direction=1) +
    scale_fill_manual(values = color_mapping) +
    theme_void()
  
  pdf(paste0("./result/Fig4B_pie_",unique(cor_mat_cancer_prop$pathway)[l],".pdf"))
  print(b)
  dev.off()
}

rm(list=ls())

#### Figure 4C ####
ten_pathways = c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                 "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                 "MP18.Interferon.MHC.II..II.",
                 "HALLMARK_INFLAMMATORY_RESPONSE",
                 "HALLMARK_APOPTOSIS",
                 "MP1..Cell.Cycle...G2.M",
                 "MP20.MYC",
                 "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                 "MP30.PDAC.classical",
                 "Mesenchymal")

add.flag <- function(pheatmap,kept.labels) {
  
  heatmap <- pheatmap$gtable
  
  
  leg_idx <- which(heatmap$layout$name == "annotation_legend")
  
  leg <- heatmap$grobs[[leg_idx]]
  distance_map = c(`0`="Proximal", `1`="Distal")
  cluster_map  = c(`1`="High_Low", `2`="Low_High")
  
  for (i in seq_along(leg$children)) {
    g <- leg$children[[i]]
    
    if (!inherits(g, "text")) next
    
    
    if (length(g$label) == 1 && identical(g$label, "Distance")) {
      g$label <- "Distance"
    }
    
    
    if (length(g$label) == 2 && all(g$label %in% c("0","1")) && !all(g$label %in% c("1","2"))) {
      g$label <- unname(distance_map[g$label])
    }
    
    
    if (length(g$label) == 1 && identical(g$label, "Cluster")) {
      g$label <- "Cluster"
    }
    
    
    if (length(g$label) == 2 && all(g$label %in% c("1","2"))) {
      new <- unname(cluster_map[g$label])
      new <- new[order(match(new, c(cluster_map[["1"]], cluster_map[["2"]])))]
      g$label <- new
    }
    
    leg$children[[i]] <- g
  }
  
  heatmap$grobs[[leg_idx]] <- leg
  
  
  rn_idx <- which(heatmap$layout$name == "row_names")
  new.label <- heatmap$grobs[[rn_idx]]
  
  paths <- pseudotime_heatmap_obj$tree_row$labels
  cl_id <- cutree(pseudotime_heatmap_obj$tree_row, k = 2)[paths]
  cluster_map <- c(`1`="High_Low", `2`="Low_High")
  cl_name <- setNames(cluster_map[as.character(cl_id)], paths)
  
  high_low_vec <- c(
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "MP18.Interferon.MHC.II..II.",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_APOPTOSIS"
  )
  
  low_high_vec <- c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "Mesenchymal",
    "MP20.MYC",
    "MP1..Cell.Cycle...G2.M",
    "MP30.PDAC.classical"
  )
  
  
  kept_vec <- c(
    high_low_vec[cl_name[high_low_vec] == "High_Low"],
    low_high_vec[cl_name[low_high_vec] == "Low_High"]
  )
  
  
  kept.labels <- intersect(kept.labels, kept_vec)
  
  new.label$label <- ifelse(new.label$label %in% kept.labels, new.label$label, "")
  new.label$label <- ifelse(new.label$label == "HALLMARK_INTERFERON_ALPHA_RESPONSE", "Interferon-alpha response", new.label$label)
  new.label$label <- ifelse(new.label$label == "HALLMARK_INTERFERON_GAMMA_RESPONSE", "Interferon-gamma response", new.label$label)
  new.label$label <- ifelse(new.label$label == "MP18.Interferon.MHC.II..II.", "MHC-II", new.label$label)
  new.label$label <- ifelse(new.label$label == "HALLMARK_INFLAMMATORY_RESPONSE", "Inflammatory", new.label$label)
  new.label$label <- ifelse(new.label$label == "HALLMARK_APOPTOSIS", "Apoptosis", new.label$label)
  
  new.label$label <- ifelse(new.label$label == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "EMT", new.label$label)
  new.label$label <- ifelse(new.label$label == "Mesenchymal", "Mesenchymal", new.label$label)
  new.label$label <- ifelse(new.label$label == "MP20.MYC", "MYC", new.label$label)
  new.label$label <- ifelse(new.label$label == "MP1..Cell.Cycle...G2.M", "G2M Cell Cycle", new.label$label)
  new.label$label <- ifelse(new.label$label == "MP30.PDAC.classical", "KRAS signaling", new.label$label)
  
  heatmap$grobs[[rn_idx]] <- new.label
  
  
  grid::grid.newpage()
  grid::grid.draw(heatmap)
  invisible(heatmap)
}




for(j in c(8,3,7,6)){
  cancer_index = j
  cancer_names <- c("Bladder", "Brain", "Breast", "ccRCC", "Colorectal", "Gastric", "Liver", "Lung", "OSCC", "Ovarian", "Pancreatic")
  
  tls_object = "IT"
  
  load(paste0("./data/","ST_",cancer_index,"_","combined_sig_",tls_object,"_objects.RData"))
  sig_variable_name <- paste0("combined_sig_", tls_object)
  meta_data_variable_name <- paste0("combined_meta_data_", tls_object)
  
  expr_mat <- get(sig_variable_name)
  cell_metadata_df <- get(meta_data_variable_name)
  
  gene_metadata <- data.frame(gene_short_name = rownames(expr_mat), row.names = rownames(expr_mat))
  cell_metadata <- cell_metadata_df
  pd <- new("AnnotatedDataFrame", data = cell_metadata)
  fd <- new("AnnotatedDataFrame", data = gene_metadata)
  cds <- newCellDataSet(expr_mat,
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0,
                        expressionFamily = uninormal())
  
  cds$Pseudotime = cds@phenoData@data[,grep('r_dist_',colnames(cds@phenoData@data))]
  model_formula_str <- "~ Pseudotime"
  
  diff_test_res <- differentialGeneTest(cds,
                                        fullModelFormulaStr = model_formula_str,
                                        cores = 1)
  
  sig_genes <- subset(diff_test_res, qval < 0.05)
  significant_genes_sorted <- sig_genes[order(sig_genes$qval), ]
  
  sig_genes_sorted = row.names(significant_genes_sorted)
  cds_subset = cds[sig_genes_sorted,]
  
  trend_formula = "~sm.ns(Pseudotime, df=3)"
  pseudocount <- 1
  norm_method = c("log")
  scale_max = 3
  scale_min = -3
  cores=1
  pseudotime_grid_df <- data.frame(Pseudotime = seq(min(pData(cds)$Pseudotime),
                                         max(pData(cds_subset)$Pseudotime), length.out = 100))
  m <- genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula,
                       relative_expr = T, new_data = pseudotime_grid_df)
  m = m[!apply(m, 1, sum) == 0, ]
  
  if(norm_method == "log") {
    m = log10(m + pseudocount)
  }
  m = m[!apply(m, 1, sd) == 0, ]
  m = Matrix::t(scale(Matrix::t(m), center = TRUE))
  m = m[is.na(row.names(m)) == FALSE, ]
  m[is.nan(m)] = 0
  m[m > scale_max] = scale_max
  m[m < scale_min] = scale_min
  
  diff_test_res_sig = read.csv(file=paste0('./data/',"ST_",cancer_index,"_", tls_object,'_sig_trend.csv'),row.names = 1)
  
  pseudotime_grid_df2 = pseudotime_grid_df
  colnames(pseudotime_grid_df2) = 'Distance'
  
  cancer_map <- c(
    Gastric = "STAD",
    Liver   = "LIHC",
    Breast  = "BRCA",
    Lung    = "LUAD"
  )
  file_name <- cancer_map[cancer_names[j]]
  
  if(j == 3){
    breast_remove = c("MP6.Hypoxia","Hypoxia","HALLMARK_CHOLESTEROL_HOMEOSTASIS","HALLMARK_MTORC1_SIGNALING","HALLMARK_MYOGENESIS","AC","HALLMARK_UNFOLDED_PROTEIN_RESPONSE","HALLMARK_PEROXISOME",
                      "HALLMARK_P53_PATHWAY","HALLMARK_E2F_TARGETS","HALLMARK_HYPOXIA","Alveolar")
    genes_breast = rownames(diff_test_res_sig[diff_test_res_sig$trend %in% c('Low_High','High_Low'),])
    genes_breast = genes_breast[!genes_breast %in% breast_remove]
    pseudotime_heatmap_obj = plot_pseudotime_heatmap(cds[genes_breast,],
                                 num_clusters = 2,
                                 add_annotation_col = pseudotime_grid_df2,
                                 cores = 1,
                                 show_rownames = T,return_heatmap=T)
    
    pdf(paste0('./result/Fig4C_',file_name,'.pdf'))
    add.flag(pseudotime_heatmap_obj,kept.labels = ten_pathways)
    dev.off()
  }
  
  if(j == 8){
    pathway_remove = c("HALLMARK_GLYCOLYSIS","MP10.Protein.maturation","MP7.Stress..in.vitro.","HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_P53_PATHWAY","Oxphos","NPC","MP21.Respiration",
                       "HALLMARK_KRAS_SIGNALING_UP","MP6.Hypoxia")
    genes_remove = rownames(diff_test_res_sig[diff_test_res_sig$trend %in% c('Low_High','High_Low'),])
    genes_remove = genes_remove[!genes_remove %in% pathway_remove]
    
    pseudotime_heatmap_obj = plot_pseudotime_heatmap(cds[genes_remove,],
                                 num_clusters = 2,
                                 add_annotation_col = pseudotime_grid_df2,
                                 cores = 1,
                                 show_rownames = T,return_heatmap=T)
    
    pdf(paste0('./result/Fig4C_',file_name,'.pdf'))
    add.flag(pseudotime_heatmap_obj,kept.labels = ten_pathways)
    dev.off()
  }
  
  if(j == 7){
    pathway_remove = c("MP21.Respiration","HALLMARK_GLYCOLYSIS","MP10.Protein.maturation","MP7.Stress..in.vitro.","HALLMARK_OXIDATIVE_PHOSPHORYLATION","MP36.IG","HALLMARK_P53_PATHWAY",
                       "HALLMARK_ESTROGEN_RESPONSE_EARLY","Hypoxia")
    genes_pathway = rownames(diff_test_res_sig[diff_test_res_sig$trend %in% c('Low_High','High_Low'),])
    genes_pathway = genes_pathway[!genes_pathway %in% pathway_remove]
    pseudotime_heatmap_obj = plot_pseudotime_heatmap(cds[genes_pathway,],
                                 num_clusters = 2,
                                 add_annotation_col = pseudotime_grid_df2,
                                 cores = 1,
                                 show_rownames = T,return_heatmap=T)
    
    pdf(paste0('./result/Fig4C_',file_name,'.pdf'))
    add.flag(pseudotime_heatmap_obj,kept.labels = ten_pathways)
    dev.off()
  }
  
  if(j == 6){
    pathway_remove = c("HALLMARK_UNFOLDED_PROTEIN_RESPONSE","MP29.NPC.OPC","MP26.NPC.Glioma","OPC","HALLMARK_E2F_TARGETS","MP2..Cell.Cycle...G1.S","HALLMARK_UV_RESPONSE_UP",
                       "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY","Cycle","HALLMARK_WNT_BETA_CATENIN_SIGNALING","Ciliated","MP33.RBCs",
                       "Alveolar","Basal","MP32.Skin.pigmentation","pEMT","HALLMARK_MITOTIC_SPINDLE","MP6.Hypoxia","Hypoxia","HALLMARK_HEDGEHOG_SIGNALING","MP27.Oligo.Progenitor",
                       "MP4..Chromatin","HALLMARK_CHOLESTEROL_HOMEOSTASIS","HALLMARK_MYOGENESIS","MP25.Astrocytes")
    genes_pathway = rownames(diff_test_res_sig[diff_test_res_sig$trend %in% c('Low_High','High_Low'),])
    genes_pathway = genes_pathway[!genes_pathway %in% pathway_remove]
    pseudotime_heatmap_obj = plot_pseudotime_heatmap(cds[genes_pathway,],
                                 num_clusters = 2,
                                 add_annotation_col = pseudotime_grid_df2,
                                 cores = 1,
                                 show_rownames = T,return_heatmap=T)
    
    
    
    pdf(paste0('./result/Fig4C_',file_name,'.pdf'))
    add.flag(pseudotime_heatmap_obj,kept.labels = ten_pathways)
    dev.off()
  }
}

rm(list=ls())

#### Figure 4D ####
IT_high_low_pathways = read.csv('./data/IT_sig_High_Low.csv')
IT_low_high_pathways = read.csv('./data/IT_sig_Low_High.csv')

IT_high_low_pathways = IT_high_low_pathways[IT_high_low_pathways$overlap_cancer_type_cnt > 3,1:3]
IT_low_high_pathways = IT_low_high_pathways[IT_low_high_pathways$overlap_cancer_type_cnt > 3,1:3]

high_low_pathways = c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                      "HALLMARK_INTERFERON_GAMMA_RESPONSE","MP18.Interferon.MHC.II..II.",
                      "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_APOPTOSIS")
low_high_pathways = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                      "MP1..Cell.Cycle...G2.M",  "MP20.MYC",
                      "MP30.PDAC.classical","Mesenchymal")

IT_high_low_pathways = IT_high_low_pathways[IT_high_low_pathways$X %in% high_low_pathways,]
IT_low_high_pathways = IT_low_high_pathways[IT_low_high_pathways$X %in% low_high_pathways,]
all_pathways = c(IT_high_low_pathways$X,IT_low_high_pathways$X)

cancer_names <- c("Bladder", "Brain", "Breast", "ccRCC", "Colorectal", "Gastric", "Liver", "Lung", "OSCC", "Ovarian", "Pancreatic")

colors <- c("#6B3F98",brewer.pal(8, 'Paired'))
cancer_types <- c("OVCA","KIRC","BLCA","","CRC","BRCA","LIHC","STAD","LUAD")
color_mapping <- setNames(colors[1:length(cancer_types)], cancer_types)

pathway_curve_df3 = list()
for(i in 1:length(all_pathways)){
  if(i <= length(IT_high_low_pathways$X)){
    selected_cancer_names = unlist(strsplit(IT_high_low_pathways[IT_high_low_pathways$X == all_pathways[i],'overlap_cancer_types'],", "))
    selected_cancer_idx = which(cancer_names %in% selected_cancer_names)
  } else{
    selected_cancer_names = unlist(strsplit(IT_low_high_pathways[IT_low_high_pathways$X == all_pathways[i],'overlap_cancer_types'],", "))
    selected_cancer_idx = which(cancer_names %in% selected_cancer_names)
  }
  
  pathway_curve_df = data.frame()
  for(j in selected_cancer_idx){
    cancer_index = j
    cancer_names <- c("Bladder", "Brain", "Breast", "ccRCC", "Colorectal", "Gastric", "Liver", "Lung", "OSCC", "Ovarian", "Pancreatic")
    
    tls_object = "IT"
    
    load(paste0("./data/","ST_",cancer_index,"_","combined_sig_",tls_object,"_objects.RData"))
    sig_variable_name <- paste0("combined_sig_", tls_object)
    meta_data_variable_name <- paste0("combined_meta_data_", tls_object)
    
    expr_mat <- get(sig_variable_name)
    cell_metadata_df <- get(meta_data_variable_name)
    
    gene_metadata <- data.frame(gene_short_name = rownames(expr_mat), row.names = rownames(expr_mat))
    cell_metadata <- cell_metadata_df
    pd <- new("AnnotatedDataFrame", data = cell_metadata)
    fd <- new("AnnotatedDataFrame", data = gene_metadata)
    cds <- newCellDataSet(expr_mat,
                          phenoData = pd,
                          featureData = fd,
                          lowerDetectionLimit = 0,
                          expressionFamily = uninormal())
    
    cds$Pseudotime = cds@phenoData@data[,grep('r_dist_',colnames(cds@phenoData@data))]
    model_formula_str <- "~ Pseudotime"
    
    diff_test_res <- differentialGeneTest(cds,
                                          fullModelFormulaStr = model_formula_str,
                                          cores = 1)
    
    sig_genes <- subset(diff_test_res, qval < 0.05)
    significant_genes_sorted <- sig_genes[order(sig_genes$qval), ]
    
    sig_genes_sorted = row.names(significant_genes_sorted)
    cds_subset = cds[sig_genes_sorted,]
    
    trend_formula = "~sm.ns(Pseudotime, df=3)"
    pseudocount <- 1
    norm_method = c("log")
    scale_max = 3
    scale_min = -3
    cores=1
    pseudotime_grid_df <- data.frame(Pseudotime = seq(min(pData(cds)$Pseudotime),
                                           max(pData(cds_subset)$Pseudotime), length.out = 100))
    m <- genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula,
                         relative_expr = T, new_data = pseudotime_grid_df)
    m = m[!apply(m, 1, sum) == 0, ]
    
    if(norm_method == "log") {
      m = log10(m + pseudocount)
    }
    m = m[!apply(m, 1, sd) == 0, ]
    m = Matrix::t(scale(Matrix::t(m), center = TRUE))
    m = m[is.na(row.names(m)) == FALSE, ]
    m[is.nan(m)] = 0
    m[m > scale_max] = scale_max
    m[m < scale_min] = scale_min
    
    m2 = m[rownames(m) %in% all_pathways[i],]
    df = cbind(m2,pseudotime_grid_df,cancer_names[j])
    colnames(df) = c("Value","Distance","Cancer_type")
    pathway_curve_df = rbind(pathway_curve_df,df)
  }
  
  cancer_type_map <- c(
    "Lung" = "LUAD",
    "Gastric" = "STAD",
    "CCRCC" = "KIRC",
    "Breast" = "BRCA",
    "CSCC" = "CSCC",
    "OSCC" = "OSCC",
    "Bladder" = "BLCA",
    "Pancreatic" = "PAAD",
    "Colorectal" = "CRC",
    "Ovarian" = "OVCA",
    "Liver" = "LIHC"
  )
  pathway_curve_df2 = pathway_curve_df
  pathway_curve_df2$Cancer_type = ifelse(pathway_curve_df2$Cancer_type == 'ccRCC',"CCRCC",pathway_curve_df2$Cancer_type)
  pathway_curve_df2$Cancer_type <- unname(cancer_type_map[pathway_curve_df2$Cancer_type])
  pathway_curve_df2$Cancer_type = factor(pathway_curve_df2$Cancer_type,levels = rev(cancer_types[cancer_types!='']))
  pathway_curve_df3[[i]] = pathway_curve_df2
}



all_pathways_rename = all_pathways
all_pathways_rename = gsub('HALLMARK_','',all_pathways_rename)
all_pathways_rename = gsub('_',' ',all_pathways_rename)

all_pathways_rename[which(all_pathways_rename=='MP20.MYC')] = 'MYC'
all_pathways_rename[which(all_pathways_rename=='MP1..Cell.Cycle...G2.M')] = 'G2M Cell Cycle'
all_pathways_rename[which(all_pathways_rename=='EPITHELIAL MESENCHYMAL TRANSITION')] = 'EMT'
all_pathways_rename[which(all_pathways_rename=="MP30.PDAC.classical")] = 'KRAS Signaling'

all_pathways_rename[which(all_pathways_rename=='INFLAMMATORY RESPONSE')] = 'Inflammatory response'
all_pathways_rename[which(all_pathways_rename=='APOPTOSIS')] = 'Apoptosis'
all_pathways_rename[which(all_pathways_rename=='MP18.Interferon.MHC.II..II.')] = 'MHC-II'
all_pathways_rename[which(all_pathways_rename=="INTERFERON ALPHA RESPONSE")] = 'Interferon-alpha response'
all_pathways_rename[which(all_pathways_rename=="INTERFERON GAMMA RESPONSE")] = 'Interferon-gamma response'

idx_myc <- which(all_pathways_rename == "MYC")
pathway_curve_df3[[idx_myc]] <-
  pathway_curve_df3[[idx_myc]][pathway_curve_df3[[idx_myc]]$Cancer_type != "BLCA", ]


for(i in 1:10){
  a <- ggplot(pathway_curve_df3[[i]], aes(x = Distance, y = Value, color = Cancer_type)) +
    geom_line(linewidth = 0.6) +
    theme_classic(base_size = 6) +
    scale_colour_manual(values = color_mapping) +
    labs(title = all_pathways_rename[i],x='',y='') +
    scale_x_continuous(breaks = c(0, 1),labels = c("Proximal", "Distal"),expand = expansion(mult = c(0.05, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    theme(
      text = element_text(size = 8, colour = "black"),
      axis.text  = element_text(size = 8, colour = "black"),
      axis.title = element_text(size = 8, colour = "black"),
      plot.title = element_text(size = 8, colour = "black"),
      plot.subtitle = element_text(size = 8, colour = "black"),
      legend.text  = element_text(size = 8, colour = "black"),
      legend.title = element_blank()
      
    )
  g <- ggplotGrob(a)
  panel_idx <- g$layout[g$layout$name == "panel", ]
  g$widths[panel_idx$l]  <- unit(120, "pt")
  g$heights[panel_idx$t] <- unit(120, "pt")
  
  
  w_in <- convertWidth(sum(g$widths),  "in", valueOnly = TRUE)
  h_in <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  
  pdf(sprintf('./result/Fig4D_%s.pdf', all_pathways_rename[i]),
      width = w_in, height = h_in, useDingbats = FALSE)
  grid.newpage(); grid.draw(g)
  dev.off()
  
}

rm(list=ls())

#### Figure 4E ####
# Before running the script below, first reconstruct the two .RData files
# by concatenating their split parts in the terminal.
#cat ./data/ST_6_combined_IT_objects.part_* > \
#./data/ST_6_combined_IT_objects.RData
#cat ./ST_8_combined_IT_objects.part_* > \
#./ST_8_combined_IT_objects.RData
IT_high_low_genes <- read.csv('./data/IT_High_Low.csv')

cancer_names <- c("Bladder", "Brain", "Breast", "ccRCC", "Colorectal", "Gastric", "Liver", "Lung", "OSCC", "Ovarian", "Pancreatic")

high_low_genes = c("TNFRSF14","CXCL16","CXCL10","CXCL11")

for(pathway in high_low_genes) {
  selected_cancer_names <- unlist(strsplit(IT_high_low_genes[IT_high_low_genes$X == pathway, 'overlap_cancer_types'], ", "))
  selected_cancer_idx <- which(cancer_names %in% selected_cancer_names)
  
  pathway_curve_df <- data.frame()
  
  for(j in selected_cancer_idx) {
    cancer_index <- j
    
    cancer_names <- c("Bladder", "Brain", "Breast", "ccRCC", "Colorectal", "Gastric", "Liver", "Lung", "OSCC", "Ovarian", "Pancreatic")
    
    tls_object <- "IT"
    
    load(paste0("./data/", "ST_", cancer_index, "_", "combined_", tls_object, "_objects.RData"))
    sig_variable_name <- paste0("combined_counts_", tls_object)
    meta_data_variable_name <- paste0("combined_meta_data_", tls_object)
    
    expr_mat <- get(sig_variable_name)
    cell_metadata_df <- get(meta_data_variable_name)
    
    gene_metadata <- data.frame(gene_short_name = rownames(expr_mat), row.names = rownames(expr_mat))
    cell_metadata <- cell_metadata_df
    pd <- new("AnnotatedDataFrame", data = cell_metadata)
    fd <- new("AnnotatedDataFrame", data = gene_metadata)
    cds <- newCellDataSet(expr_mat,
                          phenoData = pd,
                          featureData = fd,
                          lowerDetectionLimit = 0,
                          expressionFamily = uninormal())
    
    cds$Pseudotime <- cds@phenoData@data[, grep('r_dist_', colnames(cds@phenoData@data))]
    model_formula_str <- "~ Pseudotime"
    
    sig_genes_sorted = high_low_genes
    cds_subset <- cds[sig_genes_sorted, ]
    
    trend_formula <- "~ poly(Pseudotime, 3)"
    pseudocount <- 1
    norm_method <- "log"
    scale_max <- 3
    scale_min <- -3
    cores <- 1
    pseudotime_grid_df <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime),
                                           max(pData(cds_subset)$Pseudotime), length.out = 100))
    
    m <- genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula,
                         relative_expr = TRUE, new_data = pseudotime_grid_df)
    m <- m[!apply(m, 1, sum) == 0, ]
    
    if(norm_method == "log") {
      m <- log10(m + pseudocount)
    }
    m <- m[!apply(m, 1, sd) == 0, ]
    m <- Matrix::t(scale(Matrix::t(m), center = TRUE))
    m <- m[!is.na(rownames(m)), ]
    m[is.nan(m)] <- 0
    m[m > scale_max] <- scale_max
    m[m < scale_min] <- scale_min
    
    m2 <- m[rownames(m) %in% pathway, ]
    df <- cbind(m2, pseudotime_grid_df, cancer_names[j])
    colnames(df) <- c("Value", "Distance", "Cancer_type")
    pathway_curve_df <- rbind(pathway_curve_df, df)
  }
  
  cancer_type_map <- c(
    "Lung" = "LUAD",
    "Gastric" = "STAD",
    "CCRCC" = "KIRC",
    "Breast" = "BRCA",
    "CSCC" = "CSCC",
    "OSCC" = "OSCC",
    "Bladder" = "BLCA",
    "Pancreatic" = "PAAD",
    "Colorectal" = "CRC",
    "Ovarian" = "OVCA",
    "Liver" = "LIHC"
  )
  pathway_curve_df2 <- pathway_curve_df
  pathway_curve_df2$Cancer_type <- ifelse(pathway_curve_df2$Cancer_type == 'ccRCC', "CCRCC", pathway_curve_df2$Cancer_type)
  pathway_curve_df2$Cancer_type <- unname(cancer_type_map[pathway_curve_df2$Cancer_type])
  
  colors <- c("#6B3F98", brewer.pal(10, 'Paired'))
  cancer_types <- c("OVCA", "KIRC", "BLCA", "", "CRC", "BRCA", "LIHC", "STAD", "LUAD","OSCC","PAAD")
  color_mapping <- setNames(colors[1:length(cancer_types)], cancer_types)
  pathway_curve_df2$Cancer_type <- factor(pathway_curve_df2$Cancer_type, levels = rev(cancer_types[cancer_types != '']))
  
  pathway_curve_df3 <- pathway_curve_df2[pathway_curve_df2$Distance < 0.75, ]
  
  if (pathway == "CXCL16") {
    pathway_curve_df3 = pathway_curve_df3[!pathway_curve_df3$Cancer_type %in% c('LUAD','LIHC'),]
  }
  
  if (pathway == "CXCL11") {
    pathway_curve_df3 = pathway_curve_df3[!pathway_curve_df3$Cancer_type %in% c('STAD'),]
  }
  
  a <- ggplot(pathway_curve_df3, aes(x = Distance, y = Value, color = Cancer_type)) +
    geom_line(linewidth = 0.6) +
    theme_classic(base_size = 6) +
    scale_colour_manual(values = color_mapping) +
    labs(title = pathway,x='',y='') +
    scale_x_continuous(breaks = c(0, 0.75),labels = c("Proximal", "Distal"),expand = expansion(mult = c(0.05, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    theme(
      text = element_text(size = 8, colour = "black"),
      axis.text  = element_text(size = 8, colour = "black"),
      axis.title = element_text(size = 8, colour = "black"),
      plot.title = element_text(size = 8, colour = "black"),
      plot.subtitle = element_text(size = 8, colour = "black"),
      legend.text  = element_text(size = 8, colour = "black"),
      legend.title = element_blank()
      
    )
  
  if (pathway == "CXCL16") {
    a <- a + scale_y_continuous(
      breaks = c(-1.5, -0.75, 0, 0.75, 1.5)
    )
  }
  
  g <- ggplotGrob(a)
  panel_idx <- g$layout[g$layout$name == "panel", ]
  g$widths[panel_idx$l]  <- unit(100, "pt")
  g$heights[panel_idx$t] <- unit(100, "pt")
  
  
  w_in <- convertWidth(sum(g$widths),  "in", valueOnly = TRUE)
  h_in <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  
  pdf(sprintf('./result/Fig4E_%s.pdf', pathway),
      width = w_in, height = h_in, useDingbats = FALSE)
  grid.newpage(); grid.draw(g)
  dev.off()
}


rm(list=ls())

