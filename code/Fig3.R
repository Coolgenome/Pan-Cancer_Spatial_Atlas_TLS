### Figure 3 ###
### The codes are separated by Figures ###

### load packages ###
library(viridis)
library(ComplexHeatmap)
library(circlize)

#### Figure 3A ####
gene_mat <- readRDS("./data/TLS_gene_expressions.rds")

tls_gene_expr_mat = gene_mat
rownames(tls_gene_expr_mat) = c("E_TLS","P_TLS","S_TLS")
marker_annotation_df = read.csv("./data/fig3_markers.csv")
major_types = unique(marker_annotation_df$Major_type)

B_category_orders <- c("B","Bn","Breg","Bm","GCB","PC")
CAF_Macrophage_DC_category_orders <- c("TCR signaling","DC","Macrophage","Endothelial","CAF","Fibroblast","myCAF","Antigen presentation")
T_category_orders <- c("T","Tn","Tcm","Teff","Tfh","T cell exhaustion","Treg",'Th17','FDC',"Interferon response","Stress response")
Others_category_orders <- c("Cytokines and Chemokines","Proliferation")

pdf_height = 4
pdf_width = 30


for(j in 1:length(major_types)){
  major_type_markers_df = marker_annotation_df[marker_annotation_df$Major_type == major_types[j],]
  celltype_order =  get(paste0(major_types[j],"_category_orders"))
  ht = list()
  for(i in 1:length(celltype_order)){
    major_type_markers_df_subset = major_type_markers_df[major_type_markers_df$celltype == celltype_order[i],]
    
    rownames(major_type_markers_df_subset) = major_type_markers_df_subset$gene
    major_type_markers_df2 = major_type_markers_df_subset[major_type_markers_df_subset$gene %in% colnames(tls_gene_expr_mat),]
    
    rownames(major_type_markers_df2) = major_type_markers_df2$gene
    tls_gene_expr_mat_subset = tls_gene_expr_mat[,major_type_markers_df2$gene]
    
    column_split_mat <- as.matrix(major_type_markers_df2[colnames(tls_gene_expr_mat_subset), "celltype"])
    
    celltype_levels <- unique(major_type_markers_df2$celltype)
    column_split_df <- as.data.frame(column_split_mat)
    column_split_df$V1 <- factor(column_split_df$V1, levels = celltype_levels)
    if (nrow(major_type_markers_df_subset) != 1) {
      ht[[i]] <- Heatmap(tls_gene_expr_mat_subset,
                         cluster_columns = T,
                         cluster_rows =  F,
                         show_column_dend = F,
                         row_names_side = "left",
                         name = "Gene Expression",
                         column_names_gp = gpar(fontsize = 27, rot = 45),
                         row_names_gp = gpar(fontsize = 27),
                         column_split = column_split_df,
                         column_gap = unit(1.5, "mm"),
                         col = colorRampPalette(c('blue', 'yellow'))(13),
                         column_names_rot = 45)
    } else {
      ht[[i]] <- Heatmap(tls_gene_expr_mat_subset,
                         cluster_columns = T,
                         cluster_rows = F,
                         show_column_dend = F,
                         row_names_side = "left",
                         name = "Gene Expression",
                         column_names_gp = gpar(fontsize = 27, rot = 45),
                         row_names_gp = gpar(fontsize = 27),
                         column_gap = unit(1.5, "mm"),
                         col = colorRampPalette(c('blue', 'yellow'))(13),
                         column_names_rot = 45)
    }
  }
  if (major_types[j] == "T") {
    combined_heatmap <- ht[[1]] + ht[[2]] + ht[[3]] + ht[[4]]
    pdf(paste0("./result/Fig3A_", major_types[j], "_1.pdf"),
        width = pdf_width, height = pdf_height)
    draw(combined_heatmap)
    dev.off()
    
    combined_heatmap <- ht[[5]] + ht[[6]] + ht[[7]] + ht[[8]] + ht[[9]]
    pdf(paste0("./result/Fig3A_", major_types[j], "_2.pdf"),
        width = pdf_width, height = pdf_height)
    draw(combined_heatmap)
    dev.off()
  } else if (major_types[j] %in% c("B","Others")) {
    combined_heatmap <- Reduce(`+`, ht)
    
    pdf(paste0("./result/Fig3A_", major_types[j], ".pdf"),
        width = pdf_width, height = pdf_height)
    draw(combined_heatmap)
    dev.off()
  }
}

rm(list = ls())
