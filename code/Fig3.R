### Figure 3 ###
### The codes are separated by Figures ###

### load packages ###
library(viridis)
library(ComplexHeatmap)
library(circlize)

#### Figure 3A ####
gene_mat <- readRDS("./data/TLS_gene_expressions.rds")

gene_mat2 = gene_mat
rownames(gene_mat2) = c("E_TLS","P_TLS","S_TLS")
gene_category_all = read.csv("./data/fig3_markers.csv")
major_category = unique(gene_category_all$Major_type)

B_category_orders <- c("B","Bn","Breg","Bm","GCB","PC")
CAF_Macrophage_DC_category_orders <- c("TCR signaling","DC","Macrophage","Endothelial","CAF","Fibroblast","myCAF","Antigen presentation")
T_category_orders <- c("T","Tn","Tcm","Teff","Tfh","T cell exhaustion","Treg",'Th17','FDC',"Interferon response","Stress response")
Others_category_orders <- c("Cytokines and Chemokines","Proliferation")

pdf_height = 4
T_height = 80
rest_height = 30

for(j in 1:length(major_category)){
  gene_category = gene_category_all[gene_category_all$Major_type == major_category[j],]
  celltype_in_order =  get(paste0(major_category[j],"_category_orders"))
  pdf_width = ifelse(major_category[j] != "T",rest_height,rest_height)
  
  ht = list()
  for(i in 1:length(celltype_in_order)){
    gene_category_subset = gene_category[gene_category$celltype == celltype_in_order[i],]
    
    rownames(gene_category_subset) = gene_category_subset$gene
    gene_category2 = gene_category_subset[gene_category_subset$gene %in% colnames(gene_mat2),]
    
    rownames(gene_category2) = gene_category2$gene
    gene_mat2_subset = gene_mat2[,gene_category2$gene]
    
    row_anno_matrix2 <- as.matrix(gene_category2[colnames(gene_mat2_subset), "celltype"])
    
    unique_categories2 <- unique(gene_category2$celltype)
    row_anno_df2 <- as.data.frame(row_anno_matrix2)
    row_anno_df2$V1 <- factor(row_anno_df2$V1, levels = unique_categories2)
    if (nrow(gene_category_subset) != 1) {
      ht[[i]] <- Heatmap(gene_mat2_subset,
                         cluster_columns = T,
                         cluster_rows =  F,
                         show_column_dend = F,
                         row_names_side = "left",
                         name = "Gene Expression",
                         column_names_gp = gpar(fontsize = 27, rot = 45),
                         row_names_gp = gpar(fontsize = 27),
                         column_split = row_anno_df2,
                         column_gap = unit(1.5, "mm"),
                         col = colorRampPalette(c('blue', 'yellow'))(13),
                         column_names_rot = 45)
    } else {
      ht[[i]] <- Heatmap(gene_mat2_subset,
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
  if (major_category[j] == "T") {
    combined_heatmap <- ht[[1]] + ht[[2]] + ht[[3]] + ht[[4]]
    pdf(paste0("./result/Fig3A_", major_category[j], "_1.pdf"),
        width = pdf_width, height = pdf_height)
    draw(combined_heatmap)
    dev.off()
    
    combined_heatmap <- ht[[5]] + ht[[6]] + ht[[7]] + ht[[8]] + ht[[9]]
    pdf(paste0("./result/Fig3A_", major_category[j], "_2.pdf"),
        width = pdf_width, height = pdf_height)
    draw(combined_heatmap)
    dev.off()
  } else if (major_category[j] %in% c("B","Others")) {
    combined_heatmap <- Reduce(`+`, ht)
    
    pdf(paste0("./result/Fig3A_", major_category[j], ".pdf"),
        width = pdf_width, height = pdf_height)
    draw(combined_heatmap)
    dev.off()
  }
}

rm(list = ls())
