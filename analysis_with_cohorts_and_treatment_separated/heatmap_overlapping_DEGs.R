setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_and_treatment_separated")
getwd()


pl23_DEG_df <- deseq_res$AM14transfer_PL23[, c(1, 2, 8, 11, 12)]
pl23_DEG_df <- filter(pl23_DEG_df, abs(log2FoldChange) > 0.6 & padj <= 0.05)
nrow(pl23_DEG_df)
head(pl23_DEG_df)

MRLlpr_DEG_df <- deseq_res$AM14MRLlpr[, c(1, 2, 8, 11, 12)]
MRLlpr_DEG_df <- filter(MRLlpr_DEG_df, abs(log2FoldChange) > 0.6 & padj <= 0.05)
nrow(MRLlpr_DEG_df)
head(MRLlpr_DEG_df)

merged_DEG_df <- pl23_DEG_df
colnames(merged_DEG_df)  <- c("ensembl_gene_id", "external_gene_name", "PL23_LFC", "PL23_pvalue", "PL23_padj")
merged_DEG_df <- inner_join(merged_DEG_df, MRLlpr_DEG_df)
colnames(merged_DEG_df)  <- c("ensembl_gene_id", "external_gene_name", "PL23_LFC", "PL23_pvalue", "PL23_padj", "MRLlpr_LFC", "MRLlpr_pvalue", "MRLlpr_padj")
head(merged_DEG_df)
nrow(merged_DEG_df)

ID_conversion <- merged_DEG_df[, 1:2]
ID_conversion

DEGs <- merged_DEG_df$ensembl_gene_id
head(DEGs)


zscore_matrix <- function(dds){
  zscores <- counts(dds, normalized=TRUE)
  res_colnames <- colnames(zscores)
  zscores <- t(zscores)
  zscores <- scale(zscores)
  zscores <- t(zscores)
  zscores <- data.frame(rownames(zscores), zscores)
  zscores <- zscores[rownames(zscores) %in% DEGs, ]
  colnames(zscores) <- c("ensembl_gene_id", res_colnames)
  rm(res_colnames)
  return(zscores)
}

pl23_zscores <- zscore_matrix(dds_AM14trans_PL23)
head(pl23_zscores)

MRLlpr_zscores <- zscore_matrix(dds_AM14MRLlpr)
head(MRLlpr_zscores)

DEG_zscores <- inner_join(pl23_zscores, MRLlpr_zscores)
DEG_zscores <- inner_join(ID_conversion, DEG_zscores)
rownames(DEG_zscores) <- DEG_zscores$external_gene_name
DEG_zscores <- as.matrix(DEG_zscores[, c(3:ncol(DEG_zscores))])
head(DEG_zscores)
DEG_zscores

overlapping_DEGs_list <- rownames(DEG_zscores)
overlapping_DEGs_list

coldata_df <- data.frame(Treatment = c(rep("Control", 3), rep("2DG", 3), rep("Control", 6), rep("2DG", 7)),
                         Experiment = c(rep("AM14 + PL2-3", 6), rep("AM14 MRLlpr", 13)))
coldata_df

png(filename = ("Output/DEG_heatmaps/DEGs common in AM14 PL2-3 and AM14 MRLlpr.png"),
    width = 4, height = 7, units = "in", pointsize = 6, res = 600,
    bg = "white", family = "", symbolfamily="default")

htmp <- ComplexHeatmap::pheatmap(DEG_zscores,
                                 cluster_rows = TRUE, show_rownames = TRUE,
                                 cluster_cols = TRUE, show_colnames = FALSE, 
                                 annotation_col = coldata_df, 
                                 scale = "none",
                                 cutree_rows = 3, 
                                 annotation_colors = list(Experiment = c("AM14 + PL2-3" = "#0f85a0", "AM14 MRLlpr" = "#00496f"),
                                                          Treatment = c(Control = "#edd746", "2DG" = "#dd4124")),
                                 # main = "DEGs Common to AM14 Adoptive Transfer PL2-3 +/- 2DG and AM14 MRL/lpr +/- 2DG",
                                 heatmap_legend_param = list(title = "Z-score", direction = "horizontal"))
                                 
draw(htmp, legend_grouping = "original", merge_legends = TRUE, heatmap_legend_side = "top")
dev.off()

