get_pathway_genes <- function(pw_file, dds, ID_conversion){
  pw_df <- read.csv(pw_file)
  colnames(pw_df) <- c("external_gene_name")
  
  norm_cts <- counts(dds, normalized = TRUE)
  
  col_names <- c("ensembl_gene_id", colnames(norm_cts))
  norm_cts <- data.frame(rownames(norm_cts), norm_cts)
  colnames(norm_cts) <- col_names
  norm_cts <- left_join(norm_cts, ID_conversion)
  norm_cts <- norm_cts[, 2:ncol(norm_cts)]
  
  res <- left_join(pw_df, norm_cts, by = "external_gene_name")
  rownames(res) <- res$external_gene_name
  # res <- res[, 2:ncol(res)]
  res
  # norm_cts
}



GO_TLR9 <- join_all(list(
  get_pathway_genes("pathway_gene_lists/GOBP_TOLL_LIKE_RECEPTOR_9_SIGNALING_PATHWAY.txt", dds_AM14trans, gene_IDs$AM14trans[, 1:2]),
  get_pathway_genes("pathway_gene_lists/GOBP_TOLL_LIKE_RECEPTOR_9_SIGNALING_PATHWAY.txt", dds_B18trans, gene_IDs$B18trans[, 1:2]),
  get_pathway_genes("pathway_gene_lists/GOBP_TOLL_LIKE_RECEPTOR_9_SIGNALING_PATHWAY.txt", dds_AM14MRLlpr, gene_IDs$AM14MRLlpr[, 1:2])),
  by = "external_gene_name"
)
rownames(GO_TLR9) <- GO_TLR9$external_gene_name
GO_TLR9 <- GO_TLR9[, c(2:ncol(GO_TLR9))]
GO_TLR9 <- as.matrix(GO_TLR9)
head(GO_TLR9)
ncol(GO_TLR9)
GO_TLR9 <- na.omit(GO_TLR9)
nrow(GO_TLR9)

annotation_df <- data.frame(coldata$Treatment)
rownames(annotation_df) <- coldata$Sample
colnames(annotation_df) <- c("Treatment")
annotation_df
# png(filename = stri_join(c("QC_results/Heatmaps/", filename_start,
#                            " - Regularized Log Transformation.png"),
#                          collapse = ""),
#     width = 1200, height = 1200, units = "px", pointsize = 10, res = 200,
#     bg = "white", family = "", symbolfamily="default")
# pheatmap(GO_TLR9, scale = "column", color=colorRampPalette(c("navy", "white", "red"))(50),
         # cluster_rows=FALSE, show_rownames=FALSE,
         # cluster_cols=TRUE, labels_col = colnames(GO),
         # annotation_col=annotation_df)
         
# dev.off()



pheatmap(GO_TLR9, annotation_col = annotation_df)
