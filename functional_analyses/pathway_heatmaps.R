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
  res
}


annotation_df <- read.csv("raw_data/sample_info.csv")
sample_list <- annotation_df$Sample
annotation_df <- data.frame(annotation_df$Drug, annotation_df$heatmap_col1)
colnames(annotation_df) <- c("Treatment", "Experiment")
rownames(annotation_df) <- sample_list
annotation_df  

ann_colors4 <- list(Experiment = c("AM14 Transfer + PL2-3 Stimulation" = "#00496f", 
                                   "AM14 Transfer + R848 Stimulation" = "#0f85a0", 
                                   "B1-8 Transfer + NP Stimulation" = "#edd746", 
                                   "AM14 MRL/lpr" = "#dd4124"),
                    Treatment = c("2DG" = "darkgray", "Control" = "beige"))

GO_TLR9 <- join_all(list(
  get_pathway_genes("pathway_gene_lists/GOBP_TOLL_LIKE_RECEPTOR_9_SIGNALING_PATHWAY.txt", dds_AM14trans, gene_IDs$AM14trans[, 1:2]),
  get_pathway_genes("pathway_gene_lists/GOBP_TOLL_LIKE_RECEPTOR_9_SIGNALING_PATHWAY.txt", dds_B18trans, gene_IDs$B18trans[, 1:2]),
  get_pathway_genes("pathway_gene_lists/GOBP_TOLL_LIKE_RECEPTOR_9_SIGNALING_PATHWAY.txt", dds_AM14MRLlpr, gene_IDs$AM14MRLlpr[, 1:2])),
  by = "external_gene_name"
)
rownames(GO_TLR9) <- GO_TLR9$external_gene_name
GO_TLR9 <- GO_TLR9[, c(2:ncol(GO_TLR9))]
GO_TLR9 <- as.matrix(GO_TLR9)
GO_TLR9 <- na.omit(GO_TLR9)

png(filename = "Output/Functional_analyses/enrichGO_plots/TLR9 Signaling Pathway.png",
    width = 1200, height = 750, units = "px", pointsize = 10, res = 150,
    bg = "white", family = "", symbolfamily="default")
ComplexHeatmap::pheatmap(GO_TLR9, 
                         annotation_col = annotation_df, scale = "row", 
                         cluster_cols = FALSE, show_colnames = FALSE,
                         annotation_colors = ann_colors4,
                         main = "TLR9 Signaling Pathway",
                         heatmap_legend_param = list(title = "Z-score"))
dev.off()



