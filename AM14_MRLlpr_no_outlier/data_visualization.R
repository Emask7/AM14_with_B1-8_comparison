setwd("./AM14_MRLlpr_no_outlier/")
getwd()

# DEG venn diagram -------------------------------------------------------------
  DEG_lists_noout <- list(
    PL23_2DG_vs_Ctrl = deg_list_for_venn(deseq_res$AM14transfer_PL23, 0.6, 0.05),
    R848_2DG_vs_Ctrl = deg_list_for_venn(deseq_res$AM14transfer_R848, 0.6, 0.05),
    B18 = deg_list_for_venn(deseq_res$B18transfer, 0.6, 0.05),
    MRLlpr = deg_list_for_venn(deseq_res_MRLlpr_noout, 0.6, 0.05)
  )
  
  all_degs <- venndetail(list("AM14 Transfer: PL2-3 + 2DG vs PL2-3" = DEG_lists_noout$PL23_2DG_vs_Ctrl$all,
                              # "AM14 Transfer: R848 + 2DG vs R848" = DEG_lists_noout$R848_2DG_vs_Ctrl$all,
                              # "B1-8 Transfer: NP + 2DG vs NP" = DEG_lists_noout$B18$all,
                              "AM14 MRL/lpr: 2DG vs Control" = DEG_lists_noout$MRLlpr$all))
  png(filename = "Output/Venn_diagrams/All_DEGs.png", width = 2000, height = 1500, res = 250)
  plot(all_degs, type = "upset")
  dev.off()

# DEG overlap heatmap ----------------------------------------------------------
  # Run "data_visualization_functions.R" script first 
  
  pl23_DEG_df <- deseq_res$AM14transfer_PL23[, c(1, 2, 8, 11, 12)]
  pl23_DEG_df <- filter(pl23_DEG_df, abs(log2FoldChange) > 0.6 & padj <= 0.05)
  nrow(pl23_DEG_df)
  head(pl23_DEG_df)
  
  MRLlpr_DEG_df_noout <- deseq_res_MRLlpr_noout[, c(1, 2, 8, 11, 12)]
  MRLlpr_DEG_df_noout <- filter(MRLlpr_DEG_df_noout, abs(log2FoldChange) > 0.6 & padj <= 0.05)
  nrow(MRLlpr_DEG_df_noout)
  head(MRLlpr_DEG_df_noout)
  
  merged_DEG_df_noout <- pl23_DEG_df
  colnames(merged_DEG_df_noout)  <- c("ensembl_gene_id", "external_gene_name", "PL23_LFC", "PL23_pvalue", "PL23_padj")
  merged_DEG_df_noout <- inner_join(merged_DEG_df_noout, MRLlpr_DEG_df_noout)
  colnames(merged_DEG_df_noout)  <- c("ensembl_gene_id", "external_gene_name", "PL23_LFC", "PL23_pvalue", "PL23_padj", "MRLlpr_LFC", "MRLlpr_pvalue", "MRLlpr_padj")
  head(merged_DEG_df_noout)
  nrow(merged_DEG_df_noout)
  
  ID_conversion <- merged_DEG_df_noout[, 1:2]
  ID_conversion
  
  DEGs <- merged_DEG_df_noout$ensembl_gene_id
  head(DEGs)
  
  pl23_zscores <- zscore_matrix(dds_AM14trans_PL23)
  head(pl23_zscores)
  
  MRLlpr_zscores_noout <- zscore_matrix(dds_AM14MRLlpr_noout)
  head(MRLlpr_zscores_noout)
  
  DEG_zscores_noout <- inner_join(pl23_zscores, MRLlpr_zscores_noout)
  DEG_zscores_noout <- inner_join(ID_conversion, DEG_zscores_noout)
  rownames(DEG_zscores_noout) <- DEG_zscores_noout$external_gene_name
  DEG_zscores_noout <- as.matrix(DEG_zscores_noout[, c(3:ncol(DEG_zscores_noout))])
  head(DEG_zscores_noout)
  DEG_zscores_noout
  
  coldata_df_noout <- data.frame(Treatment = c(rep("Control", 3), rep("2DG", 3), rep("Control", 6), rep("2DG", 6)),
                                 Experiment = c(rep("AM14 + PL2-3", 6), rep("AM14 MRLlpr", 12)))
  coldata_df_noout
  
  png(filename = ("Output/DEG_heatmaps/DEGs common in AM14 PL2-3 and AM14 MRLlpr - no outlier.png"),
      width = 2000, height = 3500, units = "px", pointsize = 8, res = 400,
      bg = "white", family = "", symbolfamily="default")
  
  htmp <- ComplexHeatmap::pheatmap(DEG_zscores_noout,
                                   cluster_rows = TRUE, show_rownames = TRUE,
                                   cluster_cols = TRUE, show_colnames = FALSE, 
                                   annotation_col = coldata_df_noout, 
                                   scale = "none",
                                   cutree_rows = 3, 
                                   annotation_colors = list(Experiment = c("AM14 + PL2-3" = "#0f85a0", "AM14 MRLlpr" = "#00496f"),
                                                            Treatment = c(Control = "#edd746", "2DG" = "#dd4124")),
                                   # main = "DEGs Common to AM14 Adoptive Transfer PL2-3 +/- 2DG and AM14 MRL/lpr +/- 2DG",
                                   heatmap_legend_param = list(title = "Z-score", direction = "horizontal"))
  
  draw(htmp, legend_grouping = "original", merge_legends = TRUE, heatmap_legend_side = "top")
  dev.off()
  
  