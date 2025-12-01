setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_and_treatment_separated/ABC_marker_analysis/")
# setwd("./analysis_with_cohorts_and_treatment_separated/ABC_marker_analysis/")
getwd()

head(gmt_sets)


head(deseq_res$AM14transfer_PL23)

pl23_genelist <- gsea_format(deseq_res$AM14transfer_PL23)
MRLlpr_genelist <- gsea_format(deseq_res$AM14MRLlpr)
r848_enelist <- gsea_format(deseq_res$AM14transfer_R848)
np_genelist <- gsea_format(deseq_res$B18transfer)



# GSEA of GMT gene sets --------------------------------------------------------
  pl23_ABC_gsea <- GSEA(pl23_genelist, TERM2GENE = gmt_sets, pvalueCutoff = 1)
  MRLlpr_ABC_gsea <- GSEA(MRLlpr_genelist, TERM2GENE = gmt_sets, pvalueCutoff = 1)
  r848_ABC_gsea <- GSEA(r848_enelist, TERM2GENE = gmt_sets, pvalueCutoff = 1)
  np_ABC_gsea <- GSEA(np_genelist, TERM2GENE = gmt_sets, pvalueCutoff = 1)
  
  head(pl23_ABC_gsea)
  nrow(pl23_ABC_gsea)
  pl23_ABC_gsea[, 1]
  pl23_ABC_gsea[, 5:7]
  pl23_ABC_gsea[8, ]
  
  gsea_dotplot(pl23_ABC_gsea, plottitle = "PL2-3 + 2DG vs PL2-3")
  ggsave(filename = "ABC markers - PL2-3.png", width = 6.5, height = 6, units = "in", dpi = 600)
  
  gsea_dotplot(MRLlpr_ABC_gsea, plottitle = "AM14 MRL/lpr: 2DG vs Control")
  ggsave(filename = "ABC markers - AM14 MRLlpr.png", width = 6.5, height = 6, units = "in", dpi = 600)
  
  gsea_dotplot(r848_ABC_gsea, plottitle = "R848 + 2DG vs R848")
  ggsave(filename = "ABC markers - R848.png", width = 6.5, height = 6, units = "in", dpi = 600)
  
  gsea_dotplot(np_ABC_gsea, plottitle = "NP + 2DG vs NP")
  ggsave(filename = "ABC markers - NP.png", width = 6.5, height = 4, units = "in", dpi = 600)


# GSEA_Histograms --------------------------------------------------------------
  gseaplot2(pl23_ABC_gsea, geneSetID = 1, title = pl23_ABC_gsea$Description[1])
  ggsave(filename = "GSEA_Histograms/Myc Targets v1 - PL2-3.png", width = 5, height = 4, units = "in", dpi = 600)
  
  gseaplot2(pl23_ABC_gsea, geneSetID = 2, title = pl23_ABC_gsea$Description[2])
  ggsave(filename = "GSEA_Histograms/Response to IFNG - PL2-3.png", width = 5, height = 4, units = "in", dpi = 600)
  
  gseaplot2(pl23_ABC_gsea, geneSetID = 3, title = pl23_ABC_gsea$Description[3])
  ggsave(filename = "GSEA_Histograms/GC vs Naive B Cell Up - PL2-3.png", width = 5, height = 4, units = "in", dpi = 600)
  
  gseaplot2(pl23_ABC_gsea, geneSetID = 4, title = pl23_ABC_gsea$Description[4])
  ggsave(filename = "GSEA_Histograms/Ag processing and presentation - PL2-3.png", width = 5, height = 4, units = "in", dpi = 600)
  
  gseaplot2(pl23_ABC_gsea, geneSetID = 5, title = pl23_ABC_gsea$Description[5])
  ggsave(filename = "GSEA_Histograms/GC B cell vs Plasma Cell Up - PL2-3.png", width = 5, height = 4, units = "in", dpi = 600)
  
  gseaplot2(MRLlpr_ABC_gsea, geneSetID = 1, title = stri_join("AM14 MRL/lpr: 2DG vs Control", MRLlpr_ABC_gsea$Description[1], sep = "\n"))
  ggsave(filename = "GSEA_Histograms/Response to IFNG - AM14 MRLlpr.png", width = 5, height = 4, units = "in", dpi = 600)
  
  
  
  p1 <- gseaplot2(pl23_ABC_gsea, geneSetID = 1, title = pl23_ABC_gsea$Description[1])
  p2 <- gseaplot2(pl23_ABC_gsea, geneSetID = 2, title = pl23_ABC_gsea$Description[2])
  p3 <- gseaplot2(pl23_ABC_gsea, geneSetID = 3, title = pl23_ABC_gsea$Description[3])
  p4 <- gseaplot2(pl23_ABC_gsea, geneSetID = 4, title = pl23_ABC_gsea$Description[4])
  plot_list(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])
  plot_list(p1, p3, p4)


