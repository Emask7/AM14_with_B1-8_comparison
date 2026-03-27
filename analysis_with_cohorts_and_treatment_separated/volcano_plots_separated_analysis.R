setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_and_treatment_separated")
getwd()

vp_data = list(
  PL23_2DG_v_Ctrl = deseq_res$AM14transfer_PL23[, c(2, 8, 12)],
  R848_2DG_v_Ctrl = deseq_res$AM14transfer_R848[, c(2, 8, 12)],
  AM14MRLlpr = deseq_res$AM14MRLlpr[, c(2, 8, 12)],
  B18 = deseq_res$B18transfer[, c(2, 8, 12)])
head(vp_data$AM14MRLlpr)


# confirm_foldchange(dds_AM14trans, "Crisp1", gene_IDs$AM14trans)
# confirm_foldchange(dds_B18trans, "Dmrt2", gene_IDs$B18trans)
# confirm_foldchange(dds_PL23vNP, "Lars2", gene_IDs$AM14trans)

volcano_wrapper(vp_data$PL23_2DG_v_Ctrl, "AM14 Transfer: PL2-3 + 2DG vs PL2-3",
                "Upregulated in PL2-3                                               Upregulated in PL2-3 + 2DG",
                0.6, 0.05, "AM14_Transfer-PL2-3_2DG_vs_Control")

volcano_wrapper(vp_data$R848_2DG_v_Ctrl, "AM14 Transfer: R848 + 2DG vs R848", 
                "Upregulated in R848                                                Upregulated in R848 + 2DG",
                0.6, 0.05, "AM14_Transfer-R848_2DG_vs_Control")

volcano_wrapper(vp_data$AM14MRLlpr, "AM14 MRL/lpr: 2DG vs Control",
                "Upregulated in Control                                                              Upregulated in 2DG",
                0.6, 0.05, "AM14_MRLlpr")

volcano_wrapper(vp_data$B18, "B1-8 Transfer: NP + 2DG vs NP",
                "Upregulated in NP                                                         Upregulated in NP + 2DG",
                0.6, 0.05, "B18_Transfer")
