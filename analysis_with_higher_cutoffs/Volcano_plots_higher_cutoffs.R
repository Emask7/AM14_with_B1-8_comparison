vp_data = list(
  AM14transfer = list(
    PL23_2DG_v_Ctrl = deseq_res$AM14transfer$PL23_2DG_v_Ctrl[, c(2, 8, 12)],
    R848_2DG_v_Ctrl = deseq_res$AM14transfer$R848_2DG_v_Ctrl[, c(2, 8, 12)],
    PL23_vs_R848 = deseq_res$AM14transfer$PL23_vs_R848[, c(2, 8, 12)]),
  B18transfer = deseq_res$B18transfer[, c(2, 8, 12)],
  AM14MRLlpr = deseq_res$AM14MRLlpr[, c(2, 8, 12)])
head(vp_data$AM14MRLlpr)


volcano_wrapper(vp_data$AM14transfer$PL23_2DG_v_Ctrl, "AM14 Transfer: PL2-3 + 2DG vs PL2-3",
                "Upregulated in PL2-3                                                        Upregulated in PL2-3 + 2DG",
                1, 0.05, "AM14 Transfer - PL2-3_2DG_vs_Control - padj 0.05")

volcano_wrapper(vp_data$AM14transfer$R848_2DG_v_Ctrl, "AM14 Transfer: R848 + 2DG vs R848",
                "Upregulated in R848                                                        Upregulated in R848 + 2DG",
                1, 0.05, "AM14 Transfer - R848_2DG_vs_Control - padj 0.05")

volcano_wrapper(vp_data$B18transfer, "B1-8 Transfer: NP + 2DG vs NP", 
                "Upregulated in NP                                                            Upregulated in NP + 2DG",
                1, 0.05, "B1-8 Transfer - padj 0.05")

volcano_wrapper(vp_data$AM14MRLlpr, "AM14 MRL/lpr: 2DG vs Control", 
                "Upregulated in Control                                                               Upregulated in 2DG",
                1, 0.05, "AM14 MRLlpr - padj 0.05")



volcano_wrapper(vp_data$AM14transfer$PL23_2DG_v_Ctrl, "AM14 Transfer: PL2-3 + 2DG vs PL2-3",
                "Upregulated in PL2-3                                                        Upregulated in PL2-3 + 2DG",
                1, 0.01, "AM14 Transfer - PL2-3_2DG_vs_Control - padj 0.01")

volcano_wrapper(vp_data$AM14transfer$R848_2DG_v_Ctrl, "AM14 Transfer: R848 + 2DG vs R848",
                "Upregulated in R848                                                        Upregulated in R848 + 2DG",
                1, 0.01, "AM14 Transfer - R848_2DG_vs_Control - padj 0.01")

volcano_wrapper(vp_data$B18transfer, "B1-8 Transfer: NP + 2DG vs NP", 
                "Upregulated in NP                                                            Upregulated in NP + 2DG",
                1, 0.01, "B1-8 Transfer - padj 0.01")

volcano_wrapper(vp_data$AM14MRLlpr, "AM14 MRL/lpr: 2DG vs Control", 
                "Upregulated in Control                                                               Upregulated in 2DG",
                1, 0.01, "AM14 MRLlpr - padj 0.01")

