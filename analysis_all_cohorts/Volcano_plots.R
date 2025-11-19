# Run "data_visualization_functions.R" script first 

vp_data = list(
  AM14transfer = list(
    PL23_2DG_v_Ctrl = deseq_res$AM14transfer$PL23_2DG_v_Ctrl[, c(2, 6, 10)],
    R848_2DG_v_Ctrl = deseq_res$AM14transfer$R848_2DG_v_Ctrl[, c(2, 6, 10)],
    PL23_vs_R848 = deseq_res$AM14transfer$PL23_vs_R848[, c(2, 6, 10)]),
  B18transfer = deseq_res$B18transfer[, c(2, 6, 10)],
  PL23vNP = deseq_res$PL23_v_NP[, c(2, 6, 10)],
  AM14MRLlpr = deseq_res$AM14MRLlpr[, c(2, 6, 10)])
head(vp_data$AM14MRLlpr)

volcano_wrapper(vp_data$AM14transfer$PL23_2DG_v_Ctrl, 1, "AM14 Transfer: PL2-3 + 2DG vs PL2-3",
                "Upregulated in PL2-3                                                                                  Upregulated in PL2-3 + 2DG",
                c(-10, 5), c(0, 8), "AM14_Transfer-PL2-3_2DG_vs_Control")

volcano_wrapper(vp_data$AM14transfer$R848_2DG_v_Ctrl, 1, "AM14 Transfer: R848 + 2DG vs R848", 
                "Upregulated in R848                                                                                  Upregulated in R848 + 2DG",
                c(-7.5, 5), c(0, 5), "AM14_Transfer-R848_2DG_vs_Control")

volcano_wrapper(vp_data$AM14transfer$PL23_vs_R848, 1, "AM14 Transfer: PL2-3 vs R848",
                "Upregulated in R848                                                                                            Upregulated in PL2-3",
                c(-12, 5), c(0, 6), "AM14_Transfer-PL2-3_vs_R848")

volcano_wrapper(vp_data$B18transfer, 1, "B1-8 Transfer: NP + 2DG vs NP", 
                "Upregulated in NP                                                                                            Upregulated in NP + 2DG",
                c(-5, 3), c(0, 8), "B1-8_Transfer-NP_2DG_vs_NP")

volcano_wrapper(vp_data$AM14MRLlpr, 1, "AM14 MRL/lpr: 2DG vs Control",
                "Upregulated in Control                                                                                         Upregulated in 2DG",
                c(-5, 2.5), c(0, 15), "AM14_MRLlpr")

volcano_wrapper(vp_data$PL23vNP, 1, "AM14 Transfer PL2-3 vs B1-8 Transfer NP",
                "Upregulated in NP                                                                                            Upregulated in PL2-3",
                c(-6.5, 12), c(0, 182), "PL23_vs_NP")





volcano_wrapper(vp_data$AM14transfer$PL23_2DG_v_Ctrl, 0.5, "AM14 Transfer: PL2-3 + 2DG vs PL2-3",
                "Upregulated in PL2-3                                   Upregulated in PL2-3 + 2DG",
                c(-10, 5), c(0, 8), "AM14_Transfer-PL2-3_2DG_vs_Control - 0.5 LFC")

volcano_wrapper(vp_data$AM14transfer$R848_2DG_v_Ctrl, 0.5, "AM14 Transfer: R848 + 2DG vs R848", 
                "Upregulated in R848                                   Upregulated in R848 + 2DG",
                c(-7.5, 5), c(0, 5), "AM14_Transfer-R848_2DG_vs_Control - 0.5 LFC")

volcano_wrapper(vp_data$AM14transfer$PL23_vs_R848, 0.5, "AM14 Transfer: PL2-3 vs R848",
                "Upregulated in R848                                        Upregulated in PL2-3",
                c(-12, 5), c(0, 6), "AM14_Transfer-PL2-3_vs_R848 - 0.5 LFC")

volcano_wrapper(vp_data$B18transfer, 0.5, "B1-8 Transfer: NP + 2DG vs NP",  
                "Upregulated in NP                                   Upregulated in NP + 2DG",
                c(-5, 3), c(0, 17.5), "B1-8_Transfer-NP_2DG_vs_NP - 0.5 LFC")

volcano_wrapper(vp_data$AM14MRLlpr, 0.5, "AM14 MRL/lpr: 2DG vs Control", 
                "Upregulated in Control                                   Upregulated in 2DG",
                c(-5, 2.5), c(0, 15), "AM14_MRLlpr - 0.5 LFC")

volcano_wrapper(vp_data$PL23vNP, 0.5, "AM14 Transfer PL2-3 vs B1-8 Transfer NP",
                "Upregulated in NP                                   Upregulated in PL2-3",
                c(-6.5, 12), c(0, 182), "PL23_vs_NP - 0.5 LFC")
