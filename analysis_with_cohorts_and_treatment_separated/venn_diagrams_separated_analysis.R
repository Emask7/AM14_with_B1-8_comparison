setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_and_treatment_separated")
getwd()

DEG_lists <- list(
  PL23_2DG_vs_Ctrl = deg_list_for_venn(deseq_res$AM14transfer_PL23, 1, 0.05),
  R848_2DG_vs_Ctrl = deg_list_for_venn(deseq_res$AM14transfer_R848, 1, 0.05),
  B18 = deg_list_for_venn(deseq_res$B18transfer, 1, 0.05),
  MRLlpr = deg_list_for_venn(deseq_res$AM14MRLlpr, 1, 0.05)
)

DEG_lists_lfc0.6 <- list(
  PL23_2DG_vs_Ctrl = deg_list_for_venn(deseq_res$AM14transfer_PL23, 0.6, 0.05),
  R848_2DG_vs_Ctrl = deg_list_for_venn(deseq_res$AM14transfer_R848, 0.6, 0.05),
  B18 = deg_list_for_venn(deseq_res$B18transfer, 0.6, 0.05),
  MRLlpr = deg_list_for_venn(deseq_res$AM14MRLlpr, 0.6, 0.05)
)


# Comparing all groups ---------------------------------------------------------
  # all_degs <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG vs PL2-3" = DEG_lists$PL23_2DG_vs_Ctrl$all,
  #                             "AM14 Transfer:\nR848 + 2DG vs R848" = DEG_lists$R848_2DG_vs_Ctrl$all,
  #                             "AM14 MRL/lpr:\n2DG vs Control" = DEG_lists$MRLlpr$all))
  # png(filename = "Output/Venn_diagrams/All_DEGs.png", width = 2000, height = 1500, res = 250)
  # plot(all_degs, margin = 0.1)
  # dev.off()
  # 
  # up_degs <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG vs PL2-3" = DEG_lists$PL23_2DG_vs_Ctrl$up,
  #                             "AM14 Transfer:\nR848 + 2DG vs R848" = DEG_lists$R848_2DG_vs_Ctrl$up,
  #                             "AM14 MRL/lpr:\n2DG vs Control" = DEG_lists$MRLlpr$up))
  # png(filename = "Output/Venn_diagrams/up_DEGs.png", width = 2000, height = 1500, res = 250)
  # plot(up_degs, margin = 0.1)
  # dev.off()
  # 
  # down_degs <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG vs PL2-3" = DEG_lists$PL23_2DG_vs_Ctrl$down,
  #                            "AM14 Transfer:\nR848 + 2DG vs R848" = DEG_lists$R848_2DG_vs_Ctrl$down,
  #                            "AM14 MRL/lpr:\n2DG vs Control" = DEG_lists$MRLlpr$down))
  # png(filename = "Output/Venn_diagrams/down_DEGs.png", width = 2000, height = 1500, res = 250)
  # plot(down_degs, margin = 0.1)
  # dev.off()
  
  
  # DEGs with LFC > +/- 1 ------------------------------------------------------
    all_degs <- venndetail(list("AM14 Transfer: PL2-3 + 2DG vs PL2-3" = DEG_lists$PL23_2DG_vs_Ctrl$all,
                                 "AM14 Transfer: R848 + 2DG vs R848" = DEG_lists$R848_2DG_vs_Ctrl$all,
                                 "B1-8 Transfer: NP + 2DG vs NP" = DEG_lists$B18$all,
                                 "AM14 MRL/lpr: 2DG vs Control" = DEG_lists$MRLlpr$all))
    png(filename = "Output/Venn_diagrams/All_DEGs.png", width = 2000, height = 1500, res = 250)
    plot(all_degs, type = "upset")
    dev.off()
    
    up_degs <- venndetail(list("AM14 Transfer: PL2-3 + 2DG vs PL2-3" = DEG_lists$PL23_2DG_vs_Ctrl$up,
                               "AM14 Transfer: R848 + 2DG vs R848" = DEG_lists$R848_2DG_vs_Ctrl$up,
                               "B1-8 Transfer: NP + 2DG vs NP" = DEG_lists$B18$up,
                               "AM14 MRL/lpr: 2DG vs Control" = DEG_lists$MRLlpr$up))
    png(filename = "Output/Venn_diagrams/Upregulated_DEGs.png", width = 2000, height = 1500, res = 250)
    plot(up_degs, type = "upset")
    dev.off()
    
    down_degs <- venndetail(list("AM14 Transfer: PL2-3 + 2DG vs PL2-3" = DEG_lists$PL23_2DG_vs_Ctrl$down,
                                 "AM14 Transfer: R848 + 2DG vs R848" = DEG_lists$R848_2DG_vs_Ctrl$down,
                                 "B1-8 Transfer: NP + 2DG vs NP" = DEG_lists$B18$down,
                                 "AM14 MRL/lpr: 2DG vs Control" = DEG_lists$MRLlpr$down))
    png(filename = "Output/Venn_diagrams/Downregulated_DEGs.png", width = 2000, height = 1500, res = 250)
    plot(down_degs, type = "upset")
    dev.off()
    
    
    
    all_degs <- venndetail(list("PL2-3" = DEG_lists$PL23_2DG_vs_Ctrl$all,
                                "R848" = DEG_lists$R848_2DG_vs_Ctrl$all,
                                "MRL/lpr" = DEG_lists$MRLlpr$all,
                                "B1-8" = DEG_lists$B18$all))
    plot(all_degs, margin = 0.07, cat.cex = 0.5, cex = 1.7,
         mycol = c("dodgerblue", "darkorange1", "seagreen3", "orchid3"),
         filename = "Output/Venn_diagrams/All_DEGs_venn.png")
    
    up_degs <- venndetail(list("PL2-3" = DEG_lists$PL23_2DG_vs_Ctrl$up,
                                "R848" = DEG_lists$R848_2DG_vs_Ctrl$up,
                                "MRL/lpr" = DEG_lists$MRLlpr$up,
                                "B1-8" = DEG_lists$B18$up))
    plot(up_degs, margin = 0.07, cat.cex = 0.5, cex = 1.7,
         mycol = c("dodgerblue", "darkorange1", "seagreen3", "orchid3"),
         filename = "Output/Venn_diagrams/Upregulated_DEGs_venn.png")
    
    
    down_degs <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG\nvs PL2-3\n" = DEG_lists$PL23_2DG_vs_Ctrl$down,
                               "AM14 Transfer:\nR848 + 2DG\nvs R848\n" = DEG_lists$R848_2DG_vs_Ctrl$down,
                               "AM14 MRL/lpr:\n2DG vs Control" = DEG_lists$MRLlpr$down,
                               "B1-8 Transfer:\nNP + 2DG vs NP" = DEG_lists$B18$down))
    down_degs <- venndetail(list("PL2-3" = DEG_lists$PL23_2DG_vs_Ctrl$down,
                                 "R848" = DEG_lists$R848_2DG_vs_Ctrl$down,
                                 "MRL/lpr" = DEG_lists$MRLlpr$down))
                                 # "B1-8 Transfer:\nNP + 2DG vs NP" = DEG_lists$B18$down))
    
    plot(down_degs, margin = 0.07, cat.cex = 0.5, cex = 1.7,
         mycol = c("dodgerblue", "darkorange1", "seagreen3"),
         filename = "Output/Venn_diagrams/Downregulated_DEGs_venn.png")
    
    
    
    
  
  
  # DEGs with LFC > +/- 0.6 ----------------------------------------------------
    all_degs_0.6 <- venndetail(list("AM14 Transfer: PL2-3 + 2DG vs PL2-3" = DEG_lists_lfc0.6$PL23_2DG_vs_Ctrl$all,
                                "AM14 Transfer: R848 + 2DG vs R848" = DEG_lists_lfc0.6$R848_2DG_vs_Ctrl$all,
                                "B1-8 Transfer: NP + 2DG vs NP" = DEG_lists_lfc0.6$B18$all,
                                "AM14 MRL/lpr: 2DG vs Control" = DEG_lists_lfc0.6$MRLlpr$all))
    png(filename = "Output/Venn_diagrams/All_DEGs_LFC0.6.png", width = 2000, height = 1500, res = 250)
    plot(all_degs_0.6, type = "upset")
    dev.off()
    
    up_degs_0.6 <- venndetail(list("AM14 Transfer: PL2-3 + 2DG vs PL2-3" = DEG_lists_lfc0.6$PL23_2DG_vs_Ctrl$up,
                               "AM14 Transfer: R848 + 2DG vs R848" = DEG_lists_lfc0.6$R848_2DG_vs_Ctrl$up,
                               "B1-8 Transfer: NP + 2DG vs NP" = DEG_lists_lfc0.6$B18$up,
                               "AM14 MRL/lpr: 2DG vs Control" = DEG_lists_lfc0.6$MRLlpr$up))
    png(filename = "Output/Venn_diagrams/Upregulated_DEGs_LFC0.6.png", width = 2000, height = 1500, res = 250)
    plot(up_degs_0.6, type = "upset")
    dev.off()
    
    down_degs_0.6 <- venndetail(list("AM14 Transfer: PL2-3 + 2DG vs PL2-3" = DEG_lists_lfc0.6$PL23_2DG_vs_Ctrl$down,
                                 "AM14 Transfer: R848 + 2DG vs R848" = DEG_lists_lfc0.6$R848_2DG_vs_Ctrl$down,
                                 "B1-8 Transfer: NP + 2DG vs NP" = DEG_lists_lfc0.6$B18$down,
                                 "AM14 MRL/lpr: 2DG vs Control" = DEG_lists_lfc0.6$MRLlpr$down))
    png(filename = "Output/Venn_diagrams/Downregulated_DEGs_LFC0.6.png", width = 2000, height = 1500, res = 250)
    plot(down_degs_0.6, type = "upset")
    dev.off()
  
  
    
    all_degs_0.6 <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG\nvs PL2-3\n" = DEG_lists_lfc0.6$PL23_2DG_vs_Ctrl$all,
                                    "AM14 Transfer:\nR848 + 2DG\nvs R848\n" = DEG_lists_lfc0.6$R848_2DG_vs_Ctrl$all,
                                    "B1-8 Transfer:\nNP + 2DG vs NP" = DEG_lists_lfc0.6$B18$all,
                                    "AM14 MRL/lpr:\n2DG vs Control" = DEG_lists_lfc0.6$MRLlpr$all))
    plot(all_degs_0.6, margin = 0.07, cat.cex = 1.1, cex = 1.2,
         filename = "Output/Venn_diagrams/All_DEGs_LFC0.6_venn.png")

    up_degs_0.6 <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG\nvs PL2-3\n" = DEG_lists_lfc0.6$PL23_2DG_vs_Ctrl$up,
                                   "AM14 Transfer:\nR848 + 2DG\nvs R848\n" = DEG_lists_lfc0.6$R848_2DG_vs_Ctrl$up,
                                   "B1-8 Transfer:\nNP + 2DG vs NP" = DEG_lists_lfc0.6$B18$up,
                                   "AM14 MRL/lpr:\n2DG vs Control" = DEG_lists_lfc0.6$MRLlpr$up))
    plot(up_degs_0.6, margin = 0.07, cat.cex = 1.1, cex = 1.2,
         filename = "Output/Venn_diagrams/Upregulated_DEGs_LFC0.6_venn.png")
    
    
    down_degs_0.6 <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG\nvs PL2-3\n" = DEG_lists_lfc0.6$PL23_2DG_vs_Ctrl$down,
                                     "AM14 Transfer:\nR848 + 2DG\nvs R848\n" = DEG_lists_lfc0.6$R848_2DG_vs_Ctrl$down,
                                     "B1-8 Transfer:\nNP + 2DG vs NP" = DEG_lists_lfc0.6$B18$down,
                                     "AM14 MRL/lpr:\n2DG vs Control" = DEG_lists_lfc0.6$MRLlpr$down))
    plot(down_degs_0.6, margin = 0.07, cat.cex = 1.1, cex = 1.2,
         filename = "Output/Venn_diagrams/Downregulated_DEGs_LFC0.6_venn.png")
    