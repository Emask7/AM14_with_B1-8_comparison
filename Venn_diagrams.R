deg_list_for_venn <- function(res, lfc_cutoff){
  list(up = subset(res, res$padj <= 0.05 & res$log2FoldChange >= lfc_cutoff)[, 1], 
       down = subset(res, res$padj <= 0.05 & res$log2FoldChange <= (-1*lfc_cutoff))[, 1], 
       all = subset(res, res$padj <= 0.05 & abs(res$log2FoldChange) >= lfc_cutoff)[, 1])
}

DEG_lists <- list(
  PL23_2DG_vs_Ctrl = deg_list_for_venn(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, 0.5),
  R848_2DG_vs_Ctrl = deg_list_for_venn(deseq_res$AM14transfer$R848_2DG_v_Ctrl, 0.5),
  PL23_vs_R848 = deg_list_for_venn(deseq_res$AM14transfer$PL23_vs_R848, 0.5),
  PL23_vs_NP = deg_list_for_venn(deseq_res$PL23_v_NP, 0.5),
  B18 = deg_list_for_venn(deseq_res$B18transfer, 0.5),
  MRLlpr = deg_list_for_venn(deseq_res$AM14MRLlpr, 0.5)
)


# Comparing all groups ---------------------------------------------------------
  all_degs <- venndetail(list("AM14 Transfer: PL2-3 + 2DG vs PL2-3" = DEG_lists$PL23_2DG_vs_Ctrl$all,
                               "AM14 Transfer: R848 + 2DG vs R848" = DEG_lists$R848_2DG_vs_Ctrl$all,
                               "AM14 Transfer: PL2-3 vs R848" = DEG_lists$PL23_vs_R848$all,
                               "B1-8 Transfer: NP + 2DG vs NP" = DEG_lists$B18$all,
                               "AM14 MRL/lpr: 2DG vs Control" = DEG_lists$MRLlpr$all))
  png(filename = "Output/Venn_diagrams/All_DEGs.png", width = 2000, height = 1500, res = 250)
  plot(all_degs, type = "upset")
  dev.off()
  
  up_degs <- venndetail(list("AM14 Transfer: PL2-3 + 2DG vs PL2-3" = DEG_lists$PL23_2DG_vs_Ctrl$up,
                             "AM14 Transfer: R848 + 2DG vs R848" = DEG_lists$R848_2DG_vs_Ctrl$up,
                             "AM14 Transfer: PL2-3 vs R848" = DEG_lists$PL23_vs_R848$up,
                             "B1-8 Transfer: NP + 2DG vs NP" = DEG_lists$B18$up,
                             "AM14 MRL/lpr: 2DG vs Control" = DEG_lists$MRLlpr$up))
  png(filename = "Output/Venn_diagrams/Upregulated_DEGs.png", width = 2000, height = 1500, res = 250)
  plot(up_degs, type = "upset")
  dev.off()
  
  down_degs <- venndetail(list("AM14 Transfer: PL2-3 + 2DG vs PL2-3" = DEG_lists$PL23_2DG_vs_Ctrl$down,
                                  "AM14 Transfer: R848 + 2DG vs R848" = DEG_lists$R848_2DG_vs_Ctrl$down,
                                  "AM14 Transfer: PL2-3 vs R848" = DEG_lists$PL23_vs_R848$down,
                                  "B1-8 Transfer: NP + 2DG vs NP" = DEG_lists$B18$down,
                                  "AM14 MRL/lpr: 2DG vs Control" = DEG_lists$MRLlpr$down))
  png(filename = "Output/Venn_diagrams/Downregulated_DEGs.png", width = 2000, height = 1500, res = 250)
  plot(down_degs, type = "upset")
  dev.off()
  
  
  PL23_down_vs_others_up <- venndetail(list("AM14 Transfer: PL2-3 + 2DG vs PL2-3" = DEG_lists$PL23_2DG_vs_Ctrl$down,
                                             "AM14 Transfer: R848 + 2DG vs R848" = DEG_lists$R848_2DG_vs_Ctrl$up,
                                             "AM14 Transfer: PL2-3 vs R848" = DEG_lists$PL23_vs_R848$up,
                                             "B1-8 Transfer: NP + 2DG vs NP" = DEG_lists$B18$up,
                                             "AM14 MRL/lpr: 2DG vs Control" = DEG_lists$MRLlpr$up))
  png(filename = "Output/Venn_diagrams/PL23_down_vs_others_up.png", width = 2000, height = 1500, res = 250)
  plot(PL23_down_vs_others_up, type = "upset")
  dev.off()
  
  PL23_up_vs_others_all <- venndetail(list("AM14 Transfer: PL2-3 + 2DG vs PL2-3" = DEG_lists$PL23_2DG_vs_Ctrl$up,
                                            "AM14 Transfer: R848 + 2DG vs R848" = DEG_lists$R848_2DG_vs_Ctrl$all,
                                            "AM14 Transfer: PL2-3 vs R848" = DEG_lists$PL23_vs_R848$all,
                                            "B1-8 Transfer: NP + 2DG vs NP" = DEG_lists$B18$all,
                                            "AM14 MRL/lpr: 2DG vs Control" = DEG_lists$MRLlpr$all))
  png(filename = "Output/Venn_diagrams/PL23_up_vs_others_all.png", width = 2000, height = 1500, res = 250)
  plot(PL23_up_vs_others_all, type = "upset")
  dev.off()
  
  
  
  
# Comparing effects of 2DG on different conditions -----------------------------
  up_2DGvCtrl_venn <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG\nvs PL2-3\n" = DEG_lists$PL23_2DG_vs_Ctrl$up,
                                      "AM14 Transfer:\nR848 + 2DG\nvs R848\n" = DEG_lists$R848_2DG_vs_Ctrl$up,
                                      "AM14 MRL/lpr:\n2DG vs Control\n" = DEG_lists$MRLlpr$up,
                                      "B1-8 Transfer:\nNP+2DG vs NP\n" = DEG_lists$B18$up))
  dev.off()
  plot(up_2DGvCtrl_venn, mycol = met.brewer("Johnson", n = 8, direction = 1),
       filename = "Output/Venn_diagrams/2DG_vs_Control-Upregulated.png",
       margin = 0.05, cat.cex = 1, cex = 1.5)
  # detail(up_2DGvCtrl_venn)
  
  down_2DGvCtrl_venn <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG\nvs PL2-3\n" = DEG_lists$PL23_2DG_vs_Ctrl$down,
                                      "AM14 Transfer:\nR848 + 2DG\nvs R848\n" = DEG_lists$R848_2DG_vs_Ctrl$down,
                                      "AM14 MRL/lpr:\n2DG vs Control\n" = DEG_lists$MRLlpr$down,
                                      "B1-8 Transfer:\nNP+2DG vs NP\n" = DEG_lists$B18$down))
  dev.off()
  plot(down_2DGvCtrl_venn, mycol = met.brewer("Johnson", n = 8, direction = -1),
       filename = "Output/Venn_diagrams/2DG_vs_Control-Downregulated.png",
       margin = 0.05, cat.cex = 1, cex = 1.5)
  # detail(down_2DGvCtrl_venn)

# # Comparing only gene sets with overlaps ---------------------------------------
#   up_venn <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG\nvs PL2-3\n" = DEG_lists$PL23_2DG_vs_Ctrl$up,
#                              "AM14 Transfer:\nR848 + 2DG\nvs R848\n" = DEG_lists$R848_2DG_vs_Ctrl$up,
#                              "AM14 Transfer:\nPL2-3 vs R848\n" = DEG_lists$PL23_vs_R848$up,
#                              "AM14/B1-8 Transfers:\nPL2-3 vs NP" = DEG_lists$PL23_vs_NP$up,
#                              "AM14 MRL/lpr:\n2DG vs Control\n" = DEG_lists$MRLlpr$up))
#   dev.off()
#   plot(up_venn, mycol = met.brewer("Johnson", n = 10, direction = 1),
#        # filename = "Output/Venn_diagrams/2DG_vs_Control-Upregulated.png",
#        margin = 0.05, cat.cex = 1, cex = 1.5)
#   # detail(up_venn)
#   
#   down_venn <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG\nvs PL2-3\n" = DEG_lists$PL23_2DG_vs_Ctrl$down,
#                              "AM14 Transfer:\nR848 + 2DG\nvs R848\n" = DEG_lists$R848_2DG_vs_Ctrl$down,
#                              "AM14 Transfer:\nPL2-3 vs R848\n" = DEG_lists$PL23_vs_R848$down,
#                              "AM14/B1-8 Transfers:\nPL2-3 vs NP" = DEG_lists$PL23_vs_NP$down,
#                              "AM14 MRL/lpr:\n2DG vs Control\n" = DEG_lists$MRLlpr$down))
#   dev.off()
#   plot(down_venn, mycol = met.brewer("Johnson", n = 10, direction = -1),
#        # filename = "Output/Venn_diagrams/2DG_vs_Control-downregulated.png",
#        margin = 0.05, cat.cex = 1, cex = 1.5)
#   # detail(down_venn)
# 
# # Comparing autoAg specific response to others ---------------------------------
#   up_PL23_v_others_venn <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG\nvs PL2-3\n" = DEG_lists$PL23_2DG_vs_Ctrl$up,
#                                            "AM14 Transfer:\nPL2-3\nvs R848\n" = DEG_lists$PL23_vs_R848$up,
#                                            "AM14/B1-8 Transfers:\nPL2-3 vs NP\n" = DEG_lists$PL23_vs_NP$up))
#   dev.off()
#   plot(up_PL23_v_others_venn, mycol = met.brewer("Johnson", n = 8, direction = 1),
#        # filename = "Output/venn_diagram_2DG_vs_Ctrl_up.png",
#        margin = 0.05, cat.cex = 1, cex = 1.5)
#   # detail(up_PL23_v_others_venn)
#   
#   down_PL23_v_others_venn <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG\nvs PL2-3\n" = DEG_lists$PL23_2DG_vs_Ctrl$down,
#                                            "AM14 Transfer:\nPL2-3\nvs R848\n" = DEG_lists$PL23_vs_R848$down,
#                                            "AM14/B1-8 Transfers:\nPL2-3 vs NP\n" = DEG_lists$PL23_vs_NP$down))
#   dev.off()
#   plot(down_PL23_v_others_venn, mycol = met.brewer("Johnson", n = 8, direction = -1),
#        # filename = "Output/venn_diagram_2DG_vs_Ctrl_down.png",
#        margin = 0.05, cat.cex = 1, cex = 1.5)
#   detail(down_PL23_v_others_venn)
# 
