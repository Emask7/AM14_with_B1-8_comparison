DEG_lists_p0.05_p0.05 <- list(
  PL23_2DG_vs_Ctrl = deg_list_for_venn(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, 1, 0.05),
  R848_2DG_vs_Ctrl = deg_list_for_venn(deseq_res$AM14transfer$R848_2DG_v_Ctrl, 1, 0.05),
  PL23_vs_R848 = deg_list_for_venn(deseq_res$AM14transfer$PL23_vs_R848, 1, 0.05),
  B18 = deg_list_for_venn(deseq_res$B18transfer, 1, 0.05),
  MRLlpr = deg_list_for_venn(deseq_res$AM14MRLlpr, 1, 0.05)
)

DEG_lists_p0.05_p0.01 <- list(
  PL23_2DG_vs_Ctrl = deg_list_for_venn(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, 1, 0.01),
  R848_2DG_vs_Ctrl = deg_list_for_venn(deseq_res$AM14transfer$R848_2DG_v_Ctrl, 1, 0.01),
  PL23_vs_R848 = deg_list_for_venn(deseq_res$AM14transfer$PL23_vs_R848, 1, 0.01),
  B18 = deg_list_for_venn(deseq_res$B18transfer, 1, 0.01),
  MRLlpr = deg_list_for_venn(deseq_res$AM14MRLlpr, 1, 0.01)
)


setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_higher_cutoffs/")
getwd()

# Venn diagrams with padj cutoff of 0.05 ---------------------------------------
  all_degs_0.05 <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG vs PL2-3" = DEG_lists_p0.05$PL23_2DG_vs_Ctrl$all,
                                   "AM14 Transfer:\nR848 + 2DG vs R848" = DEG_lists_p0.05$R848_2DG_vs_Ctrl$all,
                                   "B1-8 Transfer:\nNP + 2DG vs NP" = DEG_lists_p0.05$B18$all,
                                   "AM14 MRL/lpr:\n2DG vs Control" = DEG_lists_p0.05$MRLlpr$all))
  # dev.off()
  plot(all_degs_0.05, filename = "Output/Venn_diagrams/All DEGs - padj 0.05.png",
       cat.cex = 1, cex = 1.5, margin = 0.1)
  
  
  up_degs_0.05 <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG vs PL2-3" = DEG_lists_p0.05$PL23_2DG_vs_Ctrl$up,
                                  "AM14 Transfer:\nR848 + 2DG vs R848" = DEG_lists_p0.05$R848_2DG_vs_Ctrl$up,
                                  "B1-8 Transfer:\nNP + 2DG vs NP" = DEG_lists_p0.05$B18$up,
                                  "AM14 MRL/lpr:\n2DG vs Control" = DEG_lists_p0.05$MRLlpr$up))
  # dev.off()
  plot(up_degs_0.05, filename = "Output/Venn_diagrams/Upregulated DEGs - padj 0.05.png",
       cat.cex = 1, cex = 1.5, margin = 0.1)
  
  down_degs_0.05 <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG vs PL2-3" = DEG_lists_p0.05$PL23_2DG_vs_Ctrl$down,
                                    "AM14 Transfer:\nR848 + 2DG vs R848" = DEG_lists_p0.05$R848_2DG_vs_Ctrl$down,
                                    "B1-8 Transfer:\nNP + 2DG vs NP" = DEG_lists_p0.05$B18$down,
                                    "AM14 MRL/lpr:\n2DG vs Control" = DEG_lists_p0.05$MRLlpr$down))
  # dev.off()
  plot(down_degs_0.05, filename = "Output/Venn_diagrams/Downregulated DEGs - padj 0.05.png",
       cat.cex = 1, cex = 1.5, margin = 0.1)
  
# Venn diagrams with padj cutoff of 0.01 ---------------------------------------
  all_degs_0.01 <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG vs PL2-3" = DEG_lists_p0.01$PL23_2DG_vs_Ctrl$all,
                                   "AM14 Transfer:\nR848 + 2DG vs R848" = DEG_lists_p0.01$R848_2DG_vs_Ctrl$all,
                                   "B1-8 Transfer:\nNP + 2DG vs NP" = DEG_lists_p0.01$B18$all,
                                   "AM14 MRL/lpr:\n2DG vs Control" = DEG_lists_p0.01$MRLlpr$all))
  # dev.off()
  plot(all_degs_0.01, filename = "Output/Venn_diagrams/All DEGs - padj 0.01.png",
       cat.cex = 1, cex = 1.5, margin = 0.1)
  
  
  up_degs_0.01 <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG vs PL2-3" = DEG_lists_p0.01$PL23_2DG_vs_Ctrl$up,
                                  "AM14 Transfer:\nR848 + 2DG vs R848" = DEG_lists_p0.01$R848_2DG_vs_Ctrl$up,
                                  "B1-8 Transfer:\nNP + 2DG vs NP" = DEG_lists_p0.01$B18$up,
                                  "AM14 MRL/lpr:\n2DG vs Control" = DEG_lists_p0.01$MRLlpr$up))
  # dev.off()
  plot(up_degs_0.01, filename = "Output/Venn_diagrams/Upregulated DEGs - padj 0.01.png",
       cat.cex = 1, cex = 1.5, margin = 0.1)
  
  down_degs_0.01 <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG vs PL2-3" = DEG_lists_p0.01$PL23_2DG_vs_Ctrl$down,
                                    "AM14 Transfer:\nR848 + 2DG vs R848" = DEG_lists_p0.01$R848_2DG_vs_Ctrl$down,
                                    "B1-8 Transfer:\nNP + 2DG vs NP" = DEG_lists_p0.01$B18$down,
                                    "AM14 MRL/lpr:\n2DG vs Control" = DEG_lists_p0.01$MRLlpr$down))
  # dev.off()
  plot(down_degs_0.01, filename = "Output/Venn_diagrams/Downregulated DEGs - padj 0.01.png",
       cat.cex = 1, cex = 1.5, margin = 0.1)
  