deg_list_for_venn <- function(res, lfc_cutoff){
  list(up = subset(res, res$padj <= 0.05 & res$log2FoldChange >= lfc_cutoff)[, 1], 
       down = subset(res, res$padj <= 0.05 & res$log2FoldChange <= (-1*lfc_cutoff))[, 1], 
       all = subset(res, res$padj <= 0.05 & abs(res$log2FoldChange) >= lfc_cutoff)[, 1])
}


DEG_lists <- list(
  PL23_2DG_vs_Ctrl = deg_list_for_venn(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, 1),
  R848_2DG_vs_Ctrl = deg_list_for_venn(deseq_res$AM14transfer$R848_2DG_v_Ctrl, 1),
  PL23_vs_R848 = deg_list_for_venn(deseq_res$AM14transfer$PL23_vs_R848, 1),
  B18 = deg_list_for_venn(deseq_res$B18transfer, 1),
  MRLlpr = deg_list_for_venn(deseq_res$AM14MRLlpr, 1)
)

venn_up_1 <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG\nvs PL2-3\n" = DEG_lists$PL23_2DG_vs_Ctrl$up,
                             "AM14 Transfer:\nR848 + 2DG\nvs R848\n" = DEG_lists$R848_2DG_vs_Ctrl$up,
                             "AM14 MRL/lpr:\n2DG vs Control\n" = DEG_lists$MRLlpr$up,
                             "B1-8 Transfer:\nNP+2DG vs NP\n" = DEG_lists$B18$up))
plot(venn_up_1, mycol = met.brewer("Johnson", n = 8, direction = 1),
     filename = "Output/venn_diagram_upreg.png",
     margin = 0.05, cat.cex = 1, cex = 1.5)
dev.off()
detail(venn_up_1)


venn_up_2 <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG\nvs PL2-3\n" = DEG_lists$PL23_2DG_vs_Ctrl$up,
                             "AM14 Transfer:\nR848 + 2DG\nvs R848\n" = DEG_lists$R848_2DG_vs_Ctrl$up,
                             "AM14 Transfer:\nPL2-3\nvs R848\n" = DEG_lists$PL23_vs_R848$up,
                             "AM14 MRL/lpr:\n2DG vs Control\n" = DEG_lists$MRLlpr$up,
                             "B1-8 Transfer:\nNP+2DG vs NP\n" = DEG_lists$B18$up))
venn_up_2 <- venndetail(list("AM14: PL2-3+2DG vs PL2-3" = DEG_lists$PL23_2DG_vs_Ctrl$up,
                             "AM14: R848+2DG vs R848" = DEG_lists$R848_2DG_vs_Ctrl$up,
                             "AM14: PL2-3 vs R848" = DEG_lists$PL23_vs_R848$up,
                             "AM14 MRL/lpr: 2DG vs Control" = DEG_lists$MRLlpr$up,
                             "B1-8: NP+2DG vs NP" = DEG_lists$B18$up),
                        sep = " X ")
detail(venn_up_2)

plot(venn_up_2, mycol = met.brewer("Johnson", n = 8, direction = 1),
     # filename = "Output/venn_diagram_upreg.png",
     margin = 0.05, cat.cex = 1, cex = 1.5)
dev.off()


# dev.off()
# plot(venn_up_1, mycol = moma.colors("Panton", 8, direction = -1),
#      # filename = "Output/venn_diagram_upreg.png",
#      margin = 0.05, cat.cex = 1, cex = 1.5)
# 
# dev.off()
# plot(venn_up_1, mycol = moma.colors("ustwo", 8, direction = 1),
#      # filename = "Output/venn_diagram_upreg.png",
#      margin = 0.05, cat.cex = 1, cex = 1.5)
# 
# dev.off()
# plot(venn_up_1, mycol = moma.colors("OKeeffe", 8, direction = 1),
#      # filename = "Output/venn_diagram_upreg.png",
#      margin = 0.05, cat.cex = 1, cex = 1.5)


venn_down_1 <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG\nvs PL2-3\n" = DEG_lists$PL23_2DG_vs_Ctrl$down,
                               "AM14 Transfer:\nR848 + 2DG\nvs R848\n" = DEG_lists$R848_2DG_vs_Ctrl$down,
                               "AM14 MRL/lpr:\n2DG vs Control\n" = DEG_lists$MRLlpr$down,
                               "B1-8 Transfer:\nNP + 2DG vs NP\n" = DEG_lists$B18$down))
plot(venn_down_1,  mycol = met.brewer("Johnson", n = 8, direction = -1),
     filename = "Output/venn_diagram_downreg.png",
     margin = 0.05, cat.cex = 1, cex = 1.5)
dev.off()

detail(venn_down_1)
