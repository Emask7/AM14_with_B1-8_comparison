pathway_list_for_venn <- function(res, lfc_cutoff){
  list(up = subset(Description , res$padj <= 0.05 & res$log2FoldChange >= lfc_cutoff)[, 1], 
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
all_GOBP <- venndetail(list("AM14 Transfer: PL2-3 + 2DG vs PL2-3" = enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified$Description,
                            "AM14 Transfer: R848 + 2DG vs R848" = enrichGO_all$AM14trans$R848_2DG_v_Ctrl$eGO_simplified$Description,
                            "AM14 Transfer: PL2-3 vs R848" = enrichGO_all$AM14trans$PL23_vs_R848$eGO_simplified$Description,
                            # "B1-8 Transfer: NP + 2DG vs NP" = enrichGO_all$B18trans,
                            "AM14 MRL/lpr: 2DG vs Control" = enrichGO_all$AM14MRLlpr$eGO_simplified$Description))
png(filename = "Output/Venn_diagrams/All_GOBP.png", width = 2000, height = 1500, res = 250)
plot(all_GOBP, type = "upset")
# plot(all_GOBP)
dev.off()
