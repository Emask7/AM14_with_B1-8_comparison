# KEGG pathway over-representation analysis ------------------------------------
# see https://www.genome.jp/kegg/ for information about KEGG
enrichKEGG_wrapper <- function(dds_res){
  sig <- dplyr::filter(dds_res, padj < 0.05)
  sig <- as.character(sig$entrezgene)
  enrichKEGG(gene = sig, organism = 'mmu', pvalueCutoff = 0.05)
}

enrichKEGG_results <- list(
  AM14trans = list(
    PL23_2DG_v_Ctrl = enrichKEGG_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl),
    R848_2DG_v_Ctrl = enrichKEGG_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl),
    PL23_vs_R848 = enrichKEGG_wrapper(deseq_res$AM14transfer$PL23_vs_R848)),
  B18trans = enrichKEGG_wrapper(deseq_res$B18transfer),
  PL23_v_NP = enrichKEGG_wrapper(deseq_res$PL23_v_NP),
  AM14MRLlpr = enrichKEGG_wrapper(deseq_res$AM14MRLlpr)
)

nrow(enrichKEGG_results$AM14trans$PL23_2DG_v_Ctrl)
nrow(enrichKEGG_results$AM14trans$R848_2DG_v_Ctrl)
nrow(enrichKEGG_results$AM14trans$PL23_vs_R848)
nrow(enrichKEGG_results$B18trans)
head(enrichKEGG_results$PL23_v_NP)
nrow(enrichKEGG_results$AM14MRLlpr)

# ## visualize enriched KEGG pathways
# browseKEGG(kk, 'hsa04510')
# 
# lg2FC <- sig$log2FoldChange
# names(lg2FC) <- sig$ENTREZID
# head(lg2FC)
# 
# hsa04510 <- pathview(gene.data  = lg2FC,
#                      pathway.id = "hsa04510",
#                      species    = "hsa",
#                      limit      = list(gene=5, cpd=1))
# 
# # set limit to 3 for optimal color scale - change this as needed
