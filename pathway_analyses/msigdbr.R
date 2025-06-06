install.packages("msigdbr")
library(msigdbr)
hallmark_genesets <- msigdbr(species = "mouse", collection = "H")
head(hallmark_genesets)
gobp_genesets <- msigdbr(species = "mouse", collection = "C5", subcollection = "BP")
head(gobp_genesets)

rank_genes <- function(dds_res, ID_type){
  if(missing(ID_type)) ID_type <- "ensembl"
  
  sig <- dplyr::filter(dds_res, padj < 0.05)
  sig_ordered <- sig[sig$baseMean > 50,]
  sig_ordered <- sig_ordered[order(-sig_ordered$log2FoldChange), ]
  gene_list <- sig_ordered$log2FoldChange
  if(length(gene_list) <= 0) {
    print("gene list length <= 0")
    return(NULL)
  }
  
  if(ID_type == "ensembl") names(gene_list) <- sig_ordered$ensembl_gene_id
  else if(ID_type == "entrez") names(gene_list) <- sig_ordered$entrezgene
  else if(ID_type == "symbol") names(gene_list) <- sig_ordered$external_gene_name
  else {
    print("ID_type must be ensembl, entrez, or symbol")
    return(NULL)
  }
  
  gene_list <- subset(gene_list, !duplicated(names(gene_list)))
  gene_list <- subset(gene_list, names(gene_list) != "NA")
  gene_list <- na.omit(gene_list)
  return(gene_list)
}



ranked_genes_PL23_2DG_v_Ctrl <- rank_genes(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "symbol")
head(ranked_genes_PL23_2DG_v_Ctrl)

msigdbr_t2g <- dplyr::distinct(hallmark_genesets, gs_name, gene_symbol)
head(msigdbr_t2g)
res_PL23_2DG_v_Ctrl <- enricher(gene = names(ranked_genes_PL23_2DG_v_Ctrl), TERM2GENE = msigdbr_t2g, qvalueCutoff = 0.2)
res_PL23_2DG_v_Ctrl
data.frame(res_PL23_2DG_v_Ctrl)


msigdbr_list <- split(x = hallmark_genesets$gene_symbol, f = hallmark_genesets$gs_name)
head(msigdbr_list)
fgsea_PL23_2DG_v_Ctrl <- fgsea(pathways = msigdbr_list, stats = ranked_genes_PL23_2DG_v_Ctrl, minSize = 8)
data.frame(fgsea_PL23_2DG_v_Ctrl)
subset(data.frame(fgsea_PL23_2DG_v_Ctrl), padj < 0.1)
