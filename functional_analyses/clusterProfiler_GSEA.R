rank_genes <- function(dds_res, ID_type){
  if(missing(ID_type)) ID_type <- "ensembl"
  
  sig <- dplyr::filter(dds_res, padj < 0.05)
  sig_ordered <- sig[sig$baseMean > 50,]
  sig_ordered <- sig_ordered[order(-sig_ordered$log2FoldChange), ]
  gene_list <- sig_ordered$log2FoldChange
  if(length(gene_list) <= 0) return(NULL)
  
  if(ID_type == "ensembl") names(gene_list) <- sig_ordered$ensembl_gene_id
  else if(ID_type == "entrez") names(gene_list) <- sig_ordered$entrezgene
  else if(ID_type == "symbol") names(gene_list) <- sig_ordered$external_gene_name
  else {
    print("ID_type must be ensembl, entrez, or symbol")
    return(NULL)
  }
  
  gene_list <- subset(gene_list, !duplicated(gene_list))
  gene_list <- subset(gene_list, names(gene_list) != "NA")
  gene_list <- na.omit(gene_list)
  return(gene_list)
}

gseGO_wrapper <- function(dds_res, ID_type){
  if(missing(ID_type)) ID_type <- "ensembl"
  
  ranked_list <- rank_genes(dds_res, ID_type)
  if(is.null(ranked_list)) {
    print("ranked gene list is NULL (either incorrect ID type or no significant DEGs)")
    return(NULL)
  } else {
    if(ID_type == "ensembl"){
      gse <- gseGO(ranked_list, ont = "BP", keyType = "ENSEMBL", OrgDb = "org.Mm.eg.db", eps = 1e-300)
    } else if(ID_type == "entrez"){
      gse <- gseGO(ranked_list, ont = "BP", keyType = "ENTREZID", OrgDb = "org.Mm.eg.db", eps = 1e-300)
    } else if(ID_type == "symbol"){
      gse <- gseGO(ranked_list, ont = "BP", keyType = "SYMBOL", OrgDb = "org.Mm.eg.db", eps = 1e-300)
    } else return(print("ID_type must be ensembl, entrez, or symbol"))
    
    simp <- simplify(gse, cutoff = 0.7)
    
    return(list(gse = gse, gse_simplified = simp, gse_summary = as.data.frame(gse)))
  }
}

plot_gseGO <- function(gse, plot_title, showCat){
  if(missing(showCat)) showCat <- nrow(gse)
  else if(nrow(gse) < showCat) showCat <- nrow(gse)
  
  # gse_top <- gse[1:showCat]
  # nes_max <- max(gse_top$NES)
  # nes_min <- min(gse_top$NES)
  
  ggplot(gse, showCategory = showCat,
         aes(p.adjust, Description)) +
    geom_segment(aes(xend=0, yend = Description)) +
    geom_point(aes(color=NES, size = setSize)) +
    scale_size_continuous(range=c(2, 10)) +
    theme_dose(12) +
    xlab("FDR") +
    ylab(NULL) +
    scale_y_discrete(labels = lapply(strwrap(gse$Description, width = 70, simplify = FALSE), paste, collapse="\n")) +
    ggtitle(plot_title) +
    scale_color_manual(values = pnw_palette("Bay", 8, type = "continuous"))
                          #guide=guide_colorbar(reverse=FALSE, order=1))
  
    # if(nes_max < 0) {
    #   scale_color_gradientn(colours=met.brewer("Hokusai2", direction = -1), 
    #                         guide=guide_colorbar(reverse=FALSE, order=1))
    # } else if(nes_min > 0) {
    #   scale_color_gradientn(colours=met.brewer("Tam", direction = 1), 
    #                         guide=guide_colorbar(reverse=FALSE, order=1))
    # } else {
    #   scale_color_gradientn(colours=met.brewer("Benedictus", direction = -1), 
    #                         guide=guide_colorbar(reverse=FALSE, order=1))
    # }
  
  # ggplot(gse, showCategory = showCat,
  #        aes(NES, fct_reorder(Description, NES), fill=p.adjust)) +
  #   geom_col() +
  #   scale_fill_gradientn(colours=met.brewer("Tam", direction = 1), guide=guide_colorbar(reverse=TRUE)) +
  #   theme_dose(12) +
  #   xlab("Normalized Enrichment Score") +
  #   ylab(NULL) +
  #   scale_y_discrete(labels = lapply(strwrap(gse$Description, width = 70, simplify = FALSE), paste, collapse="\n")) +
  #   ggtitle(plot_title)
  # 
  # ggplot(gse, showCategory = showCat,
  #        aes(NES, fct_reorder(Description, NES))) +
  #   geom_segment(aes(xend=0, yend = Description)) +
  #   geom_point(aes(color=p.adjust, size = setSize)) +
  #   scale_color_gradientn(colours=met.brewer("Tam", direction = 1), guide=guide_colorbar(reverse=TRUE, order=1)) +
  #   scale_size_continuous(range=c(2, 10)) +
  #   theme_dose(12) +
  #   xlab("Normalized Enrichment Score") +
  #   ylab(NULL) +
  #   scale_y_discrete(labels = lapply(strwrap(gse$Description, width = 70, simplify = FALSE), paste, collapse="\n")) +
  #   ggtitle(plot_title)
}

gsePathway_wrapper <- function(dds_res){
  ranked_list <- rank_genes(dds_res, "entrez")

  if(is.null(ranked_list)) {
    print("ranked gene list is NULL (no significant DEGs)")
    return(NULL)
  } else {
    gsePathway(ranked_list, organism = "mouse", pvalueCutoff = 0.1, pAdjustMethod = "BH", verbose = FALSE)
  }
}

# ranked_genelist <- rank_genes(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "ensembl")
# head(ranked_genelist)


# GO Gene set enrichment analysis ----------------------------------------------
  # # testing gene ID types --------------------------------------------------------
  # temp_symbol <- gseGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "symbol")
  # temp_ensembl <- gseGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "ensembl")
  # temp_entrez <- gseGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "entrez")
  # 
  # nrow(temp$gse_summary)
  # temp$gse_summary[, 2:8]
  # 
  # nrow(temp_symbol$gse_summary)
  # temp_symbol$gse_summary[, 2:8]
  # 
  # nrow(temp_ensembl$gse_summary)
  # temp_ensembl$gse_summary[, 2:8]
  # 
  # nrow(temp_entrez$gse_summary)
  # temp_entrez$gse_summary[, 2:8]


  # list object with results -----------------------------------------------------
    gseGO_results <- list(
      AM14trans = list(
        PL23_2DG_v_Ctrl = gseGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl),
        R848_2DG_v_Ctrl = gseGO_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl),
        PL23_vs_R848 = gseGO_wrapper(deseq_res$AM14transfer$PL23_vs_R848)),
      B18trans = gseGO_wrapper(deseq_res$B18transfer),
      PL23_v_NP = gseGO_wrapper(deseq_res$PL23_v_NP),
      AM14MRLlpr = gseGO_wrapper(deseq_res$AM14MRLlpr)
    )
    nrow(gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse_summary)
    nrow(gseGO_results$AM14trans$R848_2DG_v_Ctrl$gse_summary)
    nrow(gseGO_results$AM14trans$PL23_vs_R848$gse_summary)
    nrow(gseGO_results$B18trans$gse_summary)
    nrow(gseGO_results$PL23_v_NP$gse_summary)
    nrow(gseGO_results$AM14MRLlpr$gse_summary)

  # visualize results ------------------------------------------------------------
    # gseaplot(, geneSetID = 1, 
    #          title = gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse$Description[1])
    # gseaplot2(gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse, geneSetID = 1, 
    #           title = gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse$Description[1])
    # gseaplot2(gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse, geneSetID = 1:4)
    # head(gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse)
    
    
    gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse$p.adjust
    
    plot_gseGO(gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse_simplified, 
               "AM14 Transfer:\nPL2-3+2DG vs PL2-3\n(GO Biological Process)")
    
    plot_gseGO(gseGO_results$PL23_v_NP$gse_simplified, 
               "AM14 Transfer + PL2-3 vs B1-8 Transfer + NP\n(GO Biological Process)", 25)
    
    plot_gseGO(gseGO_results$AM14MRLlpr$gse, 
               "AM14 MRL/lpr: 2DG vs Control\n(GO Biological Process)")
    

# Reactome pathways ------------------------------------------------------------
  gsePathway_results <- list(
    AM14trans = list(
      PL23_2DG_v_Ctrl = gsePathway_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl),
      R848_2DG_v_Ctrl = gsePathway_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl),
      PL23_vs_R848 = gsePathway_wrapper(deseq_res$AM14transfer$PL23_vs_R848)),
    B18trans = gsePathway_wrapper(deseq_res$B18transfer),
    PL23_v_NP = gsePathway_wrapper(deseq_res$PL23_v_NP),
    AM14MRLlpr = gsePathway_wrapper(deseq_res$AM14MRLlpr)
  )
  nrow(gsePathway_results$AM14trans$PL23_2DG_v_Ctrl)
  nrow(gsePathway_results$AM14trans$R848_2DG_v_Ctrl)
  nrow(gsePathway_results$AM14trans$PL23_vs_R848)
  nrow(gsePathway_results$B18trans)
  nrow(gsePathway_results$PL23_v_NP)
  nrow(gsePathway_results$AM14MRLlpr)


# head(pathway)
# viewPathway("Signaling by GPCR", readable = TRUE, foldChange = genes_sorted)



