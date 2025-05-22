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

gseGO_wrapper <- function(dds_res, ID_type, p_cutoff){
  if(missing(ID_type)) ID_type <- "symbol"
  if(missing(p_cutoff)) p_cutoff <- 0.05

  ranked_list <- rank_genes(dds_res, ID_type)
  if(is.null(ranked_list)) {
    print("ranked gene list is NULL (either incorrect ID type or no significant DEGs)")
    return(list(gse = NULL, gse_simplified = NULL, gseGO_summary = NULL))
  } else {
    if(ID_type == "ensembl"){
      gse <- gseGO(ranked_list, ont = "BP", keyType = "ENSEMBL", OrgDb = "org.Mm.eg.db", eps = 1e-300, pvalueCutoff = p_cutoff)
    } else if(ID_type == "entrez"){
      gse <- gseGO(ranked_list, ont = "BP", keyType = "ENTREZID", OrgDb = "org.Mm.eg.db", eps = 1e-300, pvalueCutoff = p_cutoff)
    } else if(ID_type == "symbol"){
      gse <- gseGO(ranked_list, ont = "BP", keyType = "SYMBOL", OrgDb = "org.Mm.eg.db", eps = 1e-300, pvalueCutoff = p_cutoff)
      # gse <- gseGO(ranked_list, ont = "BP", keyType = "SYMBOL", OrgDb = "org.Mm.eg.db", eps = 1e-300, pvalueCutoff = p_cutoff)
    } else return(print("ID_type must be ensembl, entrez, or symbol"))
    
    if(nrow(gse) < 1){
      print("ranked gene list length is > 0 but no significant GO terms were detected")
      return(list(gse = NULL, gse_simplified = NULL, gseGO_summary = NULL))
    } else {
      simp <- simplify(gse, cutoff = 0.7)
      return(list(gse = gse, gse_simplified = simp))
    }
  }
}

gsePathway_wrapper <- function(dds_res){
  ranked_list <- rank_genes(dds_res, "entrez")

  if(is.null(ranked_list)) {
    print("ranked gene list is NULL (no significant DEGs)")
    return(NULL)
  } 
  gse <- gsePathway(ranked_list, organism = "mouse", pvalueCutoff = 0.1, pAdjustMethod = "BH", verbose = FALSE)
  
  if(nrow(gse) > 0) return(gse)
  else {
    print("No significant pathways detected")
    return(NULL)
  }
}

# gseKEGG_wrapper <- function(dds_res){
#   ranked_list <- rank_genes(dds_res, "entrez")
#   
#   if(is.null(ranked_list)) {
#     print("ranked gene list is NULL (no significant DEGs)")
#     return(NULL)
#   } 
#   gse <- gseKEGG(ranked_list, organism = "mmu", pvalueCutoff = 1, 
#                  pAdjustMethod = "BH", verbose = FALSE, keyType = "ncbi-geneid")
#   
#   if(nrow(gse) > 0) return(gse)
#   else {
#     print("No significant pathways detected")
#     return(NULL)
#   }
# }


gseGO_summary <- function(res, comparison, experiment, p_cutoff, nes_cutoff){
  if(missing(p_cutoff)) p_cutoff <- 0.05
  if(missing(nes_cutoff)) nes_cutoff <- 0

  if(is.null(res$gse)) {
    gse_up <- 0
    gse_down <- 0
  } else {
    gse_up <- nrow(filter(res$gse, p.adjust < p_cutoff & NES > nes_cutoff))
    gse_down <- nrow(filter(res$gse, p.adjust < p_cutoff & NES < (-1*nes_cutoff)))
  }
  
  if(is.null(res$gse_simplified)) {
    gseSimp_up <- 0
    gseSimp_down <- 0
  } else {
    gseSimp_up <- nrow(filter(res$gse_simplified, p.adjust < p_cutoff & NES > nes_cutoff))
    gseSimp_down <- nrow(filter(res$gse_simplified, p.adjust < p_cutoff & NES < (-1*nes_cutoff)))
  }
  

  summary <- data.frame(experiment, comparison, nes_cutoff, p_cutoff, 
                        gse_up, gse_down, gseSimp_up, gseSimp_down)
  colnames(summary) <- c("Experiment", "Comparison", "NES_cutoff", "P_adj_cutoff",  
                         "Up (Full)", "Down (Full)", "Up (Simplified)", "Down (Simplified)")
  summary
}

gsePlot <- function(gse, plot_title, file_name, showCat, w, h){
  if(missing(file_name)) file_name <- NULL
  if(missing(showCat)) showCat <- 25
  if(missing(w)) w <- 1700
  if(missing(h)) h <- 2000

  if(nrow(gse) < showCat) showCat <- nrow(gse)
  
  gse_top <- gse[1:showCat, ]
  nes_max <- max(abs(gse_top$NES))

  if(!is.null(file_name)){
    png(filename = stri_join(c("Output/Functional_analyses/", file_name, ".png"), collapse = ""),
        width = w, height = h, units = "px", pointsize = 10, res = 200,
        bg = "white", family = "", symbolfamily="default")
  }
  ggplot(gse, showCategory = showCat, 
         aes(p.adjust, fct_reorder(Description, p.adjust))) +
    geom_segment(aes(xend=0, yend = Description)) +
    geom_point(aes(color=NES, size = setSize)) +
    scale_size_continuous(range=c(2, 10)) +
    theme_dose(12) +
    xlab("FDR") +
    ylab(NULL) +
    scale_y_discrete(labels = label_wrap(50)) +
    # scale_y_discrete(labels = lapply(strwrap(gse$Description, width = 45, simplify = FALSE), paste, collapse="\n")) +
    ggtitle(plot_title) +
    scale_color_gradientn(colors = pnw_palette("Bay", 8, type = "continuous"), limits = c(-1*nes_max, nes_max))
    # scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous"), limits = c(-1*nes_max, nes_max))
    # scale_color_gradientn(colors = met.brewer("Hiroshige"), limits = c(-1*nes_max, nes_max))
}

# GO Gene set enrichment analysis ----------------------------------------------
  gseGO_results <- list(
    AM14trans = list(
      PL23_2DG_v_Ctrl = gseGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "symbol", 0.2),
      R848_2DG_v_Ctrl = gseGO_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl, "symbol", 0.2),
      PL23_vs_R848 = gseGO_wrapper(deseq_res$AM14transfer$PL23_vs_R848, "symbol", 0.2)),
    B18trans = gseGO_wrapper(deseq_res$B18transfer, "symbol", 0.2),
    PL23_v_NP = gseGO_wrapper(deseq_res$PL23_v_NP, "symbol", 0.2),
    AM14MRLlpr = gseGO_wrapper(deseq_res$AM14MRLlpr, "symbol", 0.2))

  gseGO_res_summary <- join_all(list(gseGO_summary(gseGO_results$AM14trans$PL23_2DG_v_Ctrl, "PL2-3+2DG vs PL2-3", "AM14 Adoptive Transfer"),
                                     gseGO_summary(gseGO_results$AM14trans$R848_2DG_v_Ctrl, "R848+2DG vs R848", "AM14 Adoptive Transfer"),
                                     gseGO_summary(gseGO_results$AM14trans$PL23_vs_R848, "PL2-3 vs R848", "AM14 Adoptive Transfer"),
                                     gseGO_summary(gseGO_results$B18trans, "NP+2DG vs NP", "B1-8 Adoptive Transfer"),
                                     gseGO_summary(gseGO_results$PL23_v_NP, "PL2-3 vs NP", "AM14 and B1-8 Adoptive Transfers"),
                                     gseGO_summary(gseGO_results$AM14MRLlpr, "2DG vs Control", "AM14 MRL/lpr")), 
                                type = "full")
  gseGO_res_summary 

  # Save results to an Excel file ----------------------------------------------
    wb <- createWorkbook("Output/Functional_analyses/gsea_GOBP_results.xlsx")
    
    addWorksheet(wb, "GSE GO Summaries")
    addWorksheet(wb, "AM14 PL23_2DG_vs_Ctrl - Full")
    # addWorksheet(wb, "AM14 - R848_2DG_vs_Ctrl - Full")
    # addWorksheet(wb, "AM14 - PL23_vs_R848 - Full")
    # addWorksheet(wb, "B1-8 - 2DG_vs_Ctrl - Full")
    addWorksheet(wb, "PL2-3 vs NP - Full")
    addWorksheet(wb, "AM14 MRLlpr - Full")
    
    addWorksheet(wb, "AM14 PL23_2DG_vs_Ctrl - Simp")
    # addWorksheet(wb, "AM14 - R848_2DG_vs_Ctrl - Simp")
    # addWorksheet(wb, "AM14 - PL23_vs_R848 - Simp")
    # addWorksheet(wb, "B1-8 - 2DG_vs_Ctrl - Simp")
    addWorksheet(wb, "PL2-3 vs NP - Simp")
    addWorksheet(wb, "AM14 MRLlpr - Simp")
    
    writeData(wb, "GSE GO Summaries", gseGO_res_summary)
    writeData(wb, "AM14 PL23_2DG_vs_Ctrl - Full", gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse)
    # writeData(wb, "AM14 - R848_2DG_vs_Ctrl - Full", gseGO_results$AM14trans$R848_2DG_v_Ctrl$gse)
    # writeData(wb, "AM14 - PL23_vs_R848 - Full", gseGO_results$AM14trans$PL23_vs_R848$gse)
    # writeData(wb, "B1-8 - 2DG_vs_Ctrl - Full", gseGO_results$B18trans$gse)
    writeData(wb, "PL2-3 vs NP - Full", gseGO_results$PL23_v_NP$gse)
    writeData(wb, "AM14 MRLlpr - Full", gseGO_results$AM14MRLlpr$gse)
    
    writeData(wb, "AM14 PL23_2DG_vs_Ctrl - Simp", gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse_simplified)
    # writeData(wb, "AM14 - R848_2DG_vs_Ctrl - Simp", gseGO_results$AM14trans$R848_2DG_v_Ctrl$gse_simplified)
    # writeData(wb, "AM14 - PL23_vs_R848 - Simp", gseGO_results$AM14trans$PL23_vs_R848$gse_simplified)
    # writeData(wb, "B1-8 - 2DG_vs_Ctrl - Simp", gseGO_results$B18trans$gse_simplified)
    writeData(wb, "PL2-3 vs NP - Simp", gseGO_results$PL23_v_NP$gse_simplified)
    writeData(wb, "AM14 MRLlpr - Simp", gseGO_results$AM14MRLlpr$gse_simplified)
    
    saveWorkbook(wb, "Output/Functional_analyses/gsea_GOBP_results.xlsx", overwrite = TRUE)
    rm(wb)
    

  # visualize results ----------------------------------------------------------
    # gseaplot(, geneSetID = 1, 
    #          title = gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse$Description[1])
    # gseaplot2(gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse, geneSetID = 1, 
    #           title = gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse$Description[1])
    # gseaplot2(gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse, geneSetID = 1:4)
    # head(gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse)
    
    gsePlot(gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse_simplified,
            "AM14 Transfer:\nPL2-3+2DG vs PL2-3\n(GO Biological Process)", showCat = 26,
            "gseGO_plots/AM14 Transfer - PL2-3 + 2DG vs PL2-3", w = 1400, h = 2000)
    dev.off()

    # gsePlot(gseGO_results$PL23_v_NP$gse_simplified,
    #         "AM14 Transfer + PL2-3\nvs B1-8 Transfer + NP\n(GO Biological Process)",
    #         "gseGO_plots/AM14 Transfer PL2-3 vs B1-8 Transfer NP - UP")
    
    
    temp_PL23_NP_up <- dplyr::filter(gseGO_results$PL23_v_NP$gse_simplified, NES > 0)
    max(temp_PL23_NP_up$NES)
    gsePlot(temp_PL23_NP_up,
            "AM14 Transfer + PL2-3\nvs B1-8 Transfer + NP\n(Upregulated GO Biological Process)", 
            "gseGO_plots/AM14 Transfer PL2-3 vs B1-8 Transfer NP - UP")
    dev.off()
    rm(temp_PL23_NP_up)
    
    temp_PL23_NP_down <- dplyr::filter(gseGO_results$PL23_v_NP$gse_simplified, NES < 0)
    nrow(temp_PL23_NP_down)
    gsePlot(temp_PL23_NP_down,
            "AM14 Transfer + PL2-3\nvs B1-8 Transfer + NP\n(Downregulated GO Biological Process)",
            "gseGO_plots/AM14 Transfer PL2-3 vs B1-8 Transfer NP - DOWN", w = 1200, h = 1000)
    dev.off()
    rm(temp_PL23_NP_down)
    
    
    gsePlot(gseGO_results$AM14MRLlpr$gse, 
            "AM14 MRL/lpr: 2DG vs Control\n(GO Biological Process)", 
            "gseGO_plots/AM14 MRLlpr 2DG vs Control", h = 1500)
    dev.off()
    

# Reactome pathways ------------------------------------------------------------
  gsePathway_results <- list(
    AM14trans = list(
      PL23_2DG_v_Ctrl = gsePathway_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl),
      R848_2DG_v_Ctrl = gsePathway_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl),
      PL23_vs_R848 = gsePathway_wrapper(deseq_res$AM14transfer$PL23_vs_R848)),
    B18trans = gsePathway_wrapper(deseq_res$B18transfer),
    # PL23_v_NP = gsePathway_wrapper(deseq_res$PL23_v_NP),
    AM14MRLlpr = gsePathway_wrapper(deseq_res$AM14MRLlpr)
  )
  nrow(gsePathway_results$AM14trans$PL23_2DG_v_Ctrl)
  nrow(gsePathway_results$AM14trans$R848_2DG_v_Ctrl)
  nrow(gsePathway_results$AM14trans$PL23_vs_R848)
  nrow(gsePathway_results$B18trans)
  # nrow(gsePathway_results$PL23_v_NP)
  nrow(gsePathway_results$AM14MRLlpr)
  
  head(gsePathway_results$AM14trans$PL23_2DG_v_Ctrl)
  gsePlot(gsePathway_results$AM14trans$PL23_2DG_v_Ctrl, 
          "AM14 Transfer: PL2-3+2DG vs PL2-3\n(Reactome)",
          "gsePathway_plots/AM14 Transfer - PL2-3 + 2DG vs PL2-3", w = 1200, h = 750)
  dev.off()
  
  gsePlot(gsePathway_results$AM14trans$PL23_vs_R848, 
          "AM14 Transfer: PL2-3 vs R848\n(Reactome)",
          "gsePathway_plots/AM14 Transfer - PL2-3 vs R848", w = 1200, h = 750)
  dev.off()
  
  # gsePlot(gsePathway_results$PL23_v_NP,
  #            "AM14 Transfer + PL2-3 vs B1-8 Transfer + NP\n(Reactome)",
  #            "gsePathway_plots/AM14 Transfer PL2-3 vs B1-8 Transfer NP", w = 1700, h = 2000)
  # dev.off()
  # 
  # temp_PL23_NP_up <- dplyr::filter(gsePathway_results$PL23_v_NP, NES > 0)
  # nrow(temp_PL23_NP_up)
  # gsePlot(temp_PL23_NP_up,
  #         "AM14 Transfer + PL2-3 vs B1-8 Transfer + NP\n(Upregulated Reactome Pathways)", 
  #         "gsePathway_plots/AM14 Transfer PL2-3 vs B1-8 Transfer NP - UP")
  # dev.off()
  # rm(temp_PL23_NP_up)
  # 
  # temp_PL23_NP_down <- dplyr::filter(gsePathway_results$PL23_v_NP, NES < 0)
  # nrow(temp_PL23_NP_down)
  # gsePlot(temp_PL23_NP_down,
  #         "AM14 Transfer + PL2-3 vs B1-8 Transfer + NP\n(Downregulated Reactome Pathways)",
  #         "gsePathway_plots/AM14 Transfer PL2-3 vs B1-8 Transfer NP - DOWN", w = 1700, h = 1600)
  # dev.off()
  # rm(temp_PL23_NP_down)
  
  gsePlot(gsePathway_results$AM14MRLlpr, 
             "AM14 MRL/lpr: 2DG vs Control\n(Reactome)",
             "gsePathway_plots/AM14 MRLlpr 2DG vs Control", w = 1200, h = 700)
  dev.off()
  
# head(pathway)
# viewPathway("Signaling by GPCR", readable = TRUE, foldChange = genes_sorted)


# KEGG pathways ------------------------------------------------------------
#   gseKEGG_results <- list(
#     AM14trans = list(
#       PL23_2DG_v_Ctrl = gseKEGG_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl),
#       R848_2DG_v_Ctrl = gseKEGG_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl),
#       PL23_vs_R848 = gseKEGG_wrapper(deseq_res$AM14transfer$PL23_vs_R848)),
#     B18trans = gseKEGG_wrapper(deseq_res$B18transfer),
#     PL23_v_NP = gseKEGG_wrapper(deseq_res$PL23_v_NP),
#     AM14MRLlpr = gseKEGG_wrapper(deseq_res$AM14MRLlpr)
#   )
