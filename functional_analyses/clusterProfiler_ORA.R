enrichGO_wrapper <- function(dds_res, lfc_direction, lfc_cutoff, GO_padj_cutoff){
  if(missing(lfc_direction)) lfc_direction <- NULL
  if(missing(lfc_cutoff)) lfc_cutoff <- 1
  if(missing(GO_padj_cutoff)) GO_padj_cutoff <- 0.05
  
  background <- as.character(dds_res$ensembl_gene_id)
  
  if(is.null(lfc_direction)) {
    sig <- dplyr::filter(dds_res, padj < 0.05 & abs(log2FoldChange) > abs(lfc_cutoff))
  } else if(lfc_direction == "up"){
    sig <- dplyr::filter(dds_res, padj < 0.05 & log2FoldChange > lfc_cutoff)
  } else if(lfc_direction == "down"){
    lfc_cutoff <- -1*abs(lfc_cutoff)
    sig <- dplyr::filter(dds_res, padj < 0.05 & log2FoldChange < lfc_cutoff)
  } else return(print("lfc_direction must be either up, down, or NULL"))
  
  if(nrow(sig) < 1) {
    print("0 significant DEGs")
    return(NULL)
  }
  
  sig <- as.character(sig$ensembl_gene_id)

  eGO <- enrichGO(gene = sig, universe = background, keyType = "ENSEMBL", qvalueCutoff = GO_padj_cutoff,
                  OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)

  if(nrow(eGO) < 1) {
    print("ranked gene list length is > 0 but no significant GO terms were detected")
    return(NULL)
  } else {
    simp <- simplify(eGO, cutoff = 0.7)
    return(list(eGO = eGO, eGO_simplified = simp))
  }
}

enrichGO_summary <- function(res, comparison, experiment, DEG_LFC_cutoff){
  if(missing(DEG_LFC_cutoff)) DEG_LFC_cutoff <- "NA"
  
  if(is.null(res$eGO)) eGO_sig <- 0
  else eGO_sig <- nrow(filter(res$eGO, p.adjust < 0.05))

  if(is.null(res$eGO_simplified)) eGOSimp_sig <- 0
  else eGOSimp_sig <- nrow(filter(res$eGO_simplified, p.adjust < 0.05))
  
  summary <- data.frame(experiment, comparison, 0.05, DEG_LFC_cutoff, eGO_sig, eGOSimp_sig)
  colnames(summary) <- c("Experiment", "Comparison", "P_adj_cutoff", "DEG_LFC_cutoff", "Significant (Full)", "Significant (Simplified)")
  summary
}
  

oraPlot <- function(ora, plot_title, wrap, file_name, showCat, w, h){
  if(missing(wrap)) wrap <- 40
  if(missing(file_name)) file_name <- NULL
  if(missing(showCat)) showCat <- 25
  if(missing(w)) w <- 1250
  if(missing(h)) h <- 2200
  
  if(nrow(ora) < showCat) showCat <- nrow(ora)
  
  plotting_data <- as.data.frame(ora)
  plotting_data <- plotting_data[order(-plotting_data$RichFactor), ]
  plotting_data <- plotting_data[1:showCat, ]
  
  if(!is.null(file_name)){
    png(filename = stri_join(c("Output/Functional_analyses/", file_name, ".png"), collapse = ""),
        width = w, height = h, units = "px", pointsize = 10, res = 200,
        bg = "white", family = "", symbolfamily="default")
  }
  
  ggplot(plotting_data, #showCategory = showCat,
         aes(RichFactor, fct_reorder(Description, RichFactor))) +
         # aes(p.adjust, fct_reorder(Description, p.adjust))) +
    geom_segment(aes(xend=0, yend = Description)) +
    geom_point(aes(color=p.adjust, size = Count)) +
    # geom_point(aes(color=RichFactor, size = Count)) +
    # scale_color_gradientn(colors = pnw_palette("Bay", 8, type = "continuous")) +
    scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous"),
                          transform = "reverse",
                          guide = guide_colorbar(reverse = TRUE)) +
    scale_size_continuous(range=c(3, 10)) +
    theme_dose(12) +
    xlab("RichFactor") +
    # xlab("FDR") +
    ylab(NULL) +
    scale_y_discrete(labels = label_wrap(wrap)) +
    ggtitle(plot_title)
}
  
# GO Enrichment Analysis -------------------------------------------------------
  # Up- and down-regulated DEGs ------------------------------------------------
    enrichGO_all <- list(
      AM14trans = list(
        PL23_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, NULL, 0.5),
        R848_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl, NULL, 0.5),
        PL23_vs_R848 = enrichGO_wrapper(deseq_res$AM14transfer$PL23_vs_R848, NULL, 0.5)),
      PL23_v_NP = enrichGO_wrapper(deseq_res$PL23_v_NP, NULL, 0.5),
      B18trans = enrichGO_wrapper(deseq_res$B18transfer, NULL, 0.5),
      AM14MRLlpr = enrichGO_wrapper(deseq_res$AM14MRLlpr, NULL, 0.5))
  
    head(enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified)
    enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified$Description
  
    enrichGO_all_summary <- join_all(list(enrichGO_summary(enrichGO_all$AM14trans$PL23_2DG_v_Ctrl, "PL2-3+2DG vs PL2-3", "AM14 Adoptive Transfer", "+/- 0.5"),
                                          enrichGO_summary(enrichGO_all$AM14trans$R848_2DG_v_Ctrl, "R848+2DG vs R848", "AM14 Adoptive Transfer", "+/- 0.5"),
                                          enrichGO_summary(enrichGO_all$AM14trans$PL23_vs_R848, "PL2-3 vs R848", "AM14 Adoptive Transfer", "+/- 0.5"),
                                          enrichGO_summary(enrichGO_all$B18trans, "NP+2DG vs NP", "B1-8 Adoptive Transfer", "+/- 0.5"),
                                          enrichGO_summary(enrichGO_all$PL23_v_NP, "PL2-3 vs NP", "AM14 and B1-8 Adoptive Transfers", "+/- 0.5"),
                                          enrichGO_summary(enrichGO_all$AM14MRLlpr, "2DG vs Control", "AM14 MRL/lpr", "+/- 0.5")), 
                                     type = "full")
    enrichGO_all_summary
    
    
    # test <- enrichGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "down", 0.5)
    # enrichGO_summary(test, "PL2-3+2DG vs PL2-3", "AM14 Adoptive Transfer", " ")
    # test$eGO_simplified$Description
    # rm(test)
    
    # test <- simplify(enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO, cutoff = 0.5)
    # nrow(as.data.frame(test))
    # test <- test[order(-test$RichFactor), ]
    # head(test, 25)[, 2:9]
    # test[, 2:9]
    # rm(test)
    
    oraPlot(enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified,
            "AM14 Transfer: PL2-3+2DG vs PL2-3\n(GO Biological Process)", 40,
            "enrichGO_plots/AM14 Transfer - PL2-3 + 2DG vs PL2-3", w = 1250, h = 2200)
    dev.off()  
    
    oraPlot(enrichGO_all$AM14trans$R848_2DG_v_Ctrl$eGO,
            "AM14 Transfer: R848+2DG vs R848\n(GO Biological Process)", 40,
            "enrichGO_plots/AM14 Transfer - R848 + 2DG vs R848", h = 600)
    dev.off()
    
    oraPlot(enrichGO_all$AM14trans$PL23_vs_R848$eGO_simplified,
            "AM14 Transfer: PL2-3 vs R848\n(GO Biological Process)", 40,
            "enrichGO_plots/AM14 Transfer - PL2-3 vs R848")
    dev.off()
    
    oraPlot(enrichGO_all$PL23_v_NP$eGO_simplified,
            "AM14 Transfer + PL2-3\nvs B1-8 Transfer + NP\n(GO Biological Process)", 42,
            "enrichGO_plots/AM14 Transfer PL2-3 vs B1-8 Transfer NP")
    dev.off()
    
    oraPlot(enrichGO_all$AM14MRLlpr$eGO_simplified,
            "AM14 MRL/lpr: 2DG vs Control\n(GO Biological Process)", 40,
            "enrichGO_plots/AM14 MRLlpr 2DG vs Control")
    dev.off()

    # Save results to an Excel file ----------------------------------------------
      wb <- createWorkbook("Output/Functional_analyses/enrich_GOBP_results_allDEGs.xlsx")
      
      addWorksheet(wb, "Enrich GO Summaries")
      writeData(wb, "Enrich GO Summaries", enrichGO_all_summary)
      
      addWorksheet(wb, "AM14 PL23_2DG_vs_Ctrl - Simp")
      addWorksheet(wb, "AM14 - R848_2DG_vs_Ctrl - Simp")
      addWorksheet(wb, "AM14 - PL23_vs_R848 - Simp")
      # addWorksheet(wb, "B1-8 - 2DG_vs_Ctrl - Simp")
      addWorksheet(wb, "PL2-3 vs NP - Simp")
      addWorksheet(wb, "AM14 MRLlpr - Simp")
      
      writeData(wb, "AM14 PL23_2DG_vs_Ctrl - Simp", enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified)
      writeData(wb, "AM14 - R848_2DG_vs_Ctrl - Simp", enrichGO_all$AM14trans$R848_2DG_v_Ctrl$eGO_simplified)
      writeData(wb, "AM14 - PL23_vs_R848 - Simp", enrichGO_all$AM14trans$PL23_vs_R848$eGO_simplified)
      # writeData(wb, "B1-8 - 2DG_vs_Ctrl - Simp", enrichGO_all$B18trans$eGO_simplified)
      writeData(wb, "PL2-3 vs NP - Simp", enrichGO_all$PL23_v_NP$eGO_simplified)
      writeData(wb, "AM14 MRLlpr - Simp", enrichGO_all$AM14MRLlpr$eGO_simplified)
      
      addWorksheet(wb, "AM14 PL23_2DG_vs_Ctrl - Full")
      addWorksheet(wb, "AM14 - R848_2DG_vs_Ctrl - Full")
      addWorksheet(wb, "AM14 - PL23_vs_R848 - Full")
      # addWorksheet(wb, "B1-8 - 2DG_vs_Ctrl - Full")
      addWorksheet(wb, "PL2-3 vs NP - Full")
      addWorksheet(wb, "AM14 MRLlpr - Full")
      
      writeData(wb, "AM14 PL23_2DG_vs_Ctrl - Full", enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO)
      writeData(wb, "AM14 - R848_2DG_vs_Ctrl - Full", enrichGO_all$AM14trans$R848_2DG_v_Ctrl$eGO)
      writeData(wb, "AM14 - PL23_vs_R848 - Full", enrichGO_all$AM14trans$PL23_vs_R848$eGO)
      # writeData(wb, "B1-8 - 2DG_vs_Ctrl - Full", enrichGO_all$B18trans$eGO)
      writeData(wb, "PL2-3 vs NP - Full", enrichGO_all$PL23_v_NP$eGO)
      writeData(wb, "AM14 MRLlpr - Full", enrichGO_all$AM14MRLlpr$eGO)
      
      saveWorkbook(wb, "Output/Functional_analyses/enrich_GOBP_results_allDEGs.xlsx", overwrite = TRUE)
      rm(wb)
    
  # Upregulated DEGs -----------------------------------------------------------
    enrichGO_up <- list(
      AM14trans = list(
        PL23_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "up", 0.5),
        R848_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl, "up", 0.5),
        PL23_vs_R848 = enrichGO_wrapper(deseq_res$AM14transfer$PL23_vs_R848, "up", 0.5)),
      PL23_v_NP = enrichGO_wrapper(deseq_res$PL23_v_NP, "up", 0.5),
      B18trans = enrichGO_wrapper(deseq_res$B18transfer, "up", 0.5),
      AM14MRLlpr = enrichGO_wrapper(deseq_res$AM14MRLlpr, "up", 0.5))
    
    enrichGO_up_summary <- join_all(list(enrichGO_summary(enrichGO_up$AM14trans$PL23_2DG_v_Ctrl, "PL2-3+2DG vs PL2-3", "AM14 Adoptive Transfer", " +0.5"),
                                         enrichGO_summary(enrichGO_up$AM14trans$R848_2DG_v_Ctrl, "R848+2DG vs R848", "AM14 Adoptive Transfer", " +0.5"),
                                         enrichGO_summary(enrichGO_up$AM14trans$PL23_vs_R848, "PL2-3 vs R848", "AM14 Adoptive Transfer", " +0.5"),
                                         enrichGO_summary(enrichGO_up$B18trans, "NP+2DG vs NP", "B1-8 Adoptive Transfer", " +0.5"),
                                         enrichGO_summary(enrichGO_up$PL23_v_NP, "PL2-3 vs NP", "AM14 and B1-8 Adoptive Transfers", " +0.5"),
                                         enrichGO_summary(enrichGO_up$AM14MRLlpr, "2DG vs Control", "AM14 MRL/lpr", " +0.5")), 
                                    type = "full")
    enrichGO_up_summary
  
    
    oraPlot(enrichGO_up$AM14trans$R848_2DG_v_Ctrl$eGO,
            "AM14 Transfer: R848+2DG vs R848\n(Upregulated DEGs,\nGO Biological Process)", 30,
            "enrichGO_plots/AM14 Transfer - R848 + 2DG vs R848 - UP", h = 600)
    dev.off()
    
    oraPlot(enrichGO_up$AM14trans$PL23_vs_R848$eGO,
            "AM14 Transfer: PL2-3 vs R848\n(Upregulated DEGs,\nGO Biological Process)", 30,
            "enrichGO_plots/AM14 Transfer - PL2-3 vs R848 - UP", h = 750)
    dev.off()
    
    oraPlot(enrichGO_up$PL23_v_NP$eGO_simplified,
            "AM14 Transfer + PL2-3 vs\nB1-8 Transfer + NP\n(Upregulated DEGs,\nGO Biological Process)", 42,
            "enrichGO_plots/AM14 Transfer PL2-3 vs B1-8 Transfer NP - UP")
    dev.off()
    
    oraPlot(enrichGO_up$AM14MRLlpr$eGO,
            "AM14 MRL/lpr: 2DG vs Control\n(Upregulated DEGs,\nGO Biological Process)", 40,
            "enrichGO_plots/AM14 MRLlpr 2DG vs Control - UP", , h = 800)
    dev.off()
    
    # Save results to an Excel file ----------------------------------------------
      wb <- createWorkbook("Output/Functional_analyses/enrich_GOBP_results_upregDEGs.xlsx")
      
      addWorksheet(wb, "Enrich GO Summaries")
      writeData(wb, "Enrich GO Summaries", enrichGO_up_summary)
      
      # addWorksheet(wb, "AM14 PL23_2DG_vs_Ctrl - Simp")
      addWorksheet(wb, "AM14 - R848_2DG_vs_Ctrl - Simp")
      addWorksheet(wb, "AM14 - PL23_vs_R848 - Simp")
      # addWorksheet(wb, "B1-8 - 2DG_vs_Ctrl - Simp")
      addWorksheet(wb, "PL2-3 vs NP - Simp")
      addWorksheet(wb, "AM14 MRLlpr - Simp")
      
      # writeData(wb, "AM14 PL23_2DG_vs_Ctrl - Simp", enrichGO_up$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified)
      writeData(wb, "AM14 - R848_2DG_vs_Ctrl - Simp", enrichGO_up$AM14trans$R848_2DG_v_Ctrl$eGO_simplified)
      writeData(wb, "AM14 - PL23_vs_R848 - Simp", enrichGO_up$AM14trans$PL23_vs_R848$eGO_simplified)
      # writeData(wb, "B1-8 - 2DG_vs_Ctrl - Simp", enrichGO_up$B18trans$eGO_simplified)
      writeData(wb, "PL2-3 vs NP - Simp", enrichGO_up$PL23_v_NP$eGO_simplified)
      writeData(wb, "AM14 MRLlpr - Simp", enrichGO_up$AM14MRLlpr$eGO_simplified)
      
      # addWorksheet(wb, "AM14 PL23_2DG_vs_Ctrl - Full")
      addWorksheet(wb, "AM14 - R848_2DG_vs_Ctrl - Full")
      addWorksheet(wb, "AM14 - PL23_vs_R848 - Full")
      # addWorksheet(wb, "B1-8 - 2DG_vs_Ctrl - Full")
      addWorksheet(wb, "PL2-3 vs NP - Full")
      addWorksheet(wb, "AM14 MRLlpr - Full")
      
      # writeData(wb, "AM14 PL23_2DG_vs_Ctrl - Full", enrichGO_up$AM14trans$PL23_2DG_v_Ctrl$eGO)
      writeData(wb, "AM14 - R848_2DG_vs_Ctrl - Full", enrichGO_up$AM14trans$R848_2DG_v_Ctrl$eGO)
      writeData(wb, "AM14 - PL23_vs_R848 - Full", enrichGO_up$AM14trans$PL23_vs_R848$eGO)
      # writeData(wb, "B1-8 - 2DG_vs_Ctrl - Full", enrichGO_up$B18trans$eGO)
      writeData(wb, "PL2-3 vs NP - Full", enrichGO_up$PL23_v_NP$eGO)
      writeData(wb, "AM14 MRLlpr - Full", enrichGO_up$AM14MRLlpr$eGO)
      
      saveWorkbook(wb, "Output/Functional_analyses/enrich_GOBP_results_upregDEGs.xlsx", overwrite = TRUE)
      rm(wb)
      
  # Downregulated DEGs ---------------------------------------------------------
    enrichGO_down <- list(
      AM14trans = list(
        PL23_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "down", 0.5),
        R848_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl, "down", 0.5),
      PL23_vs_R848 = enrichGO_wrapper(deseq_res$AM14transfer$PL23_vs_R848, "down", 0.5)),
      PL23_v_NP = enrichGO_wrapper(deseq_res$PL23_v_NP, "down", 0.5),
      B18trans = enrichGO_wrapper(deseq_res$B18transfer, "down", 0.5),
      AM14MRLlpr = enrichGO_wrapper(deseq_res$AM14MRLlpr, "down", 0.5))
      
    enrichGO_down_summary <- join_all(list(enrichGO_summary(enrichGO_down$AM14trans$PL23_2DG_v_Ctrl, "PL2-3+2DG vs PL2-3", "AM14 Adoptive Transfer", " -0.5"),
                                           enrichGO_summary(enrichGO_down$AM14trans$R848_2DG_v_Ctrl, "R848+2DG vs R848", "AM14 Adoptive Transfer", " -0.5"),
                                           enrichGO_summary(enrichGO_down$AM14trans$PL23_vs_R848, "PL2-3 vs R848", "AM14 Adoptive Transfer", " -0.5"),
                                           enrichGO_summary(enrichGO_down$B18trans, "NP+2DG vs NP", "B1-8 Adoptive Transfer", " -0.5"),
                                           enrichGO_summary(enrichGO_down$PL23_v_NP, "PL2-3 vs NP", "AM14 and B1-8 Adoptive Transfers", " -0.5"),
                                           enrichGO_summary(enrichGO_down$AM14MRLlpr, "2DG vs Control", "AM14 MRL/lpr", " -0.5")), 
                                      type = "full")
    enrichGO_down_summary
    
    
    # enrichGO_down_PL23_3DG_v_Ctrl_lowerSimp <- simplify(enrichGO_down$AM14trans$PL23_2DG_v_Ctrl$eGO, cutoff = 0.6)
    # nrow(enrichGO_down_PL23_3DG_v_Ctrl_lowerSimp)
    # enrichGO_down_PL23_3DG_v_Ctrl_lowerSimp$Description
    # rm(enrichGO_down_PL23_3DG_v_Ctrl_lowerSimp)
    
    oraPlot(enrichGO_down$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified,
            "AM14 Transfer: PL2-3+2DG vs PL2-3\n(Downregulated DEGs,\nGO Biological Process)", 40,
            "enrichGO_plots/AM14 Transfer - PL2-3 + 2DG vs PL2-3 - DOWN")
    dev.off()
    
    oraPlot(enrichGO_down$AM14trans$PL23_vs_R848$eGO_simplified,
            "AM14 Transfer: PL2-3 vs R848\n(Downregulated DEGs,\nGO Biological Process)", 40,
            "enrichGO_plots/AM14 Transfer - PL2-3 vs R848 - DOWN")
    dev.off()
    
    oraPlot(enrichGO_down$PL23_v_NP$eGO,
            "AM14 Transfer + PL2-3 vs B1-8 Transfer + NP\n(Downregulated DEGs,\nGO Biological Process)", 40,
            "enrichGO_plots/AM14 Transfer PL2-3 vs B1-8 Transfer NP - DOWN", h = 800)
    dev.off()
    
    oraPlot(enrichGO_down$AM14MRLlpr$eGO_simplified,
            "AM14 MRL/lpr: 2DG vs Control\n(Downregulated DEGs,\nGO Biological Process)", 40,
            "enrichGO_plots/AM14 MRLlpr 2DG vs Control - DOWN")
    dev.off()
    
    # Save results to an Excel file --------------------------------------------
      wb <- createWorkbook("Output/Functional_analyses/enrich_GOBP_results_downregDEGs.xlsx")
      
      addWorksheet(wb, "Enrich GO Summaries")
      writeData(wb, "Enrich GO Summaries", enrichGO_down_summary)
      
      addWorksheet(wb, "AM14 PL23_2DG_vs_Ctrl - Simp")
      # addWorksheet(wb, "AM14 - R848_2DG_vs_Ctrl - Simp")
      addWorksheet(wb, "AM14 - PL23_vs_R848 - Simp")
      # addWorksheet(wb, "B1-8 - 2DG_vs_Ctrl - Simp")
      addWorksheet(wb, "PL2-3 vs NP - Simp")
      addWorksheet(wb, "AM14 MRLlpr - Simp")
      
      writeData(wb, "AM14 PL23_2DG_vs_Ctrl - Simp", enrichGO_down$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified)
      # writeData(wb, "AM14 - R848_2DG_vs_Ctrl - Simp", enrichGO_down$AM14trans$R848_2DG_v_Ctrl$eGO_simplified)
      writeData(wb, "AM14 - PL23_vs_R848 - Simp", enrichGO_down$AM14trans$PL23_vs_R848$eGO_simplified)
      # writeData(wb, "B1-8 - 2DG_vs_Ctrl - Simp", enrichGO_down$B18trans$eGO_simplified)
      writeData(wb, "PL2-3 vs NP - Simp", enrichGO_down$PL23_v_NP$eGO_simplified)
      writeData(wb, "AM14 MRLlpr - Simp", enrichGO_down$AM14MRLlpr$eGO_simplified)
  
      addWorksheet(wb, "AM14 PL23_2DG_vs_Ctrl - Full")
      # addWorksheet(wb, "AM14 - R848_2DG_vs_Ctrl - Full")
      addWorksheet(wb, "AM14 - PL23_vs_R848 - Full")
      # addWorksheet(wb, "B1-8 - 2DG_vs_Ctrl - Full")
      addWorksheet(wb, "PL2-3 vs NP - Full")
      addWorksheet(wb, "AM14 MRLlpr - Full")
      
      writeData(wb, "AM14 PL23_2DG_vs_Ctrl - Full", enrichGO_down$AM14trans$PL23_2DG_v_Ctrl$eGO)
      # writeData(wb, "AM14 - R848_2DG_vs_Ctrl - Full", enrichGO_down$AM14trans$R848_2DG_v_Ctrl$eGO)
      writeData(wb, "AM14 - PL23_vs_R848 - Full", enrichGO_down$AM14trans$PL23_vs_R848$eGO)
      # writeData(wb, "B1-8 - 2DG_vs_Ctrl - Full", enrichGO_down$B18trans$eGO)
      writeData(wb, "PL2-3 vs NP - Full", enrichGO_down$PL23_v_NP$eGO)
      writeData(wb, "AM14 MRLlpr - Full", enrichGO_down$AM14MRLlpr$eGO)
      
      saveWorkbook(wb, "Output/Functional_analyses/enrich_GOBP_results_downregDEGs.xlsx", overwrite = TRUE)
      rm(wb)
      

# using the packaged "ReactomePA"

# extract EntrezID for Differentially Expressed Genes
de_genes <- sig$ENTREZID
head(de_genes)

# perform pathway enrichment analysis
pathway2 <- enrichPathway(gene=de_genes, pvalueCutoff = 0.05, readable=TRUE)
head(pathway2)
