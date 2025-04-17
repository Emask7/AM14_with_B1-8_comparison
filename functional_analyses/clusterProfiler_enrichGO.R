# GO Enrichment Analysis -------------------------------------------------------
## background info: http://geneontology.org/docs/ontology-documentation/
  enrichGO_wrapper <- function(dds_res, lfc_direction){
    # Create background dataset for hypergeometric testing using all genes tested for significance in the results
    background <- as.character(dds_res$ensembl_gene_id)
    
    if(missing(lfc_direction)) lfc_direction <- NULL
    
    if(is.null(lfc_direction)) sig <- dplyr::filter(dds_res, padj < 0.05)
    else if(lfc_direction == "up"){
      sig <- dplyr::filter(dds_res, padj < 0.05 & log2FoldChange > 1)
    } else if(lfc_direction == "down"){
      sig <- dplyr::filter(dds_res, padj < 0.05 & log2FoldChange < -1)
    } else return(print("lfc_direction must be either up, down, or NULL"))
    
    sig <- as.character(sig$ensembl_gene_id)
    
    eGO <- enrichGO(gene = sig, universe = background, keyType = "ENSEMBL",
                    OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, readable = TRUE)
    
    return(list(eGO = eGO, summary = data.frame(eGO)))
  }
  
  enrichGO_results <- list(
    all_DEGs = list(
      AM14trans = list(
        PL23_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, NULL),
        R848_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl, NULL),
        PL23_vs_R848 = enrichGO_wrapper(deseq_res$AM14transfer$PL23_vs_R848, NULL)),
      PL23_v_NP = enrichGO_wrapper(deseq_res$PL23_v_NP, NULL),
      B18trans = enrichGO_wrapper(deseq_res$B18transfer, NULL),
      AM14MRLlpr = enrichGO_wrapper(deseq_res$AM14MRLlpr, NULL)), 
    upregulated = list(
      AM14trans = list(
        PL23_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "up"),
        R848_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl, "up"),
        PL23_vs_R848 = enrichGO_wrapper(deseq_res$AM14transfer$PL23_vs_R848, "up")),
      PL23_v_NP = enrichGO_wrapper(deseq_res$PL23_v_NP, "up"),
      B18trans = enrichGO_wrapper(deseq_res$B18transfer, "up"),
      AM14MRLlpr = enrichGO_wrapper(deseq_res$AM14MRLlpr, "up")), 
    downregulated = list(
      AM14trans = list(
        PL23_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "down"),
        R848_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl, "down"),
        PL23_vs_R848 = enrichGO_wrapper(deseq_res$AM14transfer$PL23_vs_R848, "down")),
      PL23_v_NP = enrichGO_wrapper(deseq_res$PL23_v_NP, "down"),
      B18trans = enrichGO_wrapper(deseq_res$B18transfer, "down"),
      AM14MRLlpr = enrichGO_wrapper(deseq_res$AM14MRLlpr, "down"))
  )
  
  # simp_test <- simplify(enrichGO_results$all_DEGs$PL23_v_NP$eGO, cutoff = 0.7)
  # nrow(simp_test)
  # nrow(enrichGO_results$all_DEGs$PL23_v_NP$eGO)
  # 
  # data.frame(simp_test)
  
  head(enrichGO_results$downregulated$AM14trans$PL23_2DG_v_Ctrl$summary)
  
  wb <- createWorkbook("Output/Functional_analyses/GO_enrichment_allDEGs.xlsx")
  addWorksheet(wb, "AM14 - PL23_2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 - R848_2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 - PL23_vs_R848")
  addWorksheet(wb, "B1-8 - 2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 PL2-3 vs B1-8 NP")
  addWorksheet(wb, "AM14 MRLlpr - 2DG_vs_Ctrl")
  writeData(wb, "AM14 - PL23_2DG_vs_Ctrl", enrichGO_results$all_DEGs$AM14trans$PL23_2DG_v_Ctrl$summary)
  writeData(wb, "AM14 - R848_2DG_vs_Ctrl", enrichGO_results$all_DEGs$AM14trans$R848_2DG_v_Ctrl$summary)
  writeData(wb, "AM14 - PL23_vs_R848", enrichGO_results$all_DEGs$AM14trans$PL23_vs_R848$summary)
  writeData(wb, "B1-8 - 2DG_vs_Ctrl", enrichGO_results$all_DEGs$B18trans$summary)
  writeData(wb, "AM14 PL2-3 vs B1-8 NP", enrichGO_results$all_DEGs$PL23_v_NP$summary)
  writeData(wb, "AM14 MRLlpr - 2DG_vs_Ctrl", enrichGO_results$all_DEGs$AM14MRLlpr$summary)
  saveWorkbook(wb, "Output/Functional_analyses/GO_enrichment_allDEGs.xlsx", overwrite = TRUE)
  rm(wb)
  
  wb <- createWorkbook("Output/Functional_analyses/GO_enrichment_Upregulated.xlsx")
  addWorksheet(wb, "AM14 - PL23_2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 - R848_2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 - PL23_vs_R848")
  addWorksheet(wb, "B1-8 - 2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 PL2-3 vs B1-8 NP")
  addWorksheet(wb, "AM14 MRLlpr - 2DG_vs_Ctrl")
  writeData(wb, "AM14 - PL23_2DG_vs_Ctrl", enrichGO_results$upregulated$AM14trans$PL23_2DG_v_Ctrl$summary)
  writeData(wb, "AM14 - R848_2DG_vs_Ctrl", enrichGO_results$upregulated$AM14trans$R848_2DG_v_Ctrl$summary)
  writeData(wb, "AM14 - PL23_vs_R848", enrichGO_results$upregulated$AM14trans$PL23_vs_R848$summary)
  writeData(wb, "B1-8 - 2DG_vs_Ctrl", enrichGO_results$upregulated$B18trans$summary)
  writeData(wb, "AM14 PL2-3 vs B1-8 NP", enrichGO_results$upregulated$PL23_v_NP$summary)
  writeData(wb, "AM14 MRLlpr - 2DG_vs_Ctrl", enrichGO_results$upregulated$AM14MRLlpr$summary)
  saveWorkbook(wb, "Output/Functional_analyses/GO_enrichment_Upregulated.xlsx", overwrite = TRUE)
  rm(wb)
  
  wb <- createWorkbook("Output/Functional_analyses/GO_enrichment_Downregulated.xlsx")
  addWorksheet(wb, "AM14 - PL23_2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 - R848_2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 - PL23_vs_R848")
  addWorksheet(wb, "B1-8 - 2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 PL2-3 vs B1-8 NP")
  addWorksheet(wb, "AM14 MRLlpr - 2DG_vs_Ctrl")
  writeData(wb, "AM14 - PL23_2DG_vs_Ctrl", enrichGO_results$downregulated$AM14trans$PL23_2DG_v_Ctrl$summary)
  writeData(wb, "AM14 - R848_2DG_vs_Ctrl", enrichGO_results$downregulated$AM14trans$R848_2DG_v_Ctrl$summary)
  writeData(wb, "AM14 - PL23_vs_R848", enrichGO_results$downregulated$AM14trans$PL23_vs_R848$summary)
  writeData(wb, "B1-8 - 2DG_vs_Ctrl", enrichGO_results$downregulated$B18trans$summary)
  writeData(wb, "AM14 PL2-3 vs B1-8 NP", enrichGO_results$downregulated$PL23_v_NP$summary)
  writeData(wb, "AM14 MRLlpr - 2DG_vs_Ctrl", enrichGO_results$downregulated$AM14MRLlpr$summary)
  saveWorkbook(wb, "Output/Functional_analyses/GO_enrichment_Downregulated.xlsx", overwrite = TRUE)
  rm(wb)
  
  

### Visualizing clusterProfiler results
# plot(barplot(enrichGO_results$AM14trans$PL23_2DG_v_Ctrl$eGO, showCateGOry=25))
# goplot(enrichGO_results$AM14trans$PL23_2DG_v_Ctrl$eGO)

# Dotplot Wrapper Function -----------------------------------------------------
  dotplot_wrapper <- function(eGO_res, n, plot_title, file_name, h){
    if(missing(n)) n <- 25
    if(missing(h)) h <- 2200
    if(nrow(eGO_res) < n) n <- nrow(eGO_res)
    
    eGO_top <- eGO_res[order(eGO_res$p.adjust, decreasing = FALSE), ]
    eGO_top <- eGO_top[1:n, ]$Description
    
    if(!missing(file_name)){
      if(file_name != FALSE & !is.null(file_name)){
        png(filename = stri_join(c("Output/Functional_analyses/", file_name, ".png"),
                                 collapse = ""),
            width = 1200, height = h, units = "px", pointsize = 10, res = 200,
            bg = "white", family = "", symbolfamily="default")
        dotplot(eGO_res, showCategory = eGO_top, title = plot_title)
      } else dotplot(eGO_res, showCategory = eGO_top, title = plot_title)
    } else dotplot(eGO_res, showCategory = eGO_top, title = plot_title)
    # NOTE: GeneRatio = # genes in your input list associated with the given GO term/total # of input genes.
  }

  # All DEGs
    dotplot_wrapper(enrichGO_results$AM14trans$PL23_2DG_v_Ctrl$eGO, 25, 
                    "AM14 Adoptive Transfer:\nPL2-3 + 2DG vs PL2-3",
                    "enrichGO_plots/AM14 Transfer - PL2-3 + 2DG vs PL2-3", 1600)
    dev.off()
    
    dotplot_wrapper(enrichGO_results$AM14trans$PL23_vs_R848$eGO, 25, 
                    "AM14 Adoptive Transfer:\nPL2-3 vs R848",
                    "enrichGO_plots/AM14 Transfer - PL2-3 vs R848")
    dev.off()
    
    dotplot_wrapper(enrichGO_results$PL23_v_NP$eGO, 25, 
                    "AM14 Transfer + PL2-3 vs\nB1-8 Transfer + NP",
                    "enrichGO_plots/AM14 Transfer PL2-3 vs B1-8 Transfer NP")
    dev.off()
    
    dotplot_wrapper(enrichGO_results$AM14MRLlpr$eGO, 25, 
                    "AM14 MRLlpr:\n2DG vs Control",
                    "enrichGO_plots/AM14 MRLlpr")
    dev.off()

  # Upregulated DEGs
    # Pl2-3 vs NP
    dotplot_wrapper(enrichGO_results$upregulated$PL23_v_NP$eGO, 25, 
                    "AM14 Transfer + PL2-3 vs B1-8 Transfer\n+ NP\n(Upregulated)",
                    "enrichGO_plots/AM14 Transfer PL2-3 vs B1-8 Transfer NP - Upregulated", 2500)
    dev.off()
  
  # Downregulated DEGs
    # AM14 Pl23 +/- 2DG, Pl23 v R848, AM14 MRLlpr
    

# ## Enrichmap clusters the 25 most significant (by padj) GO terms to visualize relationships between terms
# eGO2 <- pairwise_termsim(enrichGO_results$AM14trans$PL23_2DG_v_Ctrl$eGO) 
# emapplot(eGO2)

