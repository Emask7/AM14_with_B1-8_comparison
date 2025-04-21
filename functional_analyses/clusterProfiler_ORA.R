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
    
    if(length(sig) < 1) {
      return(NULL)
      print("0 significant DEGs")
    } 
    sig <- as.character(sig$ensembl_gene_id)
      
    eGO <- enrichGO(gene = sig, universe = background, keyType = "ENSEMBL",
                    OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, readable = TRUE)
    
    if(nrow(eGO) < 1) {
      print("ranked gene list length is > 0 but no significant GO terms were detected")
      return(NULL)
    } else {
      simp <- simplify(eGO, cutoff = 0.7)
      return(list(eGO = eGO, eGO_simplified = simp, summary = data.frame(eGO)))
    }
  }

  
  oraPlot <- function(ora, plot_title, file_name, showCat, w, h){
    if(missing(file_name)) file_name <- NULL
    if(missing(showCat)) showCat <- 25
    if(missing(w)) w <- 1700
    if(missing(h)) h <- 2000
    
    if(nrow(ora) < showCat) showCat <- nrow(ora)
    
    if(!is.null(file_name)){
      png(filename = stri_join(c("Output/Functional_analyses/", file_name, ".png"), collapse = ""),
          width = w, height = h, units = "px", pointsize = 10, res = 200,
          bg = "white", family = "", symbolfamily="default")
    }
    ggplot(ora, showCategory = showCat,
           aes(p.adjust, fct_reorder(Description, p.adjust))) +
      geom_segment(aes(xend=0, yend = Description)) +
      geom_point(aes(color=RichFactor, size = Count)) +
      scale_color_gradientn(colors = pnw_palette("Bay", 8, type = "continuous")) +
      scale_size_continuous(range=c(2, 10)) +
      theme_dose(12) +
      xlab("FDR") +
      ylab(NULL) +
      scale_y_discrete(labels = label_wrap(50)) +
      ggtitle(plot_title)
  }
  
  enrichGO_all <- list(
    AM14trans = list(
      PL23_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, NULL),
      R848_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl, NULL),
      PL23_vs_R848 = enrichGO_wrapper(deseq_res$AM14transfer$PL23_vs_R848, NULL)),
    PL23_v_NP = enrichGO_wrapper(deseq_res$PL23_v_NP, NULL),
    B18trans = enrichGO_wrapper(deseq_res$B18transfer, NULL),
    AM14MRLlpr = enrichGO_wrapper(deseq_res$AM14MRLlpr, NULL))

  nrow(enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified)
  nrow(enrichGO_all$AM14trans$R848_2DG_v_Ctrl$eGO_simplified)
  nrow(enrichGO_all$AM14trans$PL23_vs_R848$eGO_simplified)
  nrow(enrichGO_all$PL23_v_NP$eGO_simplified)
  nrow(enrichGO_all$B18trans$eGO_simplified)
  nrow(enrichGO_all$AM14MRLlpr$eGO_simplified)
  
  
  
  oraPlot(enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified, 
          "AM14 Transfer: PL2-3+2DG vs PL2-3\n(GO Biological Process)",
          "enrichGO_plots/AM14 Transfer - PL2-3 + 2DG vs PL2-3", w = 1400, h = 2100)
  dev.off()  
  
  oraPlot(enrichGO_all$AM14trans$R848_2DG_v_Ctrl$eGO_simplified, 
          "AM14 Transfer: R848+2DG vs R848\n(GO Biological Process)",
          "enrichGO_plots/AM14 Transfer - R848 + 2DG vs R848", w = 1500, h = 600)
  dev.off()
  
  oraPlot(enrichGO_all$AM14trans$PL23_vs_R848$eGO_simplified, 
          "AM14 Transfer: PL2-3 vs R848\n(GO Biological Process)",
          "enrichGO_plots/AM14 Transfer - PL2-3 vs R848", w = 1400, h = 2200)
  dev.off()
  
  oraPlot(enrichGO_all$PL23_v_NP$eGO_simplified,
          "AM14 Transfer + PL2-3\nvs B1-8 Transfer + NP\n(GO Biological Process)",
          "enrichGO_plots/AM14 Transfer PL2-3 vs B1-8 Transfer NP", w = 1500, h = 2000)
  dev.off()
  
  oraPlot(enrichGO_all$B18trans$eGO_simplified,
          "B1-8 Transfer: NP + 2DG vs NP\n(GO Biological Process)",
          "enrichGO_plots/B1-8 Transfer NP + 2DG vs NP", w = 1500, h = 800)
  dev.off()
  
  oraPlot(enrichGO_all$AM14MRLlpr$eGO_simplified, 
          "AM14 MRL/lpr: 2DG vs Control\n(GO Biological Process)",
          "enrichGO_plots/AM14 MRLlpr 2DG vs Control", w = 1500, h = 2000)
  dev.off()
  
  

  
  enrichGO_up <- list(
    AM14trans = list(
      PL23_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "up"),
      R848_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl, "up"),
      PL23_vs_R848 = enrichGO_wrapper(deseq_res$AM14transfer$PL23_vs_R848, "up")),
    PL23_v_NP = enrichGO_wrapper(deseq_res$PL23_v_NP, "up"),
    B18trans = enrichGO_wrapper(deseq_res$B18transfer, "up"),
    AM14MRLlpr = enrichGO_wrapper(deseq_res$AM14MRLlpr, "up"))
  
  nrow(enrichGO_up$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified)
  nrow(enrichGO_up$AM14trans$R848_2DG_v_Ctrl$eGO_simplified)
  nrow(enrichGO_up$AM14trans$PL23_vs_R848$eGO_simplified)
  nrow(enrichGO_up$PL23_v_NP$eGO_simplified)
  nrow(enrichGO_up$B18trans$eGO_simplified)
  nrow(enrichGO_up$AM14MRLlpr$eGO_simplified)
  
  # oraPlot(enrichGO_up$PL23_v_NP$eGO_simplified,
  #         "AM14 Transfer + PL2-3 vs B1-8 Transfer + NP\n(Upregulated DEGs, GO Biological Process)",
  #         "enrichGO_plots/AM14 Transfer PL2-3 vs B1-8 Transfer NP - UP")
  # dev.off()
  
  
  enrichGO_down <- list(
    AM14trans = list(
      PL23_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "down"),
      R848_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl, "down"),
      PL23_vs_R848 = enrichGO_wrapper(deseq_res$AM14transfer$PL23_vs_R848, "down")),
    PL23_v_NP = enrichGO_wrapper(deseq_res$PL23_v_NP, "down"),
    B18trans = enrichGO_wrapper(deseq_res$B18transfer, "down"),
    AM14MRLlpr = enrichGO_wrapper(deseq_res$AM14MRLlpr, "down"))

  nrow(enrichGO_down$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified)
  nrow(enrichGO_down$AM14trans$R848_2DG_v_Ctrl$eGO_simplified)
  nrow(enrichGO_down$AM14trans$PL23_vs_R848$eGO_simplified)
  nrow(enrichGO_down$PL23_v_NP$eGO_simplified)
  nrow(enrichGO_down$B18trans$eGO_simplified)
  nrow(enrichGO_down$AM14MRLlpr$eGO_simplified)
  

  oraPlot(enrichGO_down$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified, 
          "AM14 Transfer: PL2-3+2DG vs PL2-3\n(Downregulated DEGs, GO Biological Process)",
          "enrichGO_plots/AM14 Transfer - PL2-3 + 2DG vs PL2-3 - DOWN", w = 1600, h = 1800)
  dev.off()
  
  oraPlot(enrichGO_down$AM14trans$PL23_vs_R848$eGO_simplified, 
          "AM14 Transfer: PL2-3 vs R848\n(Downregulated DEGs,\nGO Biological Process)",
          "enrichGO_plots/AM14 Transfer - PL2-3 vs R848 - DOWN", h = 1200)
  dev.off()
  
  oraPlot(enrichGO_down$AM14MRLlpr$eGO_simplified, 
          "AM14 MRL/lpr: 2DG vs Control\n(Downregulated DEGs,\nGO Biological Process)",
          "enrichGO_plots/AM14 MRLlpr 2DG vs Control - DOWN")
  dev.off()
  



# using the packaged "ReactomePA"

# extract EntrezID for Differentially Expressed Genes
de_genes <- sig$ENTREZID
head(de_genes)

# perform pathway enrichment analysis
pathway2 <- enrichPathway(gene=de_genes, pvalueCutoff = 0.05, readable=TRUE)
head(pathway2)
