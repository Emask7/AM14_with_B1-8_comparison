setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/ABC_gsea/")
getwd()

# Set up genesets from GMT file ------------------------------------------------
  # gmt_sets <- read.gmt("Bcell_subset.human.jack.gmt")
  gmt_sets <- read.gmt("Bcell_subset.human.jack - Copy.gmt")
  gmt_sets <- gmt_sets[gmt_sets$gene != "", ]
  head(gmt_sets)
  # gmt_sets
  nrow(gmt_sets)
  
  # Filter out less relevant gene sets------------------------------------------
    excluded_terms = c("BIOCARTA_STRESS_PATHWAY", "ABC signature in Malaria", "ABC signature in HIV",
                       "GOBP_CELL_CYCLE_G1_S_PHASE_TRANSITION", "GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION",
                       "GOBP_LEUKOCYTE_CELL_CELL_ADHESION")                
    gmt_sets <- gmt_sets[!(gmt_sets$term %in% excluded_terms), ]
    
    
    # gmt_sets <- gmt_sets[gmt_sets$term != "BIOCARTA_STRESS_PATHWAY", ]
    # gmt_sets <- gmt_sets[gmt_sets$term != "Cell cycle G1/S phase transition (MSigDB: GOBP_CELL_CYCLE_G1_S_PHASE_TRANSITION)", ]
    # gmt_sets <- gmt_sets[gmt_sets$term != "Cell cycle G2/M phase transition (MSigDB: GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION)", ]
    # gmt_sets <- gmt_sets[gmt_sets$term != "Myc targets V1 (MSigDB: HALLMARK_MYC_TARGETS_V1)", ]
    # gmt_sets <- gmt_sets[gmt_sets$term != "Myc targets V2 (MSigDB: HALLMARK_MYC_TARGETS_V2)", ]
    # gmt_sets <- gmt_sets[gmt_sets$term != "ABC signature in HIV", ]
    # gmt_sets <- gmt_sets[gmt_sets$term != "ABC signature in Malaria", ]
    # gmt_sets <- gmt_sets[gmt_sets$term != "BIOCARTA_STRESS_PATHWAY", ]
    # gmt_sets <- gmt_sets[gmt_sets$term != "GOBP_RESPONSE_TO_INTERFERON_GAMMA", ]
    # gmt_sets <- gmt_sets[gmt_sets$term != "GOBP_RESPONSE_TO_TYPE_I_INTERFERON", ]
    # gmt_sets <- gmt_sets[gmt_sets$term != "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION", ]
    # gmt_sets <- gmt_sets[gmt_sets$term != "GOBP_LEUKOCYTE_CELL_CELL_ADHESION", ]
    # gmt_sets <- gmt_sets[gmt_sets$term != "", ]
    # gmt_sets <- gmt_sets[gmt_sets$term != "", ]
    
    # gmt_sets <- gmt_sets[gmt_sets$term == "ABC signature in Lupus (ref: Wu et al.)", ]

    # gmt_sets_names <- gmt_sets$term
    # gmt_sets_names <- gmt_sets_names[!duplicated(gmt_sets_names)]
    # gmt_sets_names

# Input DESeq2 results and filter LFC and pvalue data for genes in sets --------
  gsea_format <- function(dds_res){
    res_filtered <- dds_res[dds_res$pvalue <= 0.05, ]
    res_filtered <- res_filtered[res_filtered$hgnc_symbol != "", ]
    genelist <- res_filtered$log2FoldChange
    names(genelist) <- res_filtered$hgnc_symbol
    genelist <- sort(genelist, decreasing = TRUE)
    genelist <- genelist[!duplicated(names(genelist))]
    return(genelist)
  }
  
# Create dotplot from ABC GSEA results -----------------------------------------
  gsea_dotplot <- function(gsea_res, plottitle, showCat, w, h){
    if(missing(plottitle)) plottitle <- ""
    if(missing(showCat)) showCat <- 30
    if(missing(w)) w <- 1700
    if(missing(h)) h <- 2000
    
    ggplot(data = gsea_res, showCategory = showCat,
           # aes(x = pvalue, y = fct_reorder(Description, pvalue))) +
           aes(x = p.adjust, y = fct_reorder(Description, p.adjust))) +
      geom_point(aes(size = setSize, color = NES)) + 
      theme_bw() +
      theme(panel.grid.minor=element_blank(),
            axis.text = element_text(color = "black"),
            legend.text = element_text(size = 12)) +
      # scale_colour_gradient2(low = "blue", high = "red") +
      # scale_colour_gradientn(colors = wes_palette("Zissou1")) +
      scale_colour_gradient2(low = "#0f85a0", mid = "white", high = "#ed8b00", 
                             midpoint = 0, na.value = "grey80", name = "Z-score",
                             guide = guide_colorbar(order = 1)) +
      ylab("") + 
      # scale_x_continuous(name = "p.adjust") +
      geom_vline(xintercept = 0.05, linetype = "dashed", color = "black", size = 0.5) +
      # scale_y_discrete(labels = lapply(strwrap(pl23_ABC_gsea$Description, width = 70, simplify = FALSE), paste, collapse="\n")) +
      labs(title = plottitle)
  }
