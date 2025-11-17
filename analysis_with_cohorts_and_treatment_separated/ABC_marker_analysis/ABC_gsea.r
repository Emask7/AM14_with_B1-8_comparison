setwd("./analysis_with_cohorts_and_treatment_separated/ABC_marker_analysis")
getwd()

# gmt_sets <- read.gmt("Bcell_subset.human.jack.gmt")
gmt_sets <- read.gmt("Bcell_subset.human.jack - Copy.gmt")
gmt_sets <- gmt_sets[gmt_sets$gene != "", ]
head(gmt_sets)
# gmt_sets
nrow(gmt_sets)


# Filter out less relevant gene sets
  excluded_terms = c("BIOCARTA_STRESS_PATHWAY", "ABC signature in Malaria", "ABC signature in HIV",
                     "GOBP_CELL_CYCLE_G1_S_PHASE_TRANSITION", "GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION",
                     "GOBP_LEUKOCYTE_CELL_CELL_ADHESION")                
  gmt_sets <- gmt_sets[!(gmt_sets$term %in% excluded_terms), ]
  
  
  gmt_sets <- gmt_sets[gmt_sets$term != "BIOCARTA_STRESS_PATHWAY", ]
  gmt_sets <- gmt_sets[gmt_sets$term != "Cell cycle G1/S phase transition (MSigDB: GOBP_CELL_CYCLE_G1_S_PHASE_TRANSITION)", ]
  gmt_sets <- gmt_sets[gmt_sets$term != "Cell cycle G2/M phase transition (MSigDB: GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION)", ]
  # gmt_sets <- gmt_sets[gmt_sets$term != "Myc targets V1 (MSigDB: HALLMARK_MYC_TARGETS_V1)", ]
  # gmt_sets <- gmt_sets[gmt_sets$term != "Myc targets V2 (MSigDB: HALLMARK_MYC_TARGETS_V2)", ]
  gmt_sets <- gmt_sets[gmt_sets$term != "ABC signature in HIV", ]
  gmt_sets <- gmt_sets[gmt_sets$term != "ABC signature in Malaria", ]
  gmt_sets <- gmt_sets[gmt_sets$term != "BIOCARTA_STRESS_PATHWAY", ]
  # gmt_sets <- gmt_sets[gmt_sets$term != "GOBP_RESPONSE_TO_INTERFERON_GAMMA", ]
  # gmt_sets <- gmt_sets[gmt_sets$term != "GOBP_RESPONSE_TO_TYPE_I_INTERFERON", ]
  # gmt_sets <- gmt_sets[gmt_sets$term != "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION", ]
  gmt_sets <- gmt_sets[gmt_sets$term != "GOBP_LEUKOCYTE_CELL_CELL_ADHESION", ]
  # gmt_sets <- gmt_sets[gmt_sets$term != "", ]
  # gmt_sets <- gmt_sets[gmt_sets$term != "", ]
  
  # gmt_sets <- gmt_sets[gmt_sets$term == "ABC signature in Lupus (ref: Wu et al.)", ]


gmt_sets_names <- gmt_sets$term
gmt_sets_names <- gmt_sets_names[!duplicated(gmt_sets_names)]
gmt_sets_names

gsea_format <- function(dds_res){
  res_filtered <- dds_res[dds_res$pvalue <= 0.05, ]
  res_filtered <- res_filtered[res_filtered$hgnc_symbol != "", ]
  genelist <- res_filtered$log2FoldChange
  names(genelist) <- res_filtered$hgnc_symbol
  genelist <- sort(genelist, decreasing = TRUE)
  genelist <- genelist[!duplicated(names(genelist))]
  return(genelist)
}

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


head(deseq_res$AM14transfer_PL23)

pl23_genelist <- gsea_format(deseq_res$AM14transfer_PL23)
MRLlpr_genelist <- gsea_format(deseq_res$AM14MRLlpr)
r848_enelist <- gsea_format(deseq_res$AM14transfer_R848)
np_genelist <- gsea_format(deseq_res$B18transfer)



# GSEA of GMT gene sets --------------------------------------------------------
  pl23_ABC_gsea <- GSEA(pl23_genelist, TERM2GENE = gmt_sets, pvalueCutoff = 1)
  MRLlpr_ABC_gsea <- GSEA(MRLlpr_genelist, TERM2GENE = gmt_sets, pvalueCutoff = 1)
  r848_ABC_gsea <- GSEA(r848_enelist, TERM2GENE = gmt_sets, pvalueCutoff = 1)
  np_ABC_gsea <- GSEA(np_genelist, TERM2GENE = gmt_sets, pvalueCutoff = 1)
  
  head(pl23_ABC_gsea)
  nrow(pl23_ABC_gsea)
  pl23_ABC_gsea[, 1]
  pl23_ABC_gsea[, 5:7]
  pl23_ABC_gsea[8, ]

  gsea_dotplot(pl23_ABC_gsea, plottitle = "PL2-3 + 2DG vs PL2-3")
  ggsave(filename = "ABC markers - PL2-3.png", width = 6.5, height = 6, units = "in", dpi = 600)
  
  gsea_dotplot(MRLlpr_ABC_gsea, plottitle = "AM14 MRL/lpr: 2DG vs Control")
  ggsave(filename = "ABC markers - AM14 MRLlpr.png", width = 6.5, height = 6, units = "in", dpi = 600)
  
  gsea_dotplot(r848_ABC_gsea, plottitle = "R848 + 2DG vs R848")
  ggsave(filename = "ABC markers - R848.png", width = 6.5, height = 6, units = "in", dpi = 600)
  
  gsea_dotplot(np_ABC_gsea, plottitle = "NP + 2DG vs NP")
  ggsave(filename = "ABC markers - NP.png", width = 6.5, height = 4, units = "in", dpi = 600)
  
  
# GSEA_Histograms --------------------------------------------------------------
  gseaplot2(pl23_ABC_gsea, geneSetID = 1, title = pl23_ABC_gsea$Description[1])
  ggsave(filename = "GSEA_Histograms/Myc Targets v1 - PL2-3.png", width = 5, height = 4, units = "in", dpi = 600)
  
  gseaplot2(pl23_ABC_gsea, geneSetID = 2, title = pl23_ABC_gsea$Description[2])
  ggsave(filename = "GSEA_Histograms/Response to IFNG - PL2-3.png", width = 5, height = 4, units = "in", dpi = 600)
  
  gseaplot2(pl23_ABC_gsea, geneSetID = 3, title = pl23_ABC_gsea$Description[3])
  ggsave(filename = "GSEA_Histograms/GC vs Naive B Cell Up - PL2-3.png", width = 5, height = 4, units = "in", dpi = 600)
  
  gseaplot2(pl23_ABC_gsea, geneSetID = 4, title = pl23_ABC_gsea$Description[4])
  ggsave(filename = "GSEA_Histograms/Ag processing and presentation - PL2-3.png", width = 5, height = 4, units = "in", dpi = 600)
  
  gseaplot2(pl23_ABC_gsea, geneSetID = 5, title = pl23_ABC_gsea$Description[5])
  ggsave(filename = "GSEA_Histograms/GC B cell vs Plasma Cell Up - PL2-3.png", width = 5, height = 4, units = "in", dpi = 600)
  
  
  gseaplot2(MRLlpr_ABC_gsea, geneSetID = 1, title = stri_join("AM14 MRL/lpr: 2DG vs Control", MRLlpr_ABC_gsea$Description[1], sep = "\n"))
  ggsave(filename = "GSEA_Histograms/Response to IFNG - AM14 MRLlpr.png", width = 5, height = 4, units = "in", dpi = 600)
  
  
  
  # p1 <- gseaplot2(pl23_ABC_gsea, geneSetID = 1, title = pl23_ABC_gsea$Description[1])
  # p2 <- gseaplot(pl23_ABC_gsea, geneSetID = 2, title = pl23_ABC_gsea$Description[2])
  # p3 <- gseaplot(pl23_ABC_gsea, geneSetID = 3, title = pl23_ABC_gsea$Description[4])
  # p4 <- gseaplot(pl23_ABC_gsea, geneSetID = 8, title = pl23_ABC_gsea$Description[10])
  # cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])
  # cowplot::plot_grid(p1, p3, p4)
  # 
  # p1 
  #   # facet_wrap( ~ ID, ncol = 2) +
    
  
  
  
  # ggplot(pl23_ABC_gsea, 
  #        aes(p.adjust, fct_reorder(Description, p.adjust))) +
  #   geom_segment(aes(xend=0, yend = Description)) +
  #   geom_point(aes(color=NES, size = setSize)) +
  #   scale_size_continuous(range=c(2, 10)) +
  #   theme_dose(12) +
  #   xlab("FDR") +
  #   ylab(NULL) +
  #   scale_y_discrete(labels = label_wrap(50)) #+
  #   # scale_y_discrete(labels = lapply(strwrap(gse$Description, width = 45, simplify = FALSE), paste, collapse="\n")) +
  #   # ggtitle(plot_title) +
  #   # scale_color_gradientn(colors = pnw_palette("Bay", 8, type = "continuous"), limits = c(-1*nes_max, nes_max))
  
  
  
  
  # gsePlot <- function(gse, plot_title, file_name, showCat, w, h){
  #   if(missing(file_name)) file_name <- NULL
  #   if(missing(showCat)) showCat <- 25
  #   if(missing(w)) w <- 1700
  #   if(missing(h)) h <- 2000
  #   
  #   if(nrow(gse) < showCat) showCat <- nrow(gse)
  #   
  #   gse_top <- gse[1:showCat, ]
  #   nes_max <- max(abs(gse_top$NES))
  #   
  #   if(!is.null(file_name)){
  #     png(filename = stri_join(c("Output/Functional_analyses/", file_name, ".png"), collapse = ""),
  #         width = w, height = h, units = "px", pointsize = 10, res = 200,
  #         bg = "white", family = "", symbolfamily="default")
  #   }
  #   ggplot(gse, showCategory = showCat, 
  #          aes(p.adjust, fct_reorder(Description, p.adjust))) +
  #     geom_segment(aes(xend=0, yend = Description)) +
  #     geom_point(aes(color=NES, size = setSize)) +
  #     scale_size_continuous(range=c(2, 10)) +
  #     theme_dose(12) +
  #     xlab("FDR") +
  #     ylab(NULL) +
  #     scale_y_discrete(labels = label_wrap(50)) +
  #     # scale_y_discrete(labels = lapply(strwrap(gse$Description, width = 45, simplify = FALSE), paste, collapse="\n")) +
  #     ggtitle(plot_title) +
  #     scale_color_gradientn(colors = pnw_palette("Bay", 8, type = "continuous"), limits = c(-1*nes_max, nes_max))
  #   # scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous"), limits = c(-1*nes_max, nes_max))
  #   # scale_color_gradientn(colors = met.brewer("Hiroshige"), limits = c(-1*nes_max, nes_max))
  # }
  
  
  
