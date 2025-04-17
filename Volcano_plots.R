vp_data = list(
  AM14transfer = list(
    PL23_2DG_v_Ctrl = deseq_res$AM14transfer$PL23_2DG_v_Ctrl[, c(2, 6, 10)],
    R848_2DG_v_Ctrl = deseq_res$AM14transfer$R848_2DG_v_Ctrl[, c(2, 6, 10)],
    PL23_vs_R848 = deseq_res$AM14transfer$PL23_vs_R848[, c(2, 6, 10)]),
  B18transfer = deseq_res$B18transfer[, c(2, 6, 10)],
  PL23vNP = deseq_res$PL23_v_NP[, c(2, 6, 10)],
  AM14MRLlpr = deseq_res$AM14MRLlpr[, c(2, 6, 10)])
head(vp_data$AM14MRLlpr)

volcano_wrapper <- function(dat, plot_title, cap, x_limits, y_limits, file_name){
  if(missing(x_limits)) {
    x_min <- min(dat$log2FoldChange)
    x_max <- max(dat$log2FoldChange)
    x_limits <- c(floor(x_min), ceiling(x_max))
    print(stri_join("x_min = ", x_min))
    print(stri_join("x_max = ", x_max))
  }
  if(missing(y_limits)) {
    y_max <- -log10(min(dat$padj))
    y_limits <- c(0, ceiling(y_max))
    print(stri_join("y_max = ", y_max))
  }
  
  vp <- EnhancedVolcano(dat, lab = dat$external_gene_name, 
                        pCutoff = 0.05, FCcutoff = 1,
                        x = 'log2FoldChange', y = 'padj', 
                        xlim = x_limits, ylim = y_limits, 
                        title = plot_title, subtitle = "(Adjusted p-values)", 
                        caption = cap,
                        legendLabels = c("NS", expression(Log[2] ~ FC > 1),
                                         expression(p - value < ~ 0.05),
                                         expression(p - value < ~ 0.05 ~ and ~ log[2] ~ FC > ~ 1)),
                        drawConnectors = TRUE, min.segment.length = 1, 
                        max.overlaps = 12, labSize = 4)
  
  if(!missing(file_name)){
    if(file_name != FALSE & !is.null(file_name)){
      ggsave(stri_join("Output/Volcano_plots/Volcano plot - ", file_name, ".png"),
             units = "px", width = 3000, height = 3000, dpi = 300)
    } 
  } 
  vp
}  

confirm_foldchange <- function(dds, gene, gene_ID_table){
  ens_ID <- subset(gene_ID_table, external_gene_name == gene)
  ens_ID <- ens_ID$ensembl_gene_id
  cts <- counts(dds, normalized = TRUE)
  res <- cts[ens_ID, ]
  res
}
confirm_foldchange(dds_AM14trans, "Crisp1", gene_IDs$AM14trans)
confirm_foldchange(dds_B18trans, "Dmrt2", gene_IDs$B18trans)
confirm_foldchange(dds_PL23vNP, "Lars2", gene_IDs$AM14trans)


volcano_wrapper(vp_data$AM14transfer$PL23_2DG_v_Ctrl, "AM14 Transfer: PL2-3 + 2DG vs PL2-3",
                "Upregulated in PL2-3                                                                                  Upregulated in PL2-3 + 2DG",
                c(-20, 20), c(0, 8), "AM14_Transfer-PL2-3_2DG_vs_Control")

volcano_wrapper(vp_data$AM14transfer$R848_2DG_v_Ctrl, "AM14 Transfer: R848 + 2DG vs R848", 
                "Upregulated in R848                                                                                  Upregulated in R848 + 2DG",
                c(-27, 27), c(0, 10), "AM14_Transfer-R848_2DG_vs_Control")

volcano_wrapper(vp_data$AM14transfer$PL23_vs_R848, "AM14 Transfer: PL2-3 vs R848", 
                "Upregulated in R848                                                                                            Upregulated in PL2-3",
                c(-23, 5), c(0, 9), "AM14_Transfer-PL2-3_vs_R848")

volcano_wrapper(vp_data$B18transfer, "B1-8 Transfer: NP + 2DG vs NP", 
                "Upregulated in NP                                                                                            Upregulated in NP + 2DG",
                c(-5, 3), c(0, 8), "B1-8_Transfer-NP_2DG_vs_NP")

volcano_wrapper(vp_data$AM14MRLlpr, "AM14 MRL/lpr: 2DG vs Control", 
                "Upregulated in Control                                                                                         Upregulated in 2DG",
                c(-5, 2.5), c(0, 15), "AM14_MRLlpr")

volcano_wrapper(vp_data$PL23vNP, "AM14 Transfer PL2-3 vs B1-8 Transfer NP", 
                "Upregulated in NP                                                                                            Upregulated in PL2-3",
                c(-6.5, 12), c(0, 182), "PL23_vs_NP")

