# Venn Diagrams ----------------------------------------------------------------
  deg_list_for_venn <- function(res, lfc_cutoff, padj_cutoff){
    if(missing(lfc_cutoff)) lfc_cutoff <- 1
    if(missing(padj_cutoff)) padj_cutoff <- 0.05
    
    list(up = subset(res, res$padj <= padj_cutoff & res$log2FoldChange >= lfc_cutoff)[, 1], 
         down = subset(res, res$padj <= padj_cutoff & res$log2FoldChange <= (-1*lfc_cutoff))[, 1], 
         all = subset(res, res$padj <= padj_cutoff & abs(res$log2FoldChange) >= lfc_cutoff)[, 1])
  }

# Volcano Plots ----------------------------------------------------------------
  volcano_wrapper <- function(dat, plot_title, cap, lfc_cutoff, padj_cutoff, file_name, h, w, x_limits, y_limits){
    if(missing(plot_title)) plot_title <- "Plot Title"
    if(missing(cap)) cap <- " "
    if(missing(lfc_cutoff)) lfc_cutoff <- 1
    if(missing(padj_cutoff)) padj_cutoff <- 0.05
    if(missing(file_name)) file_name <- NULL
    if(missing(h)) h <- 2500
    if(missing(w)) w <- 2500
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
                          FCcutoff = lfc_cutoff, pCutoff = padj_cutoff,
                          x = 'log2FoldChange', y = 'padj', 
                          xlim = x_limits, ylim = y_limits, 
                          title = plot_title, subtitle = "(Adjusted p-values)", 
                          caption = cap,
                          # col = wes_palette("Zissou1", type = "continuous"),
                          drawConnectors = TRUE, min.segment.length = 1, 
                          max.overlaps = 12, labSize = 4)
    
    if(!is.null(file_name)){
        if(!file.exists("Output/")) dir.create("Output/")
        if(!file.exists("Output/Volcano_plots/")) dir.create("Output/Volcano_plots/")
      
        ggsave(stri_join("Output/Volcano_plots/", file_name, ".png"),
               units = "px", width = w, height = h, dpi = 300)
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
  