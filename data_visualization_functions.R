setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/")
getwd()

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
  
# Heatmaps ---------------------------------------------------------------------
  zscore_matrix <- function(dds){
    zscores <- counts(dds, normalized=TRUE)
    res_colnames <- colnames(zscores)
    zscores <- t(zscores)
    zscores <- scale(zscores)
    zscores <- t(zscores)
    zscores <- data.frame(rownames(zscores), zscores)
    zscores <- zscores[rownames(zscores) %in% DEGs, ]
    colnames(zscores) <- c("ensembl_gene_id", res_colnames)
    rm(res_colnames)
    return(zscores)
  }
  
# DESeq2 result bar plots ------------------------------------------------------
  barplot_counts <- function(gene_list, species_source, include_all_comparisons, alternate_IDs){
    if(missing(species_source)) species_source <- "mouse"
    if(missing(include_all_comparisons)) include_all_comparisons <- TRUE
    if(missing(alternate_IDs)) alternate_IDs <- NULL
    
    if(species_source == "human"){
      gene_list <- gene_IDs_full[gene_IDs_full$hgnc_symbol %in% gene_list, ]
      gene_list <- gene_list$external_gene_name
    }
    
    if(include_all_comparisons == TRUE){
      res <- four_way_barplot[four_way_barplot$ID %in% gene_list, ]
      res$comparison <- factor(res$comparison, levels = c("PL2-3", "MRL/lpr", "R848", "NP"))
    } else {
      res <- two_way_barplot[two_way_barplot$ID %in% gene_list, ]
      res$comparison <- factor(res$comparison, levels = c("PL2-3", "MRL/lpr"))
    }
    
    if(!is.null(alternate_IDs)){
      gene_IDs <- data.frame(ID = gene_list)
      gene_IDs <- full_join(gene_IDs, alternate_IDs)
      
      for (r in 1:nrow(gene_IDs)) {
        if(is.na(gene_IDs[r, 2])) gene_IDs[r, 2] = gene_IDs[r, 1]
      }
      
      res <- full_join(gene_IDs, res, by = "ID")
      res <- res[, c(2:6)]
      colnames(res) <- c("ID", colnames(res)[2:5])
    }
    
    return(res)
  }
  
  deseq2_expr_barplot <- function(dat, plot_title, lfc_limits, na_color, color_breaks, color_labels){
    if(missing(plot_title)) plot_title <- ""
    if(missing(lfc_limits)){
      print(stri_join("LFC min:  ", min(dat$LFC)))
      print(stri_join("LFC max:  ", max(dat$LFC)))
      max_abs_lfc <- max(abs(dat$LFC))
      lfc_limits <- c(-1*max_abs_lfc, max_abs_lfc)
    }
    if(missing(na_color)) na_color <- "gray50"
    if(missing(color_breaks)) color_breaks <- waiver()
    if(missing(color_labels)) color_labels <- waiver()
    
    ymax <- max(dat$`-log.pval`) + 5
    
    res <- ggplot(data = dat, 
                  # aes(x = factor(ID), y = `-log.pval`, fill = LFC)) +
                  aes(x = factor(reorder(ID, LFC)), y = `-log.pval`, fill = LFC)) +
      geom_col() +
      theme_bw() +
      theme(panel.grid.minor=element_blank(),
            axis.text = element_text(color = "black", size = 12),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 12)) +
      geom_text(aes(label = ifelse(abs(LFC) > max(abs(lfc_limits)), 
                                   stri_join("LFC: ", round(LFC, digits = 1)), "")), 
                size = 10 / .pt, hjust = -0.05, color = "gray20") +
      # scale_fill_gradientn(colors = wes_palette("Zissou1"), 
      scale_fill_gradientn(colors = pnw_palette("Bay"), 
                           limits = lfc_limits, na.value = na_color,
                           breaks = color_breaks, labels = color_labels) +
      ylab(bquote(~-Log[10] ~ "(p-value)")) +
      scale_y_continuous(limits = c(0, ymax)) +
      xlab("") +
      coord_flip() +
      geom_hline(yintercept = 1.3, linetype = "dashed", color = "black", size = 0.5) +
      labs(title = plot_title,
           subtitle = "(DESeq2 Differential Expression Results)") +
      facet_wrap( ~ comparison, ncol = 4) +
      theme(strip.text = element_text(size = 12))
    
    return(res)
  }
  
  
# Normalized gene count box plots ----------------------------------------------
  find_ensembl_ID <- function(search_name) {
    dat <- dplyr::filter(gene_IDs$AM14trans, external_gene_name == search_name)
    dat$ensembl_gene_id
  }
  
  coldata <- read.csv("raw_data/sample_info.csv")
  coldata <- coldata[, c(1:3, 5:10)]

  normcount_table <- function(gene_name, dds){
    ID <- find_ensembl_ID(gene_name)
    
    res <- counts(dds, normalized = TRUE)
    res <- res[rownames(res) %in% ID, ]
    res <- data.frame(Sample = names(res), counts = res)
    
    if(nrow(res) > 0){
      res <- left_join(res, coldata, by = "Sample")
      res <- res[ , c(1:2, 4, 6, 8, 9)] %>%
        group_by(Treatment, Drug, heatmap_col1)
      
      colnames(res) <- c("Sample", "counts", "Mouse_Number", "Treatment", "Drug", "Experiment")
      
      return(res)
    } else return(NULL)
    
    # summary_df <- res %>%
      # summarise(mean_counts = mean(counts), sd_counts = sd(counts),
      #           se_counts = sd_counts / sqrt(n()), .groups = 'drop')
    # 
    # anova_result <- aov(counts ~ Treatment, data = res)
    # tukey_result <- TukeyHSD(anova_result)
    # 
    # return(list(counts_table = res, summary = summary_df,
    #             ANOVA = anova_result, Tukey = tukey_result))
  }
  
  # is_outlier <- function(x){
  #   return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
  # }
  
  count_boxplot <- function(gene_name, counts_table, alt_title){
    if(missing(alt_title)) alt_title <- NULL
    
    if(is.null(alt_title)) alt_title <- gene_name
    
    ggboxplot(counts_table, x = "Drug", y = "counts", 
              # outliers = FALSE,
              add = c("jitter"), label = "Mouse_Number", repel = TRUE,
              ggtheme = theme_bw(), ylab = "counts", xlab = "Treatment") +
      # geom_point(aes(group = Treatment), position = position_jitter(width = 0.3), size = 1.5) +
      # geom_text(aes(label = Mouse_Number), vjust = -0.5) +
      labs(title = alt_title) + 
      facet_wrap( ~ Experiment, scales = 'fixed', ncol = 4)
  }
  
  multi_boxplot <- function(gene_name, alt_title, pl23_dds, mrl_dds, r848_dds, np_dds){
    if(missing(alt_title)) alt_title <- NULL
    if(missing(pl23_dds)) pl23_dds <- dds_AM14trans_PL23
    if(missing(mrl_dds)) mrl_dds <- dds_AM14MRLlpr
    if(missing(r848_dds)) r848_dds <- dds_AM14trans_R848
    if(missing(np_dds)) np_dds <- dds_B18trans
    
    pl23_normcounts <- normcount_table(gene_name, pl23_dds)
    mrl_normcounts <- normcount_table(gene_name, mrl_dds)
    r848_normcounts <- normcount_table(gene_name, r848_dds)
    np_normcounts <- normcount_table(gene_name, np_dds)
    
    if(!is.null(pl23_normcounts)) {
      full_table <- pl23_normcounts
      if(!is.null(mrl_normcounts)) full_table <- full_join(full_table, mrl_normcounts)
      if(!is.null(r848_normcounts)) full_table <- full_join(full_table, r848_normcounts)
      if(!is.null(np_normcounts)) full_table <- full_join(full_table, np_normcounts)
    } else if(!is.null(mrl_normcounts)){
      full_table <- mrl_normcounts
      if(!is.null(r848_normcounts)) full_table <- full_join(full_table, r848_normcounts)
      if(!is.null(np_normcounts)) full_table <- full_join(full_table, np_normcounts)
    } else if(!is.null(r848_normcounts)){
      full_table <- r848_normcounts
      if(!is.null(np_normcounts)) full_table <- full_join(full_table, np_normcounts)
    } else if(!is.null(np_normcounts)) {
      full_table <- np_normcounts
    } else {
      print("gene not found")
      return(NULL)
    }
    
    count_boxplot(gene_name, full_table, alt_title)
  }
  
  count_table_for_prism <- function(gene_name, dds, dds_results, cols_ctrl, cols_2DG){
    if(missing(cols_ctrl)) {
      print("missing cols_ctrl")
      return(NULL)
    }
    if(missing(cols_2DG)) {
      print("missing cols_2DG")
      return(NULL)
    }
    
    
    ID <- find_ensembl_ID(gene_name)
    gene_counts <- counts(dds, normalized = TRUE)
    gene_counts <- gene_counts[rownames(gene_counts) %in% ID, ]
    
    if(length(cols_ctrl) > length(cols_2DG)){
      length_diff <- length(cols_ctrl) - length(cols_2DG)
      
      ctrl <- gene_counts[cols_ctrl]
      dg <- c(gene_counts[cols_2DG], rep(NA, length_diff))
      gene_counts <- data.frame("Control" = ctrl, "2DG" = dg)
    } else if(length(cols_ctrl) < length(cols_2DG)){
      length_diff <- length(cols_2DG) - length(cols_ctrl)
      
      ctrl <- c(gene_counts[cols_ctrl], rep(NA, length_diff))
      dg <- gene_counts[cols_2DG]
      gene_counts <- data.frame("Control" = ctrl, "2DG" = dg)
    }
    
    rownames(gene_counts) <- c(1:nrow(gene_counts))
    return(gene_counts)
  }
  
  
  
  
  
  