setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_and_treatment_separated")
getwd()

corr_data = list(
  PL23_2DG_v_Ctrl = deseq_res$AM14transfer_PL23[, c(2, 8, 12)],
  R848_2DG_v_Ctrl = deseq_res$AM14transfer_R848[, c(2, 8, 12)],
  AM14MRLlpr = deseq_res$AM14MRLlpr[, c(2, 8, 12)],
  B18 = deseq_res$B18transfer[, c(2, 8, 12)])
head(corr_data$AM14MRLlpr)

# Function to find significant genes in two data sets and make one table -------  
  find_DEG_set_overlap <- function(deg1, colname1, deg2, colname2, p_value_cutoff, lfc_cutoff){
    if(missing(colname1)) colname1 <- "x"
    if(missing(colname2)) colname2 <- "y"
    if(missing(p_value_cutoff)) p_value_cutoff <- 0.05
    if(missing(lfc_cutoff)) lfc_cutoff <- 0
    
    filt1 <- deg1[deg1$padj < p_value_cutoff, ]
    filt1 <- filt1[abs(filt1$log2FoldChange) > lfc_cutoff, ]
    filt1 <- filt1[ , 1:2]
    colnames(filt1) <- c("Gene", colname1)
    
    filt2 <- deg2[deg2$padj < p_value_cutoff, ]
    filt2 <- filt2[abs(filt2$log2FoldChange) > lfc_cutoff, ]
    filt2 <- filt2[ , 1:2]
    colnames(filt2) <- c("Gene", colname2)
    
    merge_filt <- inner_join(filt1, filt2, by = "Gene")
    
    merge_filt$same_direction <- ifelse(
      test = merge_filt[, 2] < 0 & merge_filt[, 3] < 0, 
      yes = "down",
      no = ifelse(
        test = merge_filt[, 2] > 0 & merge_filt[, 3] > 0, 
        yes = "up", 
        no = NA
      )
    )
    
    merge_filt$ABC_direction <- ifelse(
      test = merge_filt$Gene %in% abc_up_genes,
      yes = "ABC_up_FoB_down",
      no = ifelse(
        test = merge_filt$Gene %in% abc_down_genes,
        yes = "ABC_down_FoB_up",
        no = NA
      )
    )
    
    return(merge_filt)
  }

# Generate a contingency table and run Fisher's exact test ---------------------
  contingency_test <- function(degs, extra_lfc_filter){
    if(missing(extra_lfc_filter)) extra_lfc_filter <- NULL
    
    if(!is.null(extra_lfc_filter)){
      degs <- degs[abs(degs[, 2]) >= extra_lfc_filter, ]
      degs <- degs[abs(degs[, 3]) >= extra_lfc_filter, ]
    }
    
    # 1. Classify each gene's direction in each dataset independently
      direction_table <- data.frame(
        x = factor(ifelse(degs[, 2] > 0, "up", "down"), levels = c("down", "up")),
        y = factor(ifelse(degs[, 3] > 0, "up", "down"), levels = c("down", "up"))
      )

    # 2. Build contingency table
      cont_table <- table(x = direction_table$x, 
                          y = direction_table$y)
      # print(cont_table)

    # 3. Fisher's exact test
      fisher_result <- fisher.test(cont_table)
      print(fisher_result)
      fisher_label <- paste0(
        "OR = ", round(fisher_result$estimate, 2),
        ", p = ", signif(fisher_result$p.value, 3),
        "\n(Fisher's exact, n = ", nrow(degs), ")"
      )
    
    return(fisher_label)
  }

# Correlation analyses ---------------------------------------------------------
  # AM14 PL2-3 vs AM14 MRLlpr --------------------------------------------------
    PL23_x_MRLlpr_nofilt <- find_DEG_set_overlap(corr_data$PL23_2DG_v_Ctrl, "PL23",
                                                 corr_data$AM14MRLlpr, "MRLlpr",
                                                 p_value_cutoff = 1)
    
    PL23_x_MRLlpr_corr <- cor.test(PL23_x_MRLlpr_nofilt$PL23, 
                                   PL23_x_MRLlpr_nofilt$MRLlpr, 
                                   method = "pearson")
    # PL23_x_MRLlpr_corr_r_val <- round(PL23_x_MRLlpr_corr$estimate, 3)
    # PL23_x_MRLlpr_corr_p_val <- signif(PL23_x_MRLlpr_corr$p.value, 3)
    # PL23_x_MRLlpr_corr_label <- paste0("R = ", PL23_x_MRLlpr_corr_r_val, 
    #                                    ", p = ", PL23_x_MRLlpr_corr_p_val, 
    #                                    "\n(all genes, n = ", nrow(PL23_x_MRLlpr_nofilt), ")")
    # PL23_x_MRLlpr_corr_label

  # AM14 PL2-3 vs AM14 R848 ----------------------------------------------------
    PL23_x_R848_nofilt <- find_DEG_set_overlap(corr_data$PL23_2DG_v_Ctrl, "PL23",
                                                 corr_data$R848_2DG_v_Ctrl, "R848",
                                                 p_value_cutoff = 1)
    
    PL23_x_R848_corr <- cor.test(PL23_x_R848_nofilt$PL23, 
                                   PL23_x_R848_nofilt$R848, 
                                   method = "pearson")
    # PL23_x_R848_corr_r_val <- round(PL23_x_R848_corr$estimate, 3)
    # PL23_x_R848_corr_p_val <- signif(PL23_x_R848_corr$p.value, 3)
    # PL23_x_R848_corr_label <- paste0("R = ", PL23_x_R848_corr_r_val, 
    #                                    ", p = ", PL23_x_R848_corr_p_val, 
    #                                    "\n(all genes, n = ", nrow(PL23_x_R848_nofilt), ")")
    # PL23_x_R848_corr_label
    
  # AM14 PL2-3 vs AM14 MRLlpr --------------------------------------------------
    PL23_x_NP_nofilt <- find_DEG_set_overlap(corr_data$PL23_2DG_v_Ctrl, "PL23",
                                                 corr_data$B18, "NP",
                                                 p_value_cutoff = 1)
    
    # PL23_x_NP_corr <- cor.test(PL23_x_NP_nofilt$PL23, 
    #                                PL23_x_NP_nofilt$NP, 
    #                                method = "pearson")
    # PL23_x_NP_corr_r_val <- round(PL23_x_NP_corr$estimate, 3)
    # PL23_x_NP_corr_p_val <- signif(PL23_x_NP_corr$p.value, 3)
    # PL23_x_NP_corr_label <- paste0("R = ", PL23_x_NP_corr_r_val, 
    #                                    ", p = ", PL23_x_NP_corr_p_val, 
    #                                    "\n(all genes, n = ", nrow(PL23_x_NP_nofilt), ")")
    # PL23_x_NP_corr_label
    
  # Summarize and save as an Excel file ----------------------------------------
    corr_summary <- data.frame(Comparison = c("AM14 PL2-3 +/- 2DG vs AM14 MRLlpr +/- 2DG",
                                              "AM14 PL2-3 +/- 2DG vs AM14 R848 +/- 2DG",
                                              "AM14 PL2-3 +/- 2DG vs B1-8 NP-OVA +/- 2DG"),
                               R2 = c(PL23_x_MRLlpr_corr$estimate,
                                      PL23_x_R848_corr$estimate,
                                      PL23_x_NP_corr$estimate),
                               p = c(PL23_x_MRLlpr_corr$p.value,
                                     PL23_x_R848_corr$p.value,
                                     PL23_x_NP_corr$p.value))
    corr_summary
    
    setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_and_treatment_separated")
    getwd()

    wb <- createWorkbook("Output/Correlation Scatterplots/Correlation_Results.xlsx")
    addWorksheet(wb, "Summary")
    writeData(wb, "Summary", corr_summary)
    saveWorkbook(wb, "Output/Correlation Scatterplots/Correlation_Results.xlsx", overwrite = TRUE)
    rm(wb)


    
# AM14 PL2-3 vs AM14 MRLlpr ----------------------------------------------------
  # No LFC filter --------------------------------------------------------------
    PL23_x_MRLlpr <- find_DEG_set_overlap(corr_data$PL23_2DG_v_Ctrl, "PL23",
                                          corr_data$AM14MRLlpr, "MRLlpr")
    PL23_x_MRLlpr_fish <- contingency_test(PL23_x_MRLlpr)
  
    ggplot(PL23_x_MRLlpr, aes(x = PL23, y = MRLlpr, fill = ABC_direction)) +
      geom_point(aes(color = same_direction)) +
      scale_color_manual(
        values = c("down" = "#00496f", "up" = "#F21A00", "neither" = "darkgray"),
        na.value = "gray",
        guide = NULL # guide = guide_legend(title = "title")
      ) +
      geom_hline(yintercept = c(0), linetype = "dotted", color = "#555555") +
      geom_vline(xintercept = c(0), linetype = "dotted", color = "#555555") +
      geom_text_repel(
        data = PL23_x_MRLlpr[is.na(PL23_x_MRLlpr$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = 20
      ) +
      geom_label_repel(
        data = PL23_x_MRLlpr[!is.na(PL23_x_MRLlpr$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = Inf
      ) +      
      scale_fill_manual(
        values = c("ABC_up_FoB_down" = "#ed8b00", "ABC_down_FoB_up" = "#78B7C5"),
        na.value = "white"
      ) +
      labs(
        x = "LFC: AM14 PL2-3 2DG vs Ctrl",
        y = "LFC: AM14 MRL/lpr 2DG vs Ctrl",
        title = "AM14 PL2-3 +/- 2DG vs AM14 MRLlpr +/- 2DG"
      ) +
      annotate("text",
               x = -Inf, y = Inf,
               hjust = -0.1, vjust = 1.5,
               label = PL23_x_MRLlpr_fish,
               size = 4) +
      theme_classic()
    
    ggsave(filename = "Output/Correlation Scatterplots/Contingency Test - AM14 PL2-3 vs AM14 MRLlpr - no LFC filter.png", 
           width = 8, height = 7, units = "in", dpi = 600)
  
  # |LFC| > 0.6 ----------------------------------------------------------------
    PL23_x_MRLlpr_0.6 <- find_DEG_set_overlap(corr_data$PL23_2DG_v_Ctrl, "PL23",
                                              corr_data$AM14MRLlpr, "MRLlpr",
                                              lfc_cutoff = 0.6)
    PL23_x_MRLlpr_0.6_fish <- contingency_test(PL23_x_MRLlpr_0.6)
    
    ggplot(PL23_x_MRLlpr_0.6, aes(x = PL23, y = MRLlpr, fill = ABC_direction)) +
      geom_point(aes(color = same_direction)) +
      scale_color_manual(
        values = c("down" = "#00496f", "up" = "#F21A00", "neither" = "darkgray"),
        na.value = "gray",
        guide = NULL # guide = guide_legend(title = "title")
      ) +
      geom_hline(yintercept = c(-0.6, 0.6), linetype = "dotted", color = "#555555") +
      geom_vline(xintercept = c(-0.6, 0.6), linetype = "dotted", color = "#555555") +
      geom_text_repel(
        data = PL23_x_MRLlpr_0.6[is.na(PL23_x_MRLlpr_0.6$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = 20
      ) +
      geom_label_repel(
        data = PL23_x_MRLlpr_0.6[!is.na(PL23_x_MRLlpr_0.6$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = Inf
      ) +      
      scale_fill_manual(
        values = c("ABC_up_FoB_down" = "#ed8b00", "ABC_down_FoB_up" = "#78B7C5"),
        na.value = "white"
      ) +
      labs(
        x = "LFC: AM14 PL2-3 2DG vs Ctrl",
        y = "LFC: AM14 MRL/lpr 2DG vs Ctrl",
        title = "AM14 PL2-3 +/- 2DG vs AM14 MRLlpr +/- 2DG"
      ) +
      annotate("text",
               x = -Inf, y = Inf,
               hjust = -0.1, vjust = 1.5,
               label = PL23_x_MRLlpr_0.6_fish,
               size = 4) +
      theme_classic()
    
    ggsave(filename = "Output/Correlation Scatterplots/Contingency Test - AM14 PL2-3 vs AM14 MRLlpr - LFC filter 0.6.png", 
           width = 8, height = 7, units = "in", dpi = 600)
    
  
  
# AM14 PL2-3 vs AM14 R848 ------------------------------------------------------
  # No LFC filter --------------------------------------------------------------
    PL23_x_R848 <- find_DEG_set_overlap(corr_data$PL23_2DG_v_Ctrl, "PL23",
                                        corr_data$R848_2DG_v_Ctrl, "R848")
    PL23_x_R848_fish <- contingency_test(PL23_x_R848)
    
    ggplot(PL23_x_R848, aes(x = PL23, y = R848, fill = ABC_direction)) +
      geom_point(aes(color = same_direction)) +
      scale_color_manual(
        values = c("down" = "#00496f", "up" = "#F21A00", "neither" = "darkgray"),
        na.value = "gray",
        guide = NULL # guide = guide_legend(title = "title")
      ) +
      geom_hline(yintercept = c(0), linetype = "dotted", color = "#555555") +
      geom_vline(xintercept = c(0), linetype = "dotted", color = "#555555") +
      geom_text_repel(
        data = PL23_x_R848[is.na(PL23_x_R848$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = 20
      ) +
      geom_label_repel(
        data = PL23_x_R848[!is.na(PL23_x_R848$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = Inf
      ) +      
      scale_fill_manual(
        values = c("ABC_up_FoB_down" = "#ed8b00", "ABC_down_FoB_up" = "#78B7C5"),
        na.value = "white"
      ) +
      labs(
        x = "LFC: AM14 PL2-3 2DG vs Ctrl",
        y = "LFC: AM14 R848 2DG vs Ctrl",
        title = "AM14 PL2-3 +/- 2DG vs AM14 R848 +/- 2DG"
      ) +
      annotate("text",
               x = -Inf, y = Inf,
               hjust = -0.1, vjust = 1.5,
               label = PL23_x_R848_fish,
               size = 4) +
      theme_classic()
    
    ggsave(filename = "Output/Correlation Scatterplots/Contingency Test - AM14 PL2-3 vs AM14 R848 - no LFC filter.png", 
           width = 8, height = 7, units = "in", dpi = 600)
    
  # |LFC| > 0.6 ----------------------------------------------------------------
    PL23_x_R848_0.6 <- find_DEG_set_overlap(corr_data$PL23_2DG_v_Ctrl, "PL23",
                                              corr_data$R848_2DG_v_Ctrl, "R848",
                                              lfc_cutoff = 0.6)
    PL23_x_R848_0.6_fish <- contingency_test(PL23_x_R848_0.6)
    
    ggplot(PL23_x_R848_0.6, aes(x = PL23, y = R848, fill = ABC_direction)) +
      geom_point(aes(color = same_direction)) +
      scale_color_manual(
        values = c("down" = "#00496f", "up" = "#F21A00", "neither" = "darkgray"),
        na.value = "gray",
        guide = NULL # guide = guide_legend(title = "title")
      ) +
      geom_hline(yintercept = c(-0.6, 0.6), linetype = "dotted", color = "#555555") +
      geom_vline(xintercept = c(-0.6, 0.6), linetype = "dotted", color = "#555555") +
      geom_text_repel(
        data = PL23_x_R848_0.6[is.na(PL23_x_R848_0.6$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = 20
      ) +
      geom_label_repel(
        data = PL23_x_R848_0.6[!is.na(PL23_x_R848_0.6$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = Inf
      ) +      
      scale_fill_manual(
        values = c("ABC_up_FoB_down" = "#ed8b00", "ABC_down_FoB_up" = "#78B7C5"),
        na.value = "white"
      ) +
      labs(
        x = "LFC: AM14 PL2-3 2DG vs Ctrl",
        y = "LFC: AM14 R848 2DG vs Ctrl",
        title = "AM14 PL2-3 +/- 2DG vs AM14 R848 +/- 2DG"
      ) +
      annotate("text",
               x = -Inf, y = Inf,
               hjust = -0.1, vjust = 1.5,
               label = PL23_x_R848_0.6_fish,
               size = 4) +
      theme_classic()
    
    ggsave(filename = "Output/Correlation Scatterplots/Contingency Test - AM14 PL2-3 vs AM14 R848 - LFC filter 0.6.png", 
           width = 8, height = 7, units = "in", dpi = 600)
    
  
# AM14 PL2-3 vs B1-8 NP --------------------------------------------------------
  # No LFC filter --------------------------------------------------------------
    PL23_x_NP <- find_DEG_set_overlap(corr_data$PL23_2DG_v_Ctrl, "PL23",
                                      corr_data$B18, "NP")
    PL23_x_NP_fish <- contingency_test(PL23_x_NP)
    
    ggplot(PL23_x_NP, aes(x = PL23, y = NP, fill = ABC_direction)) +
      geom_point(aes(color = same_direction)) +
      scale_color_manual(
        values = c("down" = "#00496f", "up" = "#F21A00", "neither" = "darkgray"),
        na.value = "gray",
        guide = NULL # guide = guide_legend(title = "title")
      ) +
      geom_hline(yintercept = c(0), linetype = "dotted", color = "#555555") +
      geom_vline(xintercept = c(0), linetype = "dotted", color = "#555555") +
      geom_text_repel(
        data = PL23_x_NP[is.na(PL23_x_NP$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = 20
      ) +
      geom_label_repel(
        data = PL23_x_NP[!is.na(PL23_x_NP$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = Inf
      ) +      
      scale_fill_manual(
        values = c("ABC_up_FoB_down" = "#ed8b00", "ABC_down_FoB_up" = "#78B7C5"),
        na.value = "white"
      ) +
      labs(
        x = "LFC: AM14 PL2-3 2DG vs Ctrl",
        y = "LFC: B1-8 NP-OVA 2DG vs Ctrl",
        title = "AM14 PL2-3 +/- 2DG vs B1-8 NP-OVA +/- 2DG"
      ) +
      annotate("text",
               x = -Inf, y = Inf,
               hjust = -0.1, vjust = 1.5,
               label = PL23_x_NP_fish,
               size = 4) +
      theme_classic()
    
    ggsave(filename = "Output/Correlation Scatterplots/Contingency Test - AM14 PL2-3 vs B1-8 NP - no LFC filter.png", 
           width = 8, height = 7, units = "in", dpi = 600)
    
  # |LFC| > 0.6 ----------------------------------------------------------------
    PL23_x_NP_0.6 <- find_DEG_set_overlap(corr_data$PL23_2DG_v_Ctrl, "PL23",
                                            corr_data$B18, "NP",
                                            lfc_cutoff = 0.6)
    PL23_x_NP_0.6_fish <- contingency_test(PL23_x_NP_0.6)
    
    ggplot(PL23_x_NP_0.6, aes(x = PL23, y = NP, fill = ABC_direction)) +
      geom_point(aes(color = same_direction)) +
      scale_color_manual(
        values = c("down" = "#00496f", "up" = "#F21A00", "neither" = "darkgray"),
        na.value = "gray",
        guide = NULL # guide = guide_legend(title = "title")
      ) +
      geom_hline(yintercept = c(-0.6, 0.6), linetype = "dotted", color = "#555555") +
      geom_vline(xintercept = c(-0.6, 0.6), linetype = "dotted", color = "#555555") +
      geom_text_repel(
        data = PL23_x_NP_0.6[is.na(PL23_x_NP_0.6$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = 20
      ) +
      geom_label_repel(
        data = PL23_x_NP_0.6[!is.na(PL23_x_NP_0.6$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = Inf
      ) +      
      scale_fill_manual(
        values = c("ABC_up_FoB_down" = "#ed8b00", "ABC_down_FoB_up" = "#78B7C5"),
        na.value = "white"
      ) +
      labs(
        x = "LFC: AM14 PL2-3 2DG vs Ctrl",
        y = "LFC: B1-8 NP-OVA 2DG vs Ctrl",
        title = "AM14 PL2-3 +/- 2DG vs B1-8 NP-OVA +/- 2DG"
      ) +
      annotate("text",
               x = -Inf, y = Inf,
               hjust = -0.1, vjust = 1.5,
               label = PL23_x_NP_0.6_fish,
               size = 4) +
      theme_classic()
    
    ggsave(filename = "Output/Correlation Scatterplots/Contingency Test - AM14 PL2-3 vs B1-8 NP - LFC filter 0.6.png", 
           width = 8, height = 7, units = "in", dpi = 600)
    
  
  
  


