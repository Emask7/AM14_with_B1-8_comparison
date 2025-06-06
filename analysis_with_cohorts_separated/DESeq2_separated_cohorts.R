# Import raw counts from Rosalind output ---------------------------------------
  cts_AM14trans <- import_Rosalind_data("raw_data/AM14_Adoptive_Transfer_rawCountsWithAnnotations.txt")
  cts_B18trans <- import_Rosalind_data("raw_data/B1-8_Adoptive_Transfer_rawCountsWithAnnotations.txt")
  cts_AM14MRLlpr <- import_Rosalind_data("raw_data/AM14_MRLlpr_rawCountsWithAnnotations.txt")

  # Save a table of gene IDs in case they need to be converted later -----------
    gene_IDs <- list(AM14trans = cts_AM14trans[, c(1:4)],
                     B18trans = cts_B18trans[, c(1:4)],
                     AM14MRLlpr = cts_AM14MRLlpr[, c(1:4)])
        
  # Remove extra columns from cts ----------------------------------------------
    cts_AM14trans <- cts_AM14trans[, c(5:24)]
    cts_B18trans <- cts_B18trans[, c(5:13)]
    cts_AM14MRLlpr <- cts_AM14MRLlpr[, c(5:21)]
  
    head(cts_AM14trans)
    head(cts_B18trans)
    head(cts_AM14MRLlpr)

# Set up experimental factors --------------------------------------------------
  coldata <- read.csv("raw_data/sample_info.csv")
  coldata <- coldata[, c(1:3, 5:6, 10)]
  coldata
  
  cd_AM14trans <- coldata[1:20, ]
  cd_AM14trans$Cohort <- factor(cd_AM14trans$Cohort)
  cd_AM14trans$Treatment <- factor(cd_AM14trans$Treatment)
  cd_AM14trans
  
  cd_AM14trans_co1 <- filter(cd_AM14trans, Cohort == "A1")
  cd_AM14trans_co1
  cd_AM14trans_co2 <- filter(cd_AM14trans, Cohort == "A2")
  cd_AM14trans_co2
  
  cts_AM14trans_co1 <- cts_AM14trans[, c(1:3, 6:8, 11, 12, 16, 17)]
  head(cts_AM14trans_co1)
  
  cts_AM14trans_co2 <- cts_AM14trans[, c(4, 5, 9, 10, 13:15, 18:20)]
  head(cts_AM14trans_co2)
  
  
  cd_B18trans <- coldata[21:29, ]
  cd_B18trans$Cohort <- factor(cd_B18trans$Cohort)
  cd_B18trans$Treatment <- factor(cd_B18trans$Treatment)
  cd_B18trans
  
  cd_AM14MRLlpr <- coldata[30:46, ]
  cd_AM14MRLlpr$Cohort <- factor(cd_AM14MRLlpr$Cohort)
  cd_AM14MRLlpr$Treatment <- factor(cd_AM14MRLlpr$Treatment)
  cd_AM14MRLlpr
  
  cd_AM14MRLlpr_co1 <- filter(cd_AM14MRLlpr, Cohort == "M1")
  cd_AM14MRLlpr_co1
  cd_AM14MRLlpr_co2 <- filter(cd_AM14MRLlpr, Cohort == "M2")
  cd_AM14MRLlpr_co2
  
  cts_AM14MRLlpr_co1 <- cts_AM14MRLlpr[, c(1, 2, 9, 10)]
  head(cts_AM14MRLlpr_co1)
  
  cts_AM14MRLlpr_co2 <- cts_AM14MRLlpr[, c(3:8, 11:17)]
  head(cts_AM14MRLlpr_co2)
  

  # Make sure the columns of cts and rows of coldata are in the same order -----
    for (x in 1:ncol(cts_AM14trans_co1)) {
      if(colnames(cts_AM14trans_co1)[x] == cd_AM14trans_co1$Sample[x]) print(c(x, "true"))
      else print(c(x, "false"))
    }
    for (x in 1:ncol(cts_AM14trans_co2)) {
      if(colnames(cts_AM14trans_co2)[x] == cd_AM14trans_co2$Sample[x]) print(c(x, "true"))
      else print(c(x, "false"))
    }
    for (x in 1:ncol(cts_B18trans)) {
      if(colnames(cts_B18trans)[x] == cd_B18trans$Sample[x]) print(c(x, "true"))
      else print(c(x, "false"))
    }
    for (x in 1:ncol(cts_AM14MRLlpr)) {
      if(colnames(cts_AM14MRLlpr)[x] == cd_AM14MRLlpr$Sample[x]) print(c(x, "true"))
      else print(c(x, "false"))
    }
    rm(x)
    
    for (x in 1:ncol(cts_AM14MRLlpr_co1)) {
      if(colnames(cts_AM14MRLlpr_co1)[x] == cd_AM14MRLlpr_co1$Sample[x]) print(c(x, "true"))
      else print(c(x, "false"))
    }
    for (x in 1:ncol(cts_AM14MRLlpr_co2)) {
      if(colnames(cts_AM14MRLlpr_co2)[x] == cd_AM14MRLlpr_co2$Sample[x]) print(c(x, "true"))
      else print(c(x, "false"))
    }
    

# Set up DESeq Data Set --------------------------------------------------------
    dds_AM14trans_co1 <- DESeqDataSetFromMatrix(countData = cts_AM14trans_co1, 
                                                colData = cd_AM14trans_co1,
                                                design = ~ Treatment)
    dds_AM14trans_co1
    keep <- rowSums(counts(dds_AM14trans_co1) >= 10) >= 2
    dds_AM14trans_co1 <- dds_AM14trans_co1[keep,]
    dds_AM14trans_co1 <- DESeq(dds_AM14trans_co1)
    resultsNames(dds_AM14trans_co1)
    rm(keep)
    
    dds_AM14trans_co2 <- DESeqDataSetFromMatrix(countData = cts_AM14trans_co2, 
                                                colData = cd_AM14trans_co2,
                                                design = ~ Treatment)
    dds_AM14trans_co2
    keep <- rowSums(counts(dds_AM14trans_co2) >= 10) >= 2
    dds_AM14trans_co2 <- dds_AM14trans_co2[keep,]
    dds_AM14trans_co2 <- DESeq(dds_AM14trans_co2)
    resultsNames(dds_AM14trans_co2)
    rm(keep)
    
    dds_B18trans <- DESeqDataSetFromMatrix(countData = cts_B18trans, 
                                           colData = cd_B18trans,
                                           design = ~Treatment)
    dds_B18trans
    keep <- rowSums(counts(dds_B18trans) >= 10) >= 4
    dds_B18trans <- dds_B18trans[keep,]
    dds_B18trans <- DESeq(dds_B18trans)
    resultsNames(dds_B18trans)
    rm(keep)
    
    dds_AM14MRLlpr_co1 <- DESeqDataSetFromMatrix(countData = cts_AM14MRLlpr_co1,
                                                 colData = cd_AM14MRLlpr_co1,
                                                 design = ~Treatment)
    dds_AM14MRLlpr_co1
    keep <- rowSums(counts(dds_AM14MRLlpr_co1) >= 10) >= 2
    dds_AM14MRLlpr_co1 <- dds_AM14MRLlpr_co1[keep,]
    dds_AM14MRLlpr_co1 <- DESeq(dds_AM14MRLlpr_co1)
    resultsNames(dds_AM14MRLlpr_co1)
    rm(keep)
    
    dds_AM14MRLlpr_co2 <- DESeqDataSetFromMatrix(countData = cts_AM14MRLlpr_co2,
                                                 colData = cd_AM14MRLlpr_co2,
                                                 design = ~Treatment)
    dds_AM14MRLlpr_co2
    keep <- rowSums(counts(dds_AM14MRLlpr_co2) >= 10) >= 6
    dds_AM14MRLlpr_co2 <- dds_AM14MRLlpr_co2[keep,]
    dds_AM14MRLlpr_co2 <- DESeq(dds_AM14MRLlpr_co2)
    resultsNames(dds_AM14MRLlpr_co2)
    rm(keep)
    
    
# Run QC steps -----------------------------------------------------------------
  setwd("./analysis_with_cohorts_separated/")
  getwd()
    
  QC_heatmaps(dds_AM14trans_co1, "AM14_Adoptive_Transfer_Cohort1", "AM14 Adoptive Transfer (Cohort 1)", plot_cohorts = FALSE)
  QC_heatmaps(dds_AM14trans_co2, "AM14_Adoptive_Transfer_Cohort2", "AM14 Adoptive Transfer (Cohort 2)", plot_cohorts = FALSE)
  QC_heatmaps(dds_B18trans, "B18_Adoptive_Transfer", "B1-8 Adoptive Transfer", plot_cohorts = FALSE)
  QC_heatmaps(dds_AM14MRLlpr_co1, "AM14_MRLlpr_Cohort1", "AM14 MRL/lpr +/- 2DG (Cohort 1)", plot_cohorts = FALSE)
  QC_heatmaps(dds_AM14MRLlpr_co2, "AM14_MRLlpr_Cohort2", "AM14 MRL/lpr +/- 2DG (Cohort 2)", plot_cohorts = FALSE)
  

  QC_PCAplot(dds_AM14trans_co1, "AM14_Adoptive_Transfer-Cohort1",
             "AM14 Adoptive Transfer\nCohort 1",
             batch_effect = NULL, draw_ellipse = TRUE)
  dev.off()
  
  QC_PCAplot(dds_AM14trans_co2, "AM14_Adoptive_Transfer-Cohort2",
             "AM14 Adoptive Transfer\nCohort 2",
             batch_effect = NULL, draw_ellipse = TRUE)
  dev.off()
  
  QC_PCAplot(dds_B18trans, "B1-8_Adoptive_Transfer", "B1-8 Adoptive Transfer",
             batch_effect = NULL, draw_ellipse = TRUE)
  dev.off()

  QC_PCAplot(dds_AM14MRLlpr_co1, "AM14_MRLlpr-Cohort1",
             "AM14 MRL/lpr mice +/- 2DG\nCohort 1",
             batch_effect = NULL, draw_ellipse = TRUE)
  dev.off()

  QC_PCAplot(dds_AM14MRLlpr_co2, "AM14_MRLlpr-Cohort2",
             "AM14 MRL/lpr mice +/- 2DG\nCohort 2",
             batch_effect = NULL, draw_ellipse = TRUE)
  dev.off()
  
# Differential Expression Analysis ---------------------------------------------
  resultsNames(dds_AM14trans_co1)
  resultsNames(dds_AM14trans_co2)
  resultsNames(dds_B18trans)
  resultsNames(dds_AM14MRLlpr_co1)
  resultsNames(dds_AM14MRLlpr_co2)
  
# Get DESeq2 results -----------------------------------------------------------
  deseq_res <- list(
    AM14transfer = list(
      co1 = list(
        PL23_2DG_v_Ctrl = results_wrapper(dds_AM14trans_co1, c("Treatment", "PL23_2DG", "PL23"), gene_IDs$AM14trans),
        R848_2DG_v_Ctrl = results_wrapper(dds_AM14trans_co1, c("Treatment", "R848_2DG", "R848"), gene_IDs$AM14trans),
        PL23_vs_R848 = results_wrapper(dds_AM14trans_co1, c("Treatment", "PL23", "R848"), gene_IDs$AM14trans)
      ),
      co2 = list(
        PL23_2DG_v_Ctrl = results_wrapper(dds_AM14trans_co2, c("Treatment", "PL23_2DG", "PL23"), gene_IDs$AM14trans),
        R848_2DG_v_Ctrl = results_wrapper(dds_AM14trans_co2, c("Treatment", "R848_2DG", "R848"), gene_IDs$AM14trans),
        PL23_vs_R848 = results_wrapper(dds_AM14trans_co2, c("Treatment", "PL23", "R848"), gene_IDs$AM14trans)
      )
  ),
    # B18transfer = results_wrapper(dds_B18trans, c("Treatment", "NP_2DG", "NP"), gene_IDs$B18trans),
    AM14MRLlpr = list(
      co1 = results_wrapper(dds_AM14MRLlpr_co1, c("Treatment", "2DG", "Control"), gene_IDs$AM14MRLlpr),
      co2 = results_wrapper(dds_AM14MRLlpr_co2, c("Treatment", "2DG", "Control"), gene_IDs$AM14MRLlpr)
    )
  )
  
  head(deseq_res$AM14transfer$co1$PL23_2DG_v_Ctrl)
  # head(deseq_res$B18transfer)
  # head(deseq_res$AM14MRLlpr)
  
  # DEG Heatmaps ---------------------------------------------------------------
    AM14_DEGs_co1 <- DEG_list(list(deseq_res$AM14transfer$co1$PL23_2DG_v_Ctrl,
                                   deseq_res$AM14transfer$co1$R848_2DG_v_Ctrl,
                                   deseq_res$AM14transfer$co1$PL23_vs_R848))
    head(AM14_DEGs_co1)
    
    AM14_DEGs_co2 <- DEG_list(list(deseq_res$AM14transfer$co2$PL23_2DG_v_Ctrl,
                                   deseq_res$AM14transfer$co2$R848_2DG_v_Ctrl,
                                   deseq_res$AM14transfer$co2$PL23_vs_R848))
    head(AM14_DEGs_co2)
  
    AM14_DEGs_PL23only_co1 <- DEG_list(list(deseq_res$AM14transfer$co1$PL23_2DG_v_Ctrl))
    head(AM14_DEGs_PL23only_co1)
    
    AM14_DEGs_PL23only_co2 <- DEG_list(list(deseq_res$AM14transfer$co2$PL23_2DG_v_Ctrl))
    head(AM14_DEGs_PL23only_co2)
    
    colors_AM14 <- list(Treatment = c(PL23 = "#00496f", PL23_2DG = "#0f85a0", 
                                      R848 = "#edd746", R848_2DG = "#dd4124"))
    
    DEG_heatmap(dds_AM14trans_co1, AM14_DEGs_PL23only_co1, colors_AM14,
                "AM14 Adoptive Transfer: Cohort 1\n(PL2-3 +/- 2DG Significant DEGs)",
                "AM14 adoptive transfer - Cohort 1 - PL2-3 only",
                unsupervised = TRUE, h = 1500, w = 1250)
    
    DEG_heatmap(dds_AM14trans_co2, AM14_DEGs_PL23only_co2, colors_AM14,
                "AM14 Adoptive Transfer: Cohort 2\n(PL2-3 +/- 2DG Significant DEGs)",
                "AM14 adoptive transfer - Cohort 2 - PL2-3 only",
                unsupervised = TRUE, h = 1500, w = 1250)
    

    DEG_heatmap(dds_AM14trans_co1, AM14_DEGs_co1, colors_AM14,
                "AM14 Adoptive Transfer: Cohort 1\n(All 2DG Significant DEGs)",
                "AM14 adoptive transfer - Cohort 1",
                unsupervised = TRUE, h = 1500, w = 1250)
    
    DEG_heatmap(dds_AM14trans_co2, AM14_DEGs_co2, colors_AM14,
                "AM14 Adoptive Transfer: Cohort 2\n(All Significant DEGs)",
                "AM14 adoptive transfer - Cohort 2",
                unsupervised = TRUE, h = 1500, w = 1250)
    
    DEG_heatmap(dds_B18trans, DEG_list(list(deseq_res$B18transfer)), 
                list(Treatment = c(NP = "#0f85a0", NP_2DG = "#dd4124")), 
                "B1-8 Adoptive Transfer\n(Significant DEGs)",
                "B1-8 adoptive transfer DEG heatmap",
                unsupervised = TRUE, h = 1000, w = 1000)
    
    DEG_heatmap(dds_AM14MRLlpr_co1, DEG_list(list(deseq_res$AM14MRLlpr$co1)), 
                list(Treatment = c(Control = "#0f85a0", "2DG" = "#dd4124")), 
                "AM14 MRL/lpr Mice: Cohort 1\n(Significant DEGs)",
                "AM14 MRLlpr DEG heatmap - Cohort 1", 
                unsupervised = TRUE, h = 1500, w = 1250)
    
    DEG_heatmap(dds_AM14MRLlpr_co2, DEG_list(list(deseq_res$AM14MRLlpr$co2)), 
                list(Treatment = c(Control = "#0f85a0", "2DG" = "#dd4124")), 
                "AM14 MRL/lpr Mice: Cohort 2\n(Significant DEGs)",
                "AM14 MRLlpr DEG heatmap - Cohort 2", 
                unsupervised = TRUE, h = 1500, w = 1250)
    

# Make a data.frame that summarizes the numbers of DEGs detected -------------
  full_summary <- join_all(list(summary_wrapper(deseq_res$AM14transfer$co1$PL23_2DG_v_Ctrl, "PL2-3+2DG vs PL2-3", "AM14 Adoptive Transfer: Cohort 1", 0.05, 1),
                                summary_wrapper(deseq_res$AM14transfer$co2$PL23_2DG_v_Ctrl, "PL2-3+2DG vs PL2-3", "AM14 Adoptive Transfer: Cohort 2", 0.05, 1),
                                summary_wrapper(deseq_res$AM14transfer$co1$R848_2DG_v_Ctrl, "R848+2DG vs R848", "AM14 Adoptive Transfer: Cohort 1", 0.05, 1),
                                summary_wrapper(deseq_res$AM14transfer$co2$R848_2DG_v_Ctrl, "R848+2DG vs R848", "AM14 Adoptive Transfer: Cohort 2", 0.05, 1),
                                summary_wrapper(deseq_res$AM14transfer$co1$PL23_vs_R848, "PL2-3 vs R848", "AM14 Adoptive Transfer: Cohort 1", 0.05, 1),
                                summary_wrapper(deseq_res$AM14transfer$co2$PL23_vs_R848, "PL2-3 vs R848", "AM14 Adoptive Transfer: Cohort 2", 0.05, 1),
                                summary_wrapper(deseq_res$B18transfer, "NP+2DG vs NP", "B1-8 Adoptive Transfer", 0.05, 1),
                                summary_wrapper(deseq_res$AM14MRLlpr$co1, "2DG vs Control", "AM14 MRL/lpr: Cohort 1", 0.05, 1),
                                summary_wrapper(deseq_res$AM14MRLlpr$co2, "2DG vs Control", "AM14 MRL/lpr: Cohort 2", 0.05, 1),
                                summary_wrapper(deseq_res$AM14transfer$co1$PL23_2DG_v_Ctrl, "PL2-3+2DG vs PL2-3", "AM14 Adoptive Transfer: Cohort 1", 0.05, 0.5),
                                summary_wrapper(deseq_res$AM14transfer$co2$PL23_2DG_v_Ctrl, "PL2-3+2DG vs PL2-3", "AM14 Adoptive Transfer: Cohort 2", 0.05, 0.5),
                                summary_wrapper(deseq_res$AM14transfer$co1$R848_2DG_v_Ctrl, "R848+2DG vs R848", "AM14 Adoptive Transfer: Cohort 1", 0.05, 0.5),
                                summary_wrapper(deseq_res$AM14transfer$co2$R848_2DG_v_Ctrl, "R848+2DG vs R848", "AM14 Adoptive Transfer: Cohort 2", 0.05, 0.5),
                                summary_wrapper(deseq_res$AM14transfer$co1$PL23_vs_R848, "PL2-3 vs R848", "AM14 Adoptive Transfer: Cohort 1", 0.05, 0.5),
                                summary_wrapper(deseq_res$AM14transfer$co2$PL23_vs_R848, "PL2-3 vs R848", "AM14 Adoptive Transfer: Cohort 2", 0.05, 0.5),
                                summary_wrapper(deseq_res$B18transfer, "NP+2DG vs NP", "B1-8 Adoptive Transfer", 0.05, 0.5),
                                summary_wrapper(deseq_res$AM14MRLlpr$co1, "2DG vs Control", "AM14 MRL/lpr: Cohort 1", 0.05, 0.5),
                                summary_wrapper(deseq_res$AM14MRLlpr$co2, "2DG vs Control", "AM14 MRL/lpr: Cohort 2", 0.05, 0.5)),
                           type = "full")
  full_summary

# Save results to an Excel file ----------------------------------------------
  wb <- createWorkbook("Output/DESeq2_results.xlsx")
  
  addWorksheet(wb, "DEG Summaries")
  writeData(wb, "DEG Summaries", full_summary)
  
  addWorksheet(wb, "AM14 Co1 - PL23_2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 Co1 - R848_2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 Co1 - PL23_vs_R848")
  addWorksheet(wb, "B1-8 - 2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 MRLlpr Co1 - 2DG_vs_Ctrl")
  
  addWorksheet(wb, "AM14 Co2 - PL23_2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 Co2 - R848_2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 Co2 - PL23_vs_R848")
  addWorksheet(wb, "AM14 MRLlpr Co2 - 2DG_vs_Ctrl")
  
  writeData(wb, "AM14 Co1 - PL23_2DG_vs_Ctrl", deseq_res$AM14transfer$co1$PL23_2DG_v_Ctrl)
  writeData(wb, "AM14 Co1 - R848_2DG_vs_Ctrl", deseq_res$AM14transfer$co1$R848_2DG_v_Ctrl)
  writeData(wb, "AM14 Co1 - PL23_vs_R848", deseq_res$AM14transfer$co1$PL23_vs_R848)
  writeData(wb, "B1-8 - 2DG_vs_Ctrl", deseq_res$B18transfer)
  writeData(wb, "AM14 MRLlpr Co1 - 2DG_vs_Ctrl", deseq_res$AM14MRLlpr$co1)
  
  writeData(wb, "AM14 Co2 - PL23_2DG_vs_Ctrl", deseq_res$AM14transfer$co2$PL23_2DG_v_Ctrl)
  writeData(wb, "AM14 Co2 - R848_2DG_vs_Ctrl", deseq_res$AM14transfer$co2$R848_2DG_v_Ctrl)
  writeData(wb, "AM14 Co2 - PL23_vs_R848", deseq_res$AM14transfer$co2$PL23_vs_R848)
  writeData(wb, "AM14 MRLlpr Co2 - 2DG_vs_Ctrl", deseq_res$AM14MRLlpr$co2)
  
  
  saveWorkbook(wb, "Output/DESeq2_results.xlsx", overwrite = TRUE)
  rm(wb)
  
  # wb <- createWorkbook("Output/DESeq2_gene_counts.xlsx")
  # 
  # addWorksheet(wb, "Raw_Counts - AM14 transfer")
  # addWorksheet(wb, "Norm_Counts - AM14 transfer")
  # addWorksheet(wb, "Raw_Counts - B1-8 transfer")
  # addWorksheet(wb, "Norm_Counts - B1-8 transfer")
  # addWorksheet(wb, "Raw_Counts - AM14 MRLlpr")
  # addWorksheet(wb, "Norm_Counts - AM14 MRLlpr")
  # 
  # writeData(wb, "Raw_Counts - AM14 transfer", counts(dds_AM14trans_co1, normalized = FALSE), rowNames = TRUE)
  # writeData(wb, "Norm_Counts - AM14 transfer", counts(dds_AM14trans_co1, normalized = TRUE), rowNames = TRUE)
  # writeData(wb, "Raw_Counts - B1-8 transfer", counts(dds_B18trans, normalized = FALSE), rowNames = TRUE)
  # writeData(wb, "Norm_Counts - B1-8 transfer", counts(dds_B18trans, normalized = TRUE), rowNames = TRUE)
  # writeData(wb, "Raw_Counts - AM14 MRLlpr", counts(dds_AM14MRLlpr, normalized = FALSE), rowNames = TRUE)
  # writeData(wb, "Norm_Counts - AM14 MRLlpr", counts(dds_AM14MRLlpr, normalized = TRUE), rowNames = TRUE)
  # 
  # saveWorkbook(wb, "Output/DESeq2_gene_counts.xlsx", overwrite = TRUE)
  # rm(wb)

# Make DEG lists and export to files -------------------------------------------
  # sig_DEGs <- list(
  #   AM14transfer = list(
  #     PL23_2DG_v_Ctrl = sig_DEG_table(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "all", 1, 0.05, "AM14 - PL23_2DG_vs_Ctrl"),
  #     R848_2DG_v_Ctrl = sig_DEG_table(deseq_res$AM14transfer$R848_2DG_v_Ctrl, "all", 1, 0.05, "AM14 - R848_2DG_vs_Ctrl"),
  #     PL23_vs_R848 = sig_DEG_table(deseq_res$AM14transfer$PL23_vs_R848, "all", 1, 0.05, "AM14 - PL23_vs_R848")),
  #   B18transfer = sig_DEG_table(deseq_res$B18transfer, "all", 1, 0.05, "B1-8 - 2DG_vs_Ctrl"),
  ## PL23_v_NP = sig_DEG_table(deseq_res$PL23_v_NP, "all", 1, 0.05, "AM14 PL2-3 vs B1-8 NP"),
  #   AM14MRLlpr = sig_DEG_table(deseq_res$AM14MRLlpr, "all", 1, 0.05, "AM14 MRLlpr - 2DG_vs_Ctrl")
  # )
  
  

  
  # sig_DEGs_original <- list(
  #   AM14transfer = list(
  #     PL23_2DG_v_Ctrl = sig_DEG_table(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "external_gene_name", 1, 0.05, NULL),
  #     R848_2DG_v_Ctrl = sig_DEG_table(deseq_res$AM14transfer$R848_2DG_v_Ctrl, "external_gene_name", 1, 0.05, NULL),
  #     PL23_vs_R848 = sig_DEG_table(deseq_res$AM14transfer$PL23_vs_R848, "external_gene_name", 1, 0.05, NULL)),
  #   B18transfer = sig_DEG_table(deseq_res$B18transfer, "external_gene_name", 1, 0.05, NULL),
  ##   PL23_v_NP = sig_DEG_table(deseq_res$PL23_v_NP, "external_gene_name", 1, 0.05, NULL),
  #   AM14MRLlpr = sig_DEG_table(deseq_res$AM14MRLlpr, "external_gene_name", 1, 0.05, NULL)
  # )
  # 
  # sig_DEGs_including_RNAs <- list(
  #   AM14transfer = list(
  #     PL23_2DG_v_Ctrl = sig_DEG_table(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "external_gene_name", 1, 0.05, NULL),
  #     R848_2DG_v_Ctrl = sig_DEG_table(deseq_res$AM14transfer$R848_2DG_v_Ctrl, "external_gene_name", 1, 0.05, NULL),
  #     PL23_vs_R848 = sig_DEG_table(deseq_res$AM14transfer$PL23_vs_R848, "external_gene_name", 1, 0.05, NULL)),
  #   B18transfer = sig_DEG_table(deseq_res$B18transfer, "external_gene_name", 1, 0.05, NULL),
  ##   PL23_v_NP = sig_DEG_table(deseq_res$PL23_v_NP, "external_gene_name", 1, 0.05, NULL),
  #   AM14MRLlpr = sig_DEG_table(deseq_res$AM14MRLlpr, "external_gene_name", 1, 0.05, NULL)
  # )
  # 
  # full_summary  
  # summary(comparedf(sig_DEGs_original$PL23_v_NP$up,
  #           sig_DEGs_including_RNAs$PL23_v_NP$up))
  # anti_join(sig_DEGs_original$PL23_v_NP$up,
  #           sig_DEGs_including_RNAs$PL23_v_NP$up,
  #           by = "external_gene_name")
  # anti_join(sig_DEGs_including_RNAs$PL23_v_NP$up,
  #           sig_DEGs_original$PL23_v_NP$up,
  #           by = "external_gene_name")
  # 
