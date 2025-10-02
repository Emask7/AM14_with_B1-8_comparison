setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/")
getwd()

# Import raw counts from Rosalind output ---------------------------------------
  cts_AM14trans <- import_Rosalind_data("raw_data/AM14_Adoptive_Transfer_rawCountsWithAnnotations.txt", ipa_format = TRUE)
  cts_B18trans <- import_Rosalind_data("raw_data/B1-8_Adoptive_Transfer_rawCountsWithAnnotations.txt", ipa_format = TRUE)
  cts_AM14MRLlpr <- import_Rosalind_data("raw_data/AM14_MRLlpr_rawCountsWithAnnotations.txt", ipa_format = TRUE)

  # Save a table of gene IDs in case they need to be converted later -----------
    gene_IDs <- list(AM14trans = cts_AM14trans[, c(1:6)],
                     B18trans = cts_B18trans[, c(1:6)],
                     AM14MRLlpr = cts_AM14MRLlpr[, c(1:6)])
        
  # Remove extra columns from cts ----------------------------------------------
    cts_AM14trans <- cts_AM14trans[, c(7:26)]
    cts_AM14trans_PL23 <- cts_AM14trans[, c(1:3, 6:8)]
    cts_AM14trans_R848 <- cts_AM14trans[, c(13:15, 18:20)]
    cts_B18trans <- cts_B18trans[, c(7:15)]
    cts_AM14MRLlpr <- cts_AM14MRLlpr[, c(9:14, 17:23)]
  
    head(cts_AM14trans_PL23)
    head(cts_AM14trans_R848)
    head(cts_B18trans)
    head(cts_AM14MRLlpr)

# Set up experimental factors --------------------------------------------------
  coldata <- read.csv("raw_data/sample_info.csv")
  coldata <- coldata[, c(1:3, 5:6, 10)]
  coldata
  
  cd_AM14trans_PL23 <- filter(coldata, Cohort == "A1" & grepl("PL23", Treatment))
  cd_AM14trans_PL23$Cohort <- factor(cd_AM14trans_PL23$Cohort)
  cd_AM14trans_PL23$Treatment <- factor(cd_AM14trans_PL23$Treatment)
  cd_AM14trans_PL23
  
  cd_AM14trans_R848 <- filter(coldata, Cohort == "A2" & grepl("R848", Treatment))
  cd_AM14trans_R848$Cohort <- factor(cd_AM14trans_R848$Cohort)
  cd_AM14trans_R848$Treatment <- factor(cd_AM14trans_R848$Treatment)
  cd_AM14trans_R848
  
  cd_B18trans <- coldata[21:29, ]
  cd_B18trans$Cohort <- factor(cd_B18trans$Cohort)
  cd_B18trans$Treatment <- factor(cd_B18trans$Treatment)
  cd_B18trans
  
  cd_AM14MRLlpr <- filter(coldata, Cohort == "M2")
  cd_AM14MRLlpr$Cohort <- factor(cd_AM14MRLlpr$Cohort)
  cd_AM14MRLlpr$Treatment <- factor(cd_AM14MRLlpr$Treatment)
  cd_AM14MRLlpr

  # Make sure the columns of cts and rows of coldata are in the same order -----
    for (x in 1:ncol(cts_AM14trans_PL23)) {
      if(colnames(cts_AM14trans_PL23)[x] == cd_AM14trans_PL23$Sample[x]) print(c(x, "true"))
      else print(c(x, "false"))
    }
    for (x in 1:ncol(cts_AM14trans_R848)) {
      if(colnames(cts_AM14trans_R848)[x] == cd_AM14trans_R848$Sample[x]) print(c(x, "true"))
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

# Set up DESeq Data Set --------------------------------------------------------
  dds_AM14trans_PL23 <- DESeqDataSetFromMatrix(countData = cts_AM14trans_PL23, 
                                                colData = cd_AM14trans_PL23,
                                                design = ~ Treatment)
  dds_AM14trans_PL23
  keep <- rowSums(counts(dds_AM14trans_PL23) >= 10) >= 3
  dds_AM14trans_PL23 <- dds_AM14trans_PL23[keep,]
  dds_AM14trans_PL23 <- DESeq(dds_AM14trans_PL23)
  resultsNames(dds_AM14trans_PL23)
  rm(keep)
  
  dds_AM14trans_R848 <- DESeqDataSetFromMatrix(countData = cts_AM14trans_R848, 
                                               colData = cd_AM14trans_R848,
                                               design = ~ Treatment)
  dds_AM14trans_R848
  keep <- rowSums(counts(dds_AM14trans_R848) >= 10) >= 3
  dds_AM14trans_R848 <- dds_AM14trans_R848[keep,]
  dds_AM14trans_R848 <- DESeq(dds_AM14trans_R848)
  resultsNames(dds_AM14trans_R848)
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
  
  dds_AM14MRLlpr <- DESeqDataSetFromMatrix(countData = cts_AM14MRLlpr,
                                               colData = cd_AM14MRLlpr,
                                               design = ~Treatment)
  dds_AM14MRLlpr
  keep <- rowSums(counts(dds_AM14MRLlpr) >= 10) >= 6
  dds_AM14MRLlpr <- dds_AM14MRLlpr[keep,]
  dds_AM14MRLlpr <- DESeq(dds_AM14MRLlpr)
  resultsNames(dds_AM14MRLlpr)
  rm(keep)
  
    
# Run QC steps -----------------------------------------------------------------
  setwd("./analysis_with_cohorts_and_treatment_separated/")
  getwd()
  
  QC_heatmaps(dds_AM14trans_PL23, "AM14_Adoptive_Transfer_PL23", "AM14 Adoptive Transfer:\nPL2-3 + 2DG vs PL2-3", plot_cohorts = FALSE)
  QC_heatmaps(dds_AM14trans_R848, "AM14_Adoptive_Transfer_R848", "AM14 Adoptive Transfer:\nR848 + 2DG vs R848", plot_cohorts = FALSE)
  QC_heatmaps(dds_B18trans, "B18_Adoptive_Transfer", "B1-8 Adoptive Transfer\nNP + 2DG vs NP", plot_cohorts = FALSE)
  QC_heatmaps(dds_AM14MRLlpr, "AM14_MRLlpr", "AM14 MRL/lpr mice\n2DG vs Control", plot_cohorts = FALSE)

  sample_to_sample_plot(dds_AM14trans_PL23, "AM14_Adoptive_Transfer_PL23", "AM14 Adoptive Transfer (PL2-3 + 2DG vs PL2-3)")
  sample_to_sample_plot(dds_AM14trans_R848, "AM14_Adoptive_Transfer_R848", "AM14 Adoptive Transfer (R848 + 2DG vs R848)")
  sample_to_sample_plot(dds_B18trans, "B18_Adoptive_Transfer", "B1-8 Adoptive Transfer")
  sample_to_sample_plot(dds_AM14MRLlpr, "AM14_MRLlpr", "AM14 MRL/lpr +/- 2DG")
  
  QC_PCAplot(dds_AM14trans_PL23, "AM14_Adoptive_Transfer_PL23", "PL2-3 + 2DG vs PL2-3",
             batch_effect = NULL, draw_ellipse = FALSE)
  dev.off()
  
  QC_PCAplot(dds_AM14trans_R848, "AM14_Adoptive_Transfer_R848", "R848 + 2DG vs R848",
             batch_effect = NULL, draw_ellipse = FALSE)
  dev.off()

  QC_PCAplot(dds_B18trans, "B1-8_Adoptive_Transfer", "NP + 2DG vs NP",
             batch_effect = NULL, draw_ellipse = FALSE)
  dev.off()

  QC_PCAplot(dds_AM14MRLlpr, "AM14_MRLlpr", "AM14 MRL/lpr mice: 2DG vs Control",
             batch_effect = NULL, draw_ellipse = FALSE)
  dev.off()

# Get DESeq2 results -----------------------------------------------------------
  deseq_res <- list(
    AM14transfer_PL23 = results_wrapper(dds_AM14trans_PL23, c("Treatment", "PL23_2DG", "PL23"), gene_IDs$AM14trans),
    AM14transfer_R848 = results_wrapper(dds_AM14trans_R848, c("Treatment", "R848_2DG", "R848"), gene_IDs$AM14trans),
    B18transfer = results_wrapper(dds_B18trans, c("Treatment", "NP_2DG", "NP"), gene_IDs$B18trans),
    AM14MRLlpr = results_wrapper(dds_AM14MRLlpr, c("Treatment", "2DG", "Control"), gene_IDs$AM14MRLlpr)
  )
  
  head(deseq_res$AM14transfer_PL23)
  head(deseq_res$AM14transfer_R848)
  head(deseq_res$B18transfer)
  head(deseq_res$AM14MRLlpr)
  
  # DEG Heatmaps ---------------------------------------------------------------
  # # DEG_heatmap <- function(dds, dds_results, factors, color_list, lfc_cutoff, padj_cutoff, count_type, dds_cols, plot_title, filename, h, w){
  # 
  # DEG_heatmap(dds_AM14trans, deseq_res$AM14transfer$PL23_2DG_v_Ctrl,
  #             c("Treatment", "Cohort"), 
  #             list(Treatment = c(PL23 = "#0f85a0", PL23_2DG = "#dd4124"),
  #                  Cohort = c(A1 = "#00496f", A2 = "#edd746")),
  #             1, 0.05, "norm_counts", c(1:10),
  #             "AM14 Adoptive Transfer\n(PL2-3 +/- 2DG Significant DEGs)",
  #             "AM14 transfer - PL2-3 only - DESeq2_norm")
  
# DEG Heatmaps using default DESeq2 normalization ------------------------------  
    
  DEG_heatmap(dds_AM14trans_PL23, deseq_res$AM14transfer_PL23, 
              c("Treatment"),
              list(Treatment = c(PL23 = "#0f85a0", PL23_2DG = "#dd4124")),
              1, 0.05, "norm_counts", NULL,
              "AM14 Adoptive Transfer: PL2-3 + 2DG vs PL2-3",
              "AM14 transfer PL2-3 - padj 0.05")

  DEG_heatmap(dds_AM14trans_PL23, deseq_res$AM14transfer_PL23, 
              c("Treatment"),
              list(Treatment = c(PL23 = "#0f85a0", PL23_2DG = "#dd4124")),
              1, 0.01, "norm_counts", NULL,
              "AM14 Adoptive Transfer: PL2-3 + 2DG vs PL2-3",
              "AM14 transfer PL2-3 - padj 0.01")
  
  DEG_heatmap(dds_AM14trans_R848, deseq_res$AM14transfer_R848, 
              c("Treatment"),
              list(Treatment = c(R848 = "#0f85a0", R848_2DG = "#dd4124")),
              1, 0.05, "norm_counts", NULL,
              "AM14 Adoptive Transfer: R848 + 2DG vs R848",
              "AM14 transfer R848 - padj 0.05")
  
  DEG_heatmap(dds_AM14trans_R848, deseq_res$AM14transfer_R848, 
              c("Treatment"),
              list(Treatment = c(R848 = "#0f85a0", R848_2DG = "#dd4124")),
              1, 0.01, "norm_counts", NULL,
              "AM14 Adoptive Transfer: R848 + 2DG vs R848",
              "AM14 transfer R848 - padj 0.01")
  

  DEG_heatmap(dds_B18trans, deseq_res$B18transfer, 
              c("Treatment"),
              list(Treatment = c(NP = "#0f85a0", NP_2DG = "#dd4124")),
              1, 0.05, "norm_counts", NULL,
              "B1-8 Adoptive Transfer: NP + 2DG vs NP",
              "B1-8 transfer - padj 0.05", h = 1000)
  
  DEG_heatmap(dds_B18trans, deseq_res$B18transfer, 
              c("Treatment"),
              list(Treatment = c(NP = "#0f85a0", NP_2DG = "#dd4124")),
              1, 0.01, "norm_counts", NULL,
              "B1-8 Adoptive Transfer: NP + 2DG vs NP",
              "B1-8 transfer - padj 0.01", h = 1000)
  
  
  DEG_heatmap(dds_AM14MRLlpr, deseq_res$AM14MRLlpr, 
              c("Treatment"),
              list(Treatment = c(Control = "#0f85a0", "2DG" = "#dd4124")),
              1, 0.05, "norm_counts", NULL,
              "AM14 MRL/lpr Mice: 2DG vs Control",
              "AM14 MRLlpr - padj 0.05")
  
  DEG_heatmap(dds_AM14MRLlpr, deseq_res$AM14MRLlpr, 
              c("Treatment"),
              list(Treatment = c(Control = "#0f85a0", "2DG" = "#dd4124")),
              1, 0.01, "norm_counts", NULL,
              "AM14 MRL/lpr Mice: 2DG vs Control",
              "AM14 MRLlpr - padj 0.01")
  
# # DEG Heatmaps using Rlog normalization ----------------------------------------  
#   DEG_heatmap(dds_AM14trans_PL23, deseq_res$AM14transfer_PL23, 
#               c("Treatment"),
#               list(Treatment = c(PL23 = "#0f85a0", PL23_2DG = "#dd4124")),
#               1, 0.05, "rlog", NULL,
#               "AM14 Adoptive Transfer: PL2-3 + 2DG vs PL2-3",
#               "AM14 transfer PL2-3 - padj 0.05 - Rlog_norm")
#   
#   DEG_heatmap(dds_AM14trans_PL23, deseq_res$AM14transfer_PL23, 
#               c("Treatment"),
#               list(Treatment = c(PL23 = "#0f85a0", PL23_2DG = "#dd4124")),
#               1, 0.01, "rlog", NULL,
#               "AM14 Adoptive Transfer: PL2-3 + 2DG vs PL2-3",
#               "AM14 transfer PL2-3 - padj 0.01 - Rlog_norm")
#   
#   DEG_heatmap(dds_AM14trans_R848, deseq_res$AM14transfer_R848, 
#               c("Treatment"),
#               list(Treatment = c(R848 = "#0f85a0", R848_2DG = "#dd4124")),
#               1, 0.05, "rlog", NULL,
#               "AM14 Adoptive Transfer: R848 + 2DG vs R848",
#               "AM14 transfer R848 - padj 0.05 - Rlog_norm")
#   
#   DEG_heatmap(dds_AM14trans_R848, deseq_res$AM14transfer_R848, 
#               c("Treatment"),
#               list(Treatment = c(R848 = "#0f85a0", R848_2DG = "#dd4124")),
#               1, 0.01, "rlog", NULL,
#               "AM14 Adoptive Transfer: R848 + 2DG vs R848",
#               "AM14 transfer R848 - padj 0.01 - Rlog_norm")
#   
#   DEG_heatmap(dds_B18trans, deseq_res$B18transfer, 
#               c("Treatment"),
#               list(Treatment = c(NP = "#0f85a0", NP_2DG = "#dd4124")),
#               1, 0.05, "rlog", NULL,
#               "B1-8 Adoptive Transfer: NP + 2DG vs NP",
#               "B1-8 transfer - padj 0.05 - RLog Normalization")
#   
#   DEG_heatmap(dds_B18trans, deseq_res$B18transfer, 
#               c("Treatment"),
#               list(Treatment = c(NP = "#0f85a0", NP_2DG = "#dd4124")),
#               1, 0.01, "rlog", NULL,
#               "B1-8 Adoptive Transfer: NP + 2DG vs NP",
#               "B1-8 transfer - padj 0.01 - RLog Normalization")
#   
#   
#   DEG_heatmap(dds_AM14MRLlpr, deseq_res$AM14MRLlpr, 
#               c("Treatment"),
#               list(Treatment = c(Control = "#0f85a0", "2DG" = "#dd4124")),
#               1, 0.05, "rlog", NULL,
#               "AM14 MRL/lpr Mice: 2DG vs Control",
#               "AM14 MRLlpr - padj 0.05 - RLog Normalization")
#   
#   DEG_heatmap(dds_AM14MRLlpr, deseq_res$AM14MRLlpr, 
#               c("Treatment"),
#               list(Treatment = c(Control = "#0f85a0", "2DG" = "#dd4124")),
#               1, 0.01, "rlog", NULL,
#               "AM14 MRL/lpr Mice: 2DG vs Control",
#               "AM14 MRLlpr - padj 0.01 - RLog Normalization")
  
  
  
    

# Make a data.frame that summarizes the numbers of DEGs detected -------------
  full_summary <- join_all(list(summary_wrapper(deseq_res$AM14transfer_PL23, "PL2-3+2DG vs PL2-3", "AM14 Adoptive Transfer: Cohort 1", 0.05, 1),
                                summary_wrapper(deseq_res$AM14transfer_R848, "R848+2DG vs R848", "AM14 Adoptive Transfer: Cohort 2", 0.05, 1),
                                summary_wrapper(deseq_res$B18transfer, "NP+2DG vs NP", "B1-8 Adoptive Transfer", 0.05, 1),
                                summary_wrapper(deseq_res$AM14MRLlpr, "2DG vs Control", "AM14 MRL/lpr: Cohort 2", 0.05, 1),
                                summary_wrapper(deseq_res$AM14transfer_PL23, "PL2-3+2DG vs PL2-3", "AM14 Adoptive Transfer: Cohort 1", 0.01, 1),
                                summary_wrapper(deseq_res$AM14transfer_R848, "R848+2DG vs R848", "AM14 Adoptive Transfer: Cohort 2", 0.01, 1),
                                summary_wrapper(deseq_res$B18transfer, "NP+2DG vs NP", "B1-8 Adoptive Transfer", 0.01, 1),
                                summary_wrapper(deseq_res$AM14MRLlpr, "2DG vs Control", "AM14 MRL/lpr: Cohort 2", 0.01, 1)),
                           type = "full")
  full_summary

# Save results to an Excel file ----------------------------------------------
  wb <- createWorkbook("Output/DESeq2_results.xlsx")
  
  addWorksheet(wb, "DEG Summaries")
  writeData(wb, "DEG Summaries", full_summary)
  
  addWorksheet(wb, "AM14 Co1 - PL23_2DG_vs_Ctrl")
  addWorksheet(wb, "B1-8 - 2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 Co2 - R848_2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 MRLlpr Co2 - 2DG_vs_Ctrl")
  
  writeData(wb, "AM14 Co1 - PL23_2DG_vs_Ctrl", deseq_res$AM14transfer_PL23)
  writeData(wb, "B1-8 - 2DG_vs_Ctrl", deseq_res$B18transfer)
  writeData(wb, "AM14 Co2 - R848_2DG_vs_Ctrl", deseq_res$AM14transfer_R848)
  writeData(wb, "AM14 MRLlpr Co2 - 2DG_vs_Ctrl", deseq_res$AM14MRLlpr)
  
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
