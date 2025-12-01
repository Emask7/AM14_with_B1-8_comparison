setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/")
getwd()

# Import raw counts from Rosalind output ---------------------------------------
  cts_AM14MRLlpr_noout <- import_Rosalind_data("raw_data/AM14_MRLlpr_rawCountsWithAnnotations.txt", ipa_format = TRUE)

  # Save a table of gene IDs in case they need to be converted later -----------
    gene_IDs <-  cts_AM14MRLlpr_noout[, c(1:6)]
    head(gene_IDs)
        
  # Remove extra columns and outlier from cts ----------------------------------
    cts_AM14MRLlpr_noout <- cts_AM14MRLlpr_noout[, c(9:14, 17:22)]
    head(cts_AM14MRLlpr_noout)

# Set up experimental factors --------------------------------------------------
  coldata <- read.csv("raw_data/sample_info.csv")
  coldata <- coldata[, c(1:3, 5:10)]
  coldata
  
  cd_AM14MRLlpr_noout <- filter(coldata, Cohort == "M2")
  cd_AM14MRLlpr_noout <- filter(cd_AM14MRLlpr_noout, Sample != "MRL_2DG_9")
  cd_AM14MRLlpr_noout$Cohort <- factor(cd_AM14MRLlpr_noout$Cohort)
  cd_AM14MRLlpr_noout$Treatment <- factor(cd_AM14MRLlpr_noout$Treatment)
  cd_AM14MRLlpr_noout

  # Make sure the columns of cts and rows of coldata are in the same order -----
    for (x in 1:ncol(cts_AM14MRLlpr_noout)) {
      if(colnames(cts_AM14MRLlpr_noout)[x] == cd_AM14MRLlpr_noout$Sample[x]) print(c(x, "true"))
      else print(c(x, "false"))
    }
    rm(x)

# Set up DESeq Data Set --------------------------------------------------------
  dds_AM14MRLlpr_noout <- DESeqDataSetFromMatrix(countData = cts_AM14MRLlpr_noout,
                                                 colData = cd_AM14MRLlpr_noout,
                                                 design = ~Treatment)
  dds_AM14MRLlpr_noout
  keep <- rowSums(counts(dds_AM14MRLlpr_noout) >= 10) >= 6
  dds_AM14MRLlpr_noout <- dds_AM14MRLlpr_noout[keep,]
  dds_AM14MRLlpr_noout <- DESeq(dds_AM14MRLlpr_noout)
  resultsNames(dds_AM14MRLlpr_noout)
  rm(keep)
  
    
# Run QC steps -----------------------------------------------------------------
  setwd("./AM14_MRLlpr_no_outlier/")
  getwd()

  QC_heatmaps(dds_AM14MRLlpr_noout, "AM14_MRLlpr", "AM14 MRL/lpr mice\n2DG vs Control", plot_cohorts = FALSE)

  sample_to_sample_plot(dds_AM14MRLlpr_noout, "AM14_MRLlpr", "AM14 MRL/lpr +/- 2DG")

  QC_PCAplot(dds_AM14MRLlpr_noout, "AM14_MRLlpr", "AM14 MRL/lpr mice: 2DG vs Control",
             batch_effect = NULL, draw_ellipse = FALSE)
  dev.off()

# Get DESeq2 results -----------------------------------------------------------
  deseq_res_MRLlpr_noout <- results_wrapper(dds_AM14MRLlpr_noout, c("Treatment", "2DG", "Control"), gene_IDs)
  head(deseq_res_MRLlpr_noout)
  
# DEG Heatmaps using default DESeq2 normalization ------------------------------
  DEG_heatmap(dds_AM14MRLlpr_noout, deseq_res_MRLlpr_noout,
              c("Treatment"),
              list(Treatment = c(Control = "#0f85a0", "2DG" = "#dd4124")),
              0.6, 0.05, "norm_counts", NULL,
              "AM14 MRL/lpr Mice: 2DG vs Control",
              "AM14 MRLlpr - padj 0.05")
  
  DEG_heatmap(dds_AM14MRLlpr_noout, deseq_res_MRLlpr_noout,
              c("Treatment"),
              list(Treatment = c(Control = "#0f85a0", "2DG" = "#dd4124")),
              0, 0.05, "norm_counts", NULL,
              "AM14 MRL/lpr Mice: 2DG vs Control",
              "AM14 MRLlpr - padj 0.05 - any LFC")
  
  
  # DEG_heatmap(dds_AM14MRLlpr_noout, deseq_res_MRLlpr_noout,
  #             c("Treatment"),
  #             list(Treatment = c(Control = "#0f85a0", "2DG" = "#dd4124")),
  #             1, 0.01, "norm_counts", NULL,
  #             "AM14 MRL/lpr Mice: 2DG vs Control",
  #             "AM14 MRLlpr - padj 0.01")
  
# # DEG Heatmaps using Rlog normalization ----------------------------------------  
#   DEG_heatmap(dds_AM14MRLlpr_noout, deseq_res_MRLlpr_noout, 
#               c("Treatment"),
#               list(Treatment = c(Control = "#0f85a0", "2DG" = "#dd4124")),
#               1, 0.05, "rlog", NULL,
#               "AM14 MRL/lpr Mice: 2DG vs Control",
#               "AM14 MRLlpr - padj 0.05 - RLog Normalization")
#   
#   DEG_heatmap(dds_AM14MRLlpr_noout, deseq_res_MRLlpr_noout, 
#               c("Treatment"),
#               list(Treatment = c(Control = "#0f85a0", "2DG" = "#dd4124")),
#               1, 0.01, "rlog", NULL,
#               "AM14 MRL/lpr Mice: 2DG vs Control",
#               "AM14 MRLlpr - padj 0.01 - RLog Normalization")
  
  
  

# Make a data.frame that summarizes the numbers of DEGs detected -------------
  full_summary_noout <- join_all(list(summary_wrapper(deseq_res_MRLlpr_noout, "2DG vs Control", "AM14 MRL/lpr: Cohort 2", 0.05, 1),
                                summary_wrapper(deseq_res_MRLlpr_noout, "2DG vs Control", "AM14 MRL/lpr: Cohort 2", 0.05, 0.6)),
                           type = "full")
  full_summary_noout

# Save results to an Excel file ----------------------------------------------
  # setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/AM14_MRLlpr_no_outlier/")
  # getwd()

  wb <- createWorkbook("Output/DESeq2_results - AM14 MRLlpr no outlier.xlsx")

  addWorksheet(wb, "DEG Summaries")
  writeData(wb, "DEG Summaries", full_summary_noout)

  addWorksheet(wb, "AM14 MRLlpr - no outlier")
  writeData(wb, "AM14 MRLlpr - no outlier", deseq_res_MRLlpr_noout)

  saveWorkbook(wb, "Output/DESeq2_results - AM14 MRLlpr no outlier.xlsx", overwrite = TRUE)
  rm(wb)

  wb <- createWorkbook("Output/DESeq2_normalized_gene_counts - AM14 MRLlpr no outlier.xlsx")

  addWorksheet(wb, "MRLlpr")
  writeData(wb, "MRLlpr", counts(dds_AM14MRLlpr_noout, normalized = TRUE), rowNames = TRUE)

  saveWorkbook(wb, "Output/DESeq2_normalized_gene_counts - AM14 MRLlpr no outlier.xlsx", overwrite = TRUE)
  rm(wb)
