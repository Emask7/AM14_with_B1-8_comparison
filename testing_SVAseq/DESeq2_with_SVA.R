BiocManager::install("sva")
library(sva)


# Import raw counts from Rosalind output ---------------------------------------
  cts_AM14trans <- import_Rosalind_data("raw_data/AM14_Adoptive_Transfer_rawCountsWithAnnotations.txt")
  cts_B18trans <- import_Rosalind_data("raw_data/B1-8_Adoptive_Transfer_rawCountsWithAnnotations.txt")
  # cts_AM14MRLlpr <- import_Rosalind_data("raw_data/AM14_MRLlpr_rawCountsWithAnnotations.txt")
  
  # cts_full <- join_all(list(cts_AM14trans, cts_B18trans, cts_AM14MRLlpr), type = "full")
  cts_full <- join_all(list(cts_AM14trans, cts_B18trans), type = "full")
  
  # Save a table of gene IDs in case they need to be converted later -----------
    gene_IDs <- cts_full[, 1:4]
    head(gene_IDs)
  
  rownames(cts_full) <- cts_full$ensembl_gene_id
  ncol(cts_full)
  cts_full <- cts_full[, 5:33]
  head(cts_full)
  
  # rm(cts_AM14trans, cts_B18trans, cts_AM14MRLlpr)
  rm(cts_AM14trans, cts_B18trans)
  
# Set up experimental factors --------------------------------------------------
  coldata <- read.csv("raw_data/sample_info.csv")
  # coldata <- coldata[, c(1:3, 5:6)]
  coldata <- coldata[1:29, c(1:3, 5:6)]
  coldata$Cohort <- factor(coldata$Cohort)
  coldata$Treatment <- factor(coldata$Treatment)
  coldata
  

# Set up DESeq Data Set with only treatment as the design ----------------------
  dds <- DESeqDataSetFromMatrix(countData = cts_full,
                                colData = coldata,
                                design = ~Treatment)
  dds
  dds <- DESeq(dds)
  resultsNames(dds)

# Run QC steps -----------------------------------------------------------------
  QC_PCAplot(dds, "Full-counts", NULL, batch_effect = FALSE)               
  dev.off()
  
# Use SVA ----------------------------------------------------------------------
  sva_dat  <- counts(dds, normalized = TRUE)
  idx  <- rowMeans(sva_dat) > 1
  sva_dat  <- sva_dat[idx, ]
  mod  <- model.matrix(~ Treatment, colData(dds))
  mod0 <- model.matrix(~ 1, colData(dds))
  svseq <- svaseq(sva_dat, mod, mod0, n.sv = NULL)
  svseq$sv
  
  par(mfrow = c(2, 2), mar = c(3,5,3,1))
  for (i in 1:4) {
    stripchart(svseq$sv[, i] ~ dds$Cohort, vertical = TRUE, main = paste0("SV", i))
    abline(h = 0)
  }
  dev.off()
  
  ddssva <- dds
  ddssva$SV1 <- svseq$sv[,1]
  ddssva$SV2 <- svseq$sv[,2]
  ddssva$SV3 <- svseq$sv[,3]
  design(ddssva) <- ~ SV1 + SV2 + SV3 + Treatment
  ddssva <- DESeq(ddssva)
  resultsNames(ddssva)
  
  
  
# Differential Expression Analysis ---------------------------------------------
# Get DESeq2 results -----------------------------------------------------------
  sva_deseq_res <- list(
    AM14transfer = list(
      PL23_2DG_v_Ctrl = results_wrapper(dds, c("Treatment", "PL23_2DG", "PL23"), gene_IDs),
      R848_2DG_v_Ctrl = results_wrapper(dds, c("Treatment", "R848_2DG", "R848"), gene_IDs),
      PL23_vs_R848 = results_wrapper(dds, c("Treatment", "PL23", "R848"), gene_IDs)),
    B18transfer = results_wrapper(dds, c("Treatment", "NP_2DG", "NP"), gene_IDs),
    PL23_v_NP = results_wrapper(dds, c("Treatment", "PL23", "NP"), gene_IDs)
    # AM14MRLlpr = results_wrapper(dds, c("Treatment", "2DG", "Control"), gene_IDs)
  )
  
  head(sva_deseq_res$AM14transfer$PL23_2DG_v_Ctrl)
  head(sva_deseq_res$B18transfer)
  head(sva_deseq_res$AM14_PL23_vs_B18_NP)
  # head(sva_deseq_res$AM14MRLlpr)
  
  # Make a data.frame that summarizes the numbers of DEGs detected -------------
  full_summary <- join_all(list(summary_wrapper(sva_deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "PL2-3+2DG vs PL2-3", "AM14 Adoptive Transfer"),
                                summary_wrapper(sva_deseq_res$AM14transfer$R848_2DG_v_Ctrl, "R848+2DG vs R848", "AM14 Adoptive Transfer"),
                                summary_wrapper(sva_deseq_res$AM14transfer$PL23_vs_R848, "PL2-3 vs R848", "AM14 Adoptive Transfer"),
                                summary_wrapper(sva_deseq_res$B18transfer, "NP+2DG vs NP", "B1-8 Adoptive Transfer"),
                                summary_wrapper(sva_deseq_res$PL23_v_NP, "PL2-3 vs NP", "AM14 and B1-8 Adoptive Transfer")),
                                # summary_wrapper(sva_deseq_res$AM14MRLlpr, "2DG vs Control", "AM14 MRL/lpr") 
                           type = "full")
  full_summary 
  
  # Save results to an Excel file ----------------------------------------------
  wb <- createWorkbook("testing_SVAseq/DESeq2_SVA_results.xlsx")
  
  addWorksheet(wb, "DEG Summaries")
  addWorksheet(wb, "AM14 - PL23_2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 - R848_2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 - PL23_vs_R848")
  addWorksheet(wb, "B1-8 - 2DG_vs_Ctrl")
  # addWorksheet(wb, "AM14 MRLlpr - 2DG_vs_Ctrl")
  
  writeData(wb, "DEG Summaries", full_summary)
  writeData(wb, "AM14 - PL23_2DG_vs_Ctrl", sva_deseq_res$AM14transfer$PL23_2DG_v_Ctrl)
  writeData(wb, "AM14 - R848_2DG_vs_Ctrl", sva_deseq_res$AM14transfer$R848_2DG_v_Ctrl)
  writeData(wb, "AM14 - PL23_vs_R848", sva_deseq_res$AM14transfer$PL23_vs_R848)
  writeData(wb, "B1-8 - 2DG_vs_Ctrl", sva_deseq_res$B18transfer)
  # writeData(wb, "AM14 MRLlpr - 2DG_vs_Ctrl", sva_deseq_res$AM14MRLlpr)
  
  saveWorkbook(wb, "testing_SVAseq/DESeq2_SVA_results.xlsx", overwrite = TRUE)
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
  # writeData(wb, "Raw_Counts - AM14 transfer", counts(dds, normalized = FALSE), rowNames = TRUE)
  # writeData(wb, "Norm_Counts - AM14 transfer", counts(dds, normalized = TRUE), rowNames = TRUE)
  # writeData(wb, "Raw_Counts - B1-8 transfer", counts(dds, normalized = FALSE), rowNames = TRUE)
  # writeData(wb, "Norm_Counts - B1-8 transfer", counts(dds, normalized = TRUE), rowNames = TRUE)
  # writeData(wb, "Raw_Counts - AM14 MRLlpr", counts(dds, normalized = FALSE), rowNames = TRUE)
  # writeData(wb, "Norm_Counts - AM14 MRLlpr", counts(dds, normalized = TRUE), rowNames = TRUE)
  # 
  # saveWorkbook(wb, "Output/DESeq2_gene_counts.xlsx", overwrite = TRUE)
  # rm(wb)
  
  # # Make DEG lists and export to CSV files -------------------------------------
  #   write_DEG_CSV(res_PL23, 1, "PL2-3_2DG_vs_Control")
  #   write_DEG_CSV(res_R848, 1, "R848_2DG_vs_Control")
  #   write_DEG_CSV(res_PL23vR848, 1, "PL2-3_Ctrl_vs_R848_Ctrl")
  #   write_DEG_CSV(res_2DG_PL23vR848, 1, "PL2-3_2DG_vs_R848_2DG")