# Import raw counts from Rosalind output ---------------------------------------
  cts_AM14trans <- import_Rosalind_data("raw_data/AM14_Adoptive_Transfer_rawCountsWithAnnotations.txt")
  cts_B18trans <- import_Rosalind_data("raw_data/B1-8_Adoptive_Transfer_rawCountsWithAnnotations.txt")
  cts_AM14MRLlpr <- import_Rosalind_data("raw_data/AM14_MRLlpr_rawCountsWithAnnotations.txt")
  
  cts_full <- join_all(list(cts_AM14trans, cts_B18trans, cts_AM14MRLlpr), type = "full")
  
  # Save a table of gene IDs in case they need to be converted later -----------
    gene_IDs <- cts_full[, 1:4]
    head(gene_IDs)
    
  rownames(cts_full) <- cts_full$ensembl_gene_id
  cts_full <- cts_full[, 5:50]
  head(cts_full)
  
  rm(cts_AM14trans, cts_B18trans, cts_AM14MRLlpr)
  

# Set up experimental factors --------------------------------------------------
  coldata <- read.csv("raw_data/sample_info.csv")
  coldata <- coldata[, c(1:3, 5:6)]
  coldata$Cohort <- factor(coldata$Cohort)
  coldata$Treatment <- factor(coldata$Treatment)
  coldata

  cohort <- factor(c(rep("A1", 3), rep("A2", 2), rep("A1", 3), rep("A2", 2), 
                   rep("A1", 2), rep("A2", 3), rep("A1", 2), rep("A2", 3), 
                   rep("B1", 9), 
                   rep("M1", 2), rep("M2", 6), rep("M1", 2), rep("M2", 7)))
  cohort
  
  treatment <- factor(c(rep("PL23", 5), rep("PL23_2DG", 5), rep("R848", 5), rep("R848_2DG", 5),
                        rep("NP_2DG", 5), rep("NP", 4), rep("Control", 8), rep("2DG", 9)))
  treatment
  
  d <- DataFrame(cohort, treatment)[,]
  as.data.frame(d)
  m1 <- model.matrix(~ treatment*cohort, d)
  colnames(m1)
  unname(m1)
  all.zero <- apply(m1, 2, function(x) all(x==0))
  all.zero  
  idx <- which(all.zero)
  m1 <- m1[,-idx]
  unname(m1)
  rownames(m1) <- colnames(cts_full)
  m1
  
  # # Make sure the columns of cts and rows of coldata are in the same order -----
  #   for (x in 1:ncol(cts_full)) {
  #     if(colnames(cts_full)[x] == coldata$Sample[x]) print(c(x, "true")) 
  #     else print(c(x, "false"))
  #   }
  #   rm(x)

# Set up DESeq Data Set --------------------------------------------------------
  dds <- DESeqDataSetFromMatrix(countData = cts_full,
                                colData = coldata,
                                design = ~Treatment)
  dds
  # keep <- rowSums(counts(dds) >= 10) >= 5
  # dds <- dds[keep,]
  # dds <- DESeq(dds, full = m1, reduced = model.matrix(~ treatment*cohort, d), test = "LRT")
  dds <- DESeq(dds)
  
  resultsNames(dds)
  # rm(keep)

# Run QC steps -----------------------------------------------------------------
  QC_PCAplot(dds, "AM14_Adoptive_Transfer-Before_Batch_Correction", 
             # "AM14 Adoptive Transfer\n(Before Batch Correction)", 
             NULL,
             batch_effect = FALSE)               
  dev.off()

# Differential Expression Analysis ---------------------------------------------
  resultsNames(dds)

# Get DESeq2 results -----------------------------------------------------------
  deseq_res <- list(
    AM14transfer = list(
      PL23_2DG_v_Ctrl = results_wrapper(dds, c("Treatment", "PL23_2DG", "PL23"), gene_IDs),
      R848_2DG_v_Ctrl = results_wrapper(dds, c("Treatment", "R848_2DG", "R848"), gene_IDs),
      PL23_vs_R848 = results_wrapper(dds, c("Treatment", "PL23", "R848"), gene_IDs)),
    B18transfer = results_wrapper(dds, c("Treatment", "NP_2DG", "NP"), gene_IDs),
    AM14MRLlpr = results_wrapper(dds, c("Treatment", "2DG", "Control"), gene_IDs),
    AM14_PL23_vs_B18_NP = results_wrapper(dds, c("Treatment", "PL23", "NP"), gene_IDs)
  )
  
  head(deseq_res$AM14transfer$PL23_2DG_v_Ctrl)
  head(deseq_res$B18transfer)
  head(deseq_res$AM14MRLlpr)
  head(deseq_res$AM14_PL23_vs_B18_NP)
  
  # Make a data.frame that summarizes the numbers of DEGs detected -------------
  full_summary <- join_all(list(summary_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "PL2-3+2DG vs PL2-3", "AM14 Adoptive Transfer"),
                                summary_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl, "R848+2DG vs R848", "AM14 Adoptive Transfer"),
                                summary_wrapper(deseq_res$AM14transfer$PL23_vs_R848, "PL2-3 vs R848", "AM14 Adoptive Transfer"),
                                summary_wrapper(deseq_res$B18transfer, "NP+2DG vs NP", "B1-8 Adoptive Transfer"),
                                summary_wrapper(deseq_res$AM14_PL23_vs_B18_NP, "PL2-3 vs NP", "AM14 and B1-8 Adoptive Transfer"),
                                summary_wrapper(deseq_res$AM14MRLlpr, "2DG vs Control", "AM14 MRL/lpr")), 
                           type = "full")
  full_summary 
  
  # Save results to an Excel file ----------------------------------------------
  wb <- createWorkbook("Output/DESeq2_results.xlsx")
  
  addWorksheet(wb, "DEG Summaries")
  addWorksheet(wb, "AM14 - PL23_2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 - R848_2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 - PL23_vs_R848")
  addWorksheet(wb, "B1-8 - 2DG_vs_Ctrl")
  addWorksheet(wb, "AM14 MRLlpr - 2DG_vs_Ctrl")
  
  writeData(wb, "DEG Summaries", full_summary)
  writeData(wb, "AM14 - PL23_2DG_vs_Ctrl", deseq_res$AM14transfer$PL23_2DG_v_Ctrl)
  writeData(wb, "AM14 - R848_2DG_vs_Ctrl", deseq_res$AM14transfer$R848_2DG_v_Ctrl)
  writeData(wb, "AM14 - PL23_vs_R848", deseq_res$AM14transfer$PL23_vs_R848)
  writeData(wb, "B1-8 - 2DG_vs_Ctrl", deseq_res$B18transfer)
  writeData(wb, "AM14 MRLlpr - 2DG_vs_Ctrl", deseq_res$AM14MRLlpr)
  
  saveWorkbook(wb, "Output/DESeq2_results.xlsx", overwrite = TRUE)
  rm(wb)
  
  
  wb <- createWorkbook("Output/DESeq2_gene_counts.xlsx")
  
  addWorksheet(wb, "Raw_Counts - AM14 transfer")
  addWorksheet(wb, "Norm_Counts - AM14 transfer")
  addWorksheet(wb, "Raw_Counts - B1-8 transfer")
  addWorksheet(wb, "Norm_Counts - B1-8 transfer")
  addWorksheet(wb, "Raw_Counts - AM14 MRLlpr")
  addWorksheet(wb, "Norm_Counts - AM14 MRLlpr")
  
  writeData(wb, "Raw_Counts - AM14 transfer", counts(dds, normalized = FALSE), rowNames = TRUE)
  writeData(wb, "Norm_Counts - AM14 transfer", counts(dds, normalized = TRUE), rowNames = TRUE)
  writeData(wb, "Raw_Counts - B1-8 transfer", counts(dds, normalized = FALSE), rowNames = TRUE)
  writeData(wb, "Norm_Counts - B1-8 transfer", counts(dds, normalized = TRUE), rowNames = TRUE)
  writeData(wb, "Raw_Counts - AM14 MRLlpr", counts(dds, normalized = FALSE), rowNames = TRUE)
  writeData(wb, "Norm_Counts - AM14 MRLlpr", counts(dds, normalized = TRUE), rowNames = TRUE)
  
  saveWorkbook(wb, "Output/DESeq2_gene_counts.xlsx", overwrite = TRUE)
  rm(wb)
  
  # # Make DEG lists and export to CSV files -------------------------------------
  #   write_DEG_CSV(res_PL23, 1, "PL2-3_2DG_vs_Control")
  #   write_DEG_CSV(res_R848, 1, "R848_2DG_vs_Control")
  #   write_DEG_CSV(res_PL23vR848, 1, "PL2-3_Ctrl_vs_R848_Ctrl")
  #   write_DEG_CSV(res_2DG_PL23vR848, 1, "PL2-3_2DG_vs_R848_2DG")