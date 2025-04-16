# Import raw counts from Rosalind output ---------------------------------------
  cts_AM14trans <- import_Rosalind_data("raw_data/AM14_Adoptive_Transfer_rawCountsWithAnnotations.txt")
  cts_B18trans <- import_Rosalind_data("raw_data/B1-8_Adoptive_Transfer_rawCountsWithAnnotations.txt")
  cts_AM14MRLlpr <- import_Rosalind_data("raw_data/AM14_MRLlpr_rawCountsWithAnnotations.txt")

  # Save a table of gene IDs in case they need to be converted later -----------
    gene_IDs <- list(AM14trans = cts_AM14trans[, c(1:4)],
                     B18trans = cts_B18trans[, c(1:4)],
                     AM14MRLlpr = cts_AM14MRLlpr[, c(1:4)])
        
  # Remove extra columns from cts ----------------------------------------------
    cts_PL23vNP <- full_join(cts_AM14trans[, c(1, 5:9)], cts_B18trans[, c(1, 10:13)])
    rownames(cts_PL23vNP) <- cts_PL23vNP$ensembl_gene_id
    cts_PL23vNP <- cts_PL23vNP[, 2:10]
    head(cts_PL23vNP)
    
    cts_AM14trans <- cts_AM14trans[, c(5:24)]
    cts_B18trans <- cts_B18trans[, c(5:13)]
    cts_AM14MRLlpr <- cts_AM14MRLlpr[, c(5:21)]
  
    head(cts_AM14trans)
    head(cts_B18trans)
    head(cts_AM14MRLlpr)

# Set up experimental factors --------------------------------------------------
  coldata <- read.csv("raw_data/sample_info.csv")
  coldata <- coldata[, c(1:3, 5:6)]
  coldata
  
  cd_PL23vNP <- coldata[c(1:5, 26:29), ]
  cd_PL23vNP$Cohort <- factor(cd_PL23vNP$Cohort)
  cd_PL23vNP$Treatment <- factor(cd_PL23vNP$Treatment)
  cd_PL23vNP
  
  cd_AM14trans <- coldata[1:20, ]
  cd_AM14trans$Cohort <- factor(cd_AM14trans$Cohort)
  cd_AM14trans$Treatment <- factor(cd_AM14trans$Treatment)
  cd_AM14trans
  
  cd_B18trans <- coldata[21:29, ]
  cd_B18trans$Cohort <- factor(cd_B18trans$Cohort)
  cd_B18trans$Treatment <- factor(cd_B18trans$Treatment)
  cd_B18trans
  
  cd_AM14MRLlpr <- coldata[30:46, ]
  cd_AM14MRLlpr$Cohort <- factor(cd_AM14MRLlpr$Cohort)
  cd_AM14MRLlpr$Treatment <- factor(cd_AM14MRLlpr$Treatment)
  cd_AM14MRLlpr

  # Make sure the columns of cts and rows of coldata are in the same order -----
    for (x in 1:ncol(cts_PL23vNP)) {
      if(colnames(cts_PL23vNP)[x] == cd_PL23vNP$Sample[x]) print(c(x, "true")) 
      else print(c(x, "false"))
    }
    for (x in 1:ncol(cts_AM14trans)) {
      if(colnames(cts_AM14trans)[x] == cd_AM14trans$Sample[x]) print(c(x, "true")) 
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
  dds_PL23vNP <- DESeqDataSetFromMatrix(countData = cts_PL23vNP, 
                                          colData = cd_PL23vNP,
                                          design = ~Treatment)
  dds_PL23vNP
  dds_PL23vNP <- DESeq(dds_PL23vNP)
  resultsNames(dds_PL23vNP)
  
  # Use SVA for PL2-3 vs NP ----------------------------------------------------
    sva_dat  <- counts(dds_PL23vNP, normalized = TRUE)
    idx  <- rowMeans(sva_dat) > 1
    sva_dat  <- sva_dat[idx, ]
    mod  <- model.matrix(~ Treatment, colData(dds_PL23vNP))
    mod0 <- model.matrix(~ 1, colData(dds_PL23vNP))
    svseq <- svaseq(sva_dat, mod, mod0, n.sv = NULL)
    svseq$sv
    
    # Test the number of needed SVs
    # par(mfrow = c(2, 1), mar = c(3,5,3,1))
    # for (i in 1:2) {
    #   stripchart(svseq$sv[, i] ~ dds_PL23vNP$Cohort, vertical = TRUE, main = paste0("SV", i))
    #   abline(h = 0)
    # }
    # dev.off()
    # 
    # # Test results using 1 SV --------------------------------------------------
    #   ddssva_1sv <- dds_PL23vNP
    #   ddssva_1sv$SV1 <- svseq$sv[,1]
    #   design(ddssva_1sv) <- ~ SV1 + Treatment
    #   ddssva_1sv <- DESeq(ddssva_1sv)
    #   resultsNames(ddssva_1sv)
    #   
    #   sv1_res <- results_wrapper(ddssva_1sv, c("Treatment", "PL23", "NP"), full_join(gene_IDs$B18trans, gene_IDs$AM14trans))
    #   summary_wrapper(sv1_res, "PL2-3 vs NP -- 1 SV", "AM14 and B1-8 Adoptive Transfers")
    # 
    # 
    # # Test results using 2 SVs -------------------------------------------------
    #   ddssva_2svs <- dds_PL23vNP
    #   ddssva_2svs$SV1 <- svseq$sv[,1]
    #   ddssva_2svs$SV2 <- svseq$sv[,2]
    #   design(ddssva_2svs) <- ~ SV1 + SV2 + Treatment
    #   ddssva_2svs <- DESeq(ddssva_2svs)
    #   resultsNames(ddssva_2svs)
    #   
    #   sv2_res <- results_wrapper(ddssva_2svs, c("Treatment", "PL23", "NP"), full_join(gene_IDs$B18trans, gene_IDs$AM14trans))
    #   summary_wrapper(sv2_res, "PL2-3 vs NP -- 2 SVs", "AM14 and B1-8 Adoptive Transfers")
    #   
    #   
    #   sv0_res <- results_wrapper(dds_PL23vNP, c("Treatment", "PL23", "NP"), full_join(gene_IDs$B18trans, gene_IDs$AM14trans))
    #   summary_wrapper(sv0_res, "PL2-3 vs NP -- no SVA", "AM14 and B1-8 Adoptive Transfers")
    
    # replace the original PL23 vs NP dds object with the SVA adjusted dds -----
      dds_PL23vNP$SV1 <- svseq$sv[,1]
      dds_PL23vNP$SV2 <- svseq$sv[,2]
      design(dds_PL23vNP) <- ~ SV1 + SV2 + Treatment
      dds_PL23vNP <- DESeq(dds_PL23vNP)
      resultsNames(dds_PL23vNP)
    
  dds_AM14trans <- DESeqDataSetFromMatrix(countData = cts_AM14trans, 
                                          colData = cd_AM14trans,
                                          design = ~Cohort + Treatment)
  dds_AM14trans
  dds_AM14trans <- DESeq(dds_AM14trans)
  resultsNames(dds_AM14trans)

  dds_B18trans <- DESeqDataSetFromMatrix(countData = cts_B18trans, 
                                         colData = cd_B18trans,
                                         design = ~Treatment)
  dds_B18trans
  dds_B18trans <- DESeq(dds_B18trans)
  resultsNames(dds_B18trans)

  dds_AM14MRLlpr <- DESeqDataSetFromMatrix(countData = cts_AM14MRLlpr,
                                           colData = cd_AM14MRLlpr,
                                           design = ~Cohort + Treatment)
  dds_AM14MRLlpr
  dds_AM14MRLlpr <- DESeq(dds_AM14MRLlpr)
  resultsNames(dds_AM14MRLlpr)

# Run QC steps -----------------------------------------------------------------
  QC_heatmaps(dds_PL23vNP, "PL23_vs_NP", "PL2-3 vs NP")
  QC_heatmaps(dds_AM14trans, "AM14_Adoptive_Transfer", "AM14 Adoptive Transfer")
  QC_heatmaps(dds_B18trans, "B18_Adoptive_Transfer", "B1-8 Adoptive Transfer")
  QC_heatmaps(dds_AM14MRLlpr, "AM14_MRLlpr", "AM14 MRL/lpr +/- 2DG")
  
  # QC_heatmaps_batch_corrected(dds_AM14trans, NULL, "AM14 Adoptive Transfer")
  # QC_heatmaps_batch_corrected(dds_AM14MRLlpr, NULL, "AM14 MRL/lpr +/- 2DG")
  
  # QC_PCAplot(dds_PL23vNP, NULL, "PL2-3 vs NP", batch_effect = FALSE)               
  # QC_PCAplot(dds_PL23vNP, NULL, "PL2-3 vs NP", batch_effect = TRUE)               
  QC_PCAplot(dds_PL23vNP, "PL23_vs_NP", "PL2-3 vs NP", batch_effect = FALSE)               
  dev.off()
  QC_PCAplot(dds_AM14trans, "AM14_Adoptive_Transfer-Before_Batch_Correction", 
             "AM14 Adoptive Transfer\n(Before Batch Correction)", 
             batch_effect = FALSE)               
  dev.off()
  QC_PCAplot(dds_AM14trans, "AM14_Adoptive_Transfer-After_Batch_Correction", 
             "AM14 Adoptive Transfer\n(After Batch Correction)", 
             batch_effect = TRUE)               
  dev.off()
  
  QC_PCAplot(dds_B18trans, "B1-8_Adoptive_Transfer","B1-8 Adoptive Transfer", 
             batch_effect = FALSE)               
  dev.off()
  
  QC_PCAplot(dds_AM14MRLlpr, "AM14_MRLlpr-Before_Batch_Correction", 
             "AM14 MRL/lpr mice +/- 2DG\n(Before Batch Correction)", 
             batch_effect = FALSE)               
  dev.off()
  QC_PCAplot(dds_AM14MRLlpr, "AM14_MRLlpr-After_Batch_Correction", 
             "AM14 MRL/lpr mice +/- 2DG\n(After Batch Correction)", 
             batch_effect = TRUE)               
  dev.off()
  
# Differential Expression Analysis ---------------------------------------------
  resultsNames(dds_PL23vNP)
  resultsNames(dds_AM14trans)
  resultsNames(dds_B18trans)
  resultsNames(dds_AM14MRLlpr)
  
# Get DESeq2 results -----------------------------------------------------------
  deseq_res <- list(
    AM14transfer = list(
      PL23_2DG_v_Ctrl = results_wrapper(dds_AM14trans, c("Treatment", "PL23_2DG", "PL23"), gene_IDs$AM14trans),
      R848_2DG_v_Ctrl = results_wrapper(dds_AM14trans, c("Treatment", "R848_2DG", "R848"), gene_IDs$AM14trans),
      PL23_vs_R848 = results_wrapper(dds_AM14trans, c("Treatment", "PL23", "R848"), gene_IDs$AM14trans)),
    B18transfer = results_wrapper(dds_B18trans, c("Treatment", "NP_2DG", "NP"), gene_IDs$B18trans),
    PL23_v_NP = results_wrapper(dds_PL23vNP, c("Treatment", "PL23", "NP"), full_join(gene_IDs$B18trans, gene_IDs$AM14trans)),
    AM14MRLlpr = results_wrapper(dds_AM14MRLlpr, c("Treatment", "2DG", "Control"), gene_IDs$AM14MRLlpr)
  )
  
  head(deseq_res$AM14transfer$PL23_2DG_v_Ctrl)
  head(deseq_res$B18transfer)
  head(deseq_res$AM14MRLlpr)

# Make a data.frame that summarizes the numbers of DEGs detected -------------
  full_summary <- join_all(list(summary_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "PL2-3+2DG vs PL2-3", "AM14 Adoptive Transfer"),
                                summary_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl, "R848+2DG vs R848", "AM14 Adoptive Transfer"),
                                summary_wrapper(deseq_res$AM14transfer$PL23_vs_R848, "PL2-3 vs R848", "AM14 Adoptive Transfer"),
                                summary_wrapper(deseq_res$B18transfer, "NP+2DG vs NP", "B1-8 Adoptive Transfer"),
                                summary_wrapper(deseq_res$PL23_v_NP, "PL2-3 vs NP", "AM14 and B1-8 Adoptive Transfers"),
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
  addWorksheet(wb, "AM14 PL2-3 vs B1-8 NP")
  addWorksheet(wb, "AM14 MRLlpr - 2DG_vs_Ctrl")
  
  writeData(wb, "DEG Summaries", full_summary)
  writeData(wb, "AM14 - PL23_2DG_vs_Ctrl", deseq_res$AM14transfer$PL23_2DG_v_Ctrl)
  writeData(wb, "AM14 - R848_2DG_vs_Ctrl", deseq_res$AM14transfer$R848_2DG_v_Ctrl)
  writeData(wb, "AM14 - PL23_vs_R848", deseq_res$AM14transfer$PL23_vs_R848)
  writeData(wb, "B1-8 - 2DG_vs_Ctrl", deseq_res$B18transfer)
  writeData(wb, "AM14 PL2-3 vs B1-8 NP", deseq_res$PL23_v_NP)
  writeData(wb, "AM14 MRLlpr - 2DG_vs_Ctrl", deseq_res$AM14MRLlpr)
  
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
  # writeData(wb, "Raw_Counts - AM14 transfer", counts(dds_AM14trans, normalized = FALSE), rowNames = TRUE)
  # writeData(wb, "Norm_Counts - AM14 transfer", counts(dds_AM14trans, normalized = TRUE), rowNames = TRUE)
  # writeData(wb, "Raw_Counts - B1-8 transfer", counts(dds_B18trans, normalized = FALSE), rowNames = TRUE)
  # writeData(wb, "Norm_Counts - B1-8 transfer", counts(dds_B18trans, normalized = TRUE), rowNames = TRUE)
  # writeData(wb, "Raw_Counts - AM14 MRLlpr", counts(dds_AM14MRLlpr, normalized = FALSE), rowNames = TRUE)
  # writeData(wb, "Norm_Counts - AM14 MRLlpr", counts(dds_AM14MRLlpr, normalized = TRUE), rowNames = TRUE)
  # 
  # saveWorkbook(wb, "Output/DESeq2_gene_counts.xlsx", overwrite = TRUE)
  # rm(wb)

# # Make DEG lists and export to CSV files -------------------------------------
#   write_DEG_CSV(res_PL23, 1, "PL2-3_2DG_vs_Control")
#   write_DEG_CSV(res_R848, 1, "R848_2DG_vs_Control")
#   write_DEG_CSV(res_PL23vR848, 1, "PL2-3_Ctrl_vs_R848_Ctrl")
#   write_DEG_CSV(res_2DG_PL23vR848, 1, "PL2-3_2DG_vs_R848_2DG")