# Function to remove unwanted gene IDs -----------------------------------------
  import_Rosalind_data <- function(res_file){
    temp_cts <- read.delim(res_file)
    temp_cts <- temp_cts %>%
      filter(!grepl("Igh", external_gene_name)) %>%
      filter(!grepl("Igk", external_gene_name)) %>%
      filter(!grepl("7SK", external_gene_name)) %>%
      filter(!grepl("5_8S_rRNA", external_gene_name)) %>%
      filter(!grepl("5S_rRNA", external_gene_name))
    temp_cts <- temp_cts[, c(1:4, 12:ncol(temp_cts))]
    temp_cts <- temp_cts[!duplicated(temp_cts), ]
    # keep <- rowSums(temp_cts[, 5:ncol(temp_cts)]) > 0
    # temp_cts <- temp_cts[keep,]
    temp_cts <- temp_cts[!duplicated(temp_cts$ensembl_gene_id), ]
    rownames(temp_cts) <- temp_cts$ensembl_gene_id
    return(temp_cts)
  }

# Import raw counts from Rosalind output ---------------------------------------
  cts_AM14trans <- import_Rosalind_data("raw_data/AM14_Adoptive_Transfer_rawCountsWithAnnotations.txt")
  cts_B18trans <- import_Rosalind_data("raw_data/B1-8_Adoptive_Transfer_rawCountsWithAnnotations.txt")
  cts_AM14MRLlpr <- import_Rosalind_data("raw_data/AM14_MRLlpr_rawCountsWithAnnotations.txt")

  # cts <- full_join(cts_AM14trans, cts_B18trans)
  # cts <- full_join(cts, cts_AM14MRLlpr)
  # colnames(cts)
  # head(cts)
  # nrow(cts)
  
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
    
# # Set row names to Ensembl IDs -----------------------------------------------
#   cts_dups <- cts[duplicated(cts$ensembl_gene_id), ]
#   nrow(cts_dups)
#   cts_dups
#   rowSums(cts_dups > 0)
#   
#   cts <- cts[!duplicated(cts$ensembl_gene_id), ]
#   nrow(cts)
#   cts <- cts[!is.na(cts$ensembl_gene_id), ]
#   nrow(cts)
#   rownames(cts) <- cts$ensembl_gene_id
# 
#   # Save a table of gene IDs in case they need to be converted later ---------
#     gene_IDs <- cts[, c(1:4)]
#     head(gene_IDs)
#     
#   # Remove extra columns from cts --------------------------------------------
#     cts <- cts[, c(5:50)]
# 
# # For gene symbols as row names ----------------------------------------------
#   cts_AM14trans <- cts_AM14trans[!duplicated(cts_AM14trans$external_gene_name), ]
#   nrow(cts_AM14trans)
#   cts_AM14trans <- cts_AM14trans[!is.na(cts_AM14trans$external_gene_name), ]
#   nrow(cts_AM14trans)
#   rownames(cts_AM14trans) <- cts_AM14trans$external_gene_name
#   nrow(cts_AM14trans)
# 
# # For Entrez IDs as row names ------------------------------------------------
#   cts_AM14trans <- cts_AM14trans[!duplicated(cts_AM14trans$entrezgene), ]
#   nrow(cts_AM14trans)
#   cts_AM14trans <- cts_AM14trans[!is.na(cts_AM14trans$entrezgene), ]
#   nrow(cts_AM14trans)
#   rownames(cts_AM14trans) <- cts_AM14trans$entrezgene
#   nrow(cts_AM14trans)

# Set up experimental factors --------------------------------------------------
  coldata <- read.csv("raw_data/sample_info.csv")
  coldata <- coldata[, c(1, 3, 5:6)]
  coldata
  
  cd_AM14trans <- coldata[1:20, ]
  cd_AM14trans$Cohort <- factor(cd_AM14trans$Cohort)
  cd_AM14trans$Treatment <- factor(cd_AM14trans$Treatment)
  cd_AM14trans
  
  cd_B18trans <- coldata[21:29,  c(1, 2, 4)]
  cd_B18trans$Treatment <- factor(cd_B18trans$Treatment)
  cd_B18trans
  
  cd_AM14MRLlpr <- coldata[30:46, ]
  cd_AM14MRLlpr$Cohort <- factor(cd_AM14MRLlpr$Cohort)
  cd_AM14MRLlpr$Treatment <- factor(cd_AM14MRLlpr$Treatment)
  cd_AM14MRLlpr

  # Make sure the columns of cts and rows of coldata are in the same order -----
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
  dds_AM14trans <- DESeqDataSetFromMatrix(countData = cts_AM14trans, 
                                          colData = cd_AM14trans,
                                          design = ~Cohort + Treatment)
  dds_AM14trans
  keep <- rowSums(counts(dds_AM14trans) >= 10) >= 5
  dds_AM14trans <- dds_AM14trans[keep,]
  dds_AM14trans <- DESeq(dds_AM14trans)
  resultsNames(dds_AM14trans)
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
                                           design = ~Cohort + Treatment)
  dds_AM14MRLlpr
  keep <- rowSums(counts(dds_AM14MRLlpr) >= 10) >= 8
  dds_AM14MRLlpr <- dds_AM14MRLlpr[keep,]
  dds_AM14MRLlpr <- DESeq(dds_AM14MRLlpr)
  resultsNames(dds_AM14MRLlpr)
  rm(keep)
  
# Run QC steps -----------------------------------------------------------------
  
  # col_factors = c("Treatment", "Cohort")
  
  # count_matrix_heatmap <- function(dds, folder_name, file_name_start, col_factors)
  
  count_matrix_heatmap(dds_AM14trans, "AM14_Adoptive_Transfer", c("Treatment", "Cohort"))
                                   
  
# Differential Expression Analysis ---------------------------------------------
  resultsNames(dds)
    
# Get results (default methods) ------------------------------------------------
  res_PL23 <- results(dds, contrast = c("Treatment", "PL2-3+2DG", "PL2-3"))
  res_PL23 <- data.frame(subset(res_PL23, !is.na(padj)))
  summary_v2(res_PL23, "PL2-3+2DG vs PL2-3")
  
  res_R848 <- results(dds, contrast = c("Treatment", "R848+2DG", "R848"))
  res_R848 <- data.frame(subset(res_R848, !is.na(padj)))
  summary_v2(res_R848, "R848+2DG vs R848")
  
  res_PL23vR848 <- results(dds, contrast = c("Treatment", "PL2-3", "R848"))
  res_PL23vR848 <- data.frame(subset(res_PL23vR848, !is.na(padj)))
  summary_v2(res_PL23vR848, "PL2-3 vs R848")
  
  res_2DG_PL23vR848 <- results(dds, contrast = c("Treatment", "PL2-3+2DG", "R848+2DG"))
  res_2DG_PL23vR848 <- data.frame(subset(res_2DG_PL23vR848, !is.na(padj)))
  summary_v2(res_2DG_PL23vR848, "PL2-3+2DG vs R848+2DG")
  
    # Note: There are a few reasons why a p value or padj value would be NA
    # According to the DESeq2 manual, these are the reasons:
    # If within a row, all samples have zero counts, the baseMean column will
    # be zero, and the LFC estimates, p value and padj will all be NA.
    # If a row contains a sample with an extreme count outlier then the 
    # p value and padj will be set to NA.
    # If a row is filtered by automatic independent filtering, for having a 
    # low mean normalized count, then only padj will be set to NA.
    
    
  # Make DEG lists and export to CSV files -------------------------------------
    write_DEG_CSV(res_PL23, 1, "PL2-3_2DG_vs_Control")
    write_DEG_CSV(res_R848, 1, "R848_2DG_vs_Control")
    write_DEG_CSV(res_PL23vR848, 1, "PL2-3_Ctrl_vs_R848_Ctrl")
    write_DEG_CSV(res_2DG_PL23vR848, 1, "PL2-3_2DG_vs_R848_2DG")

# # Get results (Threshold-Based Wald tests) -----------------------------------
#   res_PL23_lfc <- results(dds, contrast = c("Treatment", "PL2-3+2DG", "PL2-3"),
#                           lfcThreshold = 0.6, altHypothesis = "greaterAbs", alpha = 0.05)
#   res_R848_lfc <- results(dds, contrast = c("Treatment", "R848+2DG", "R848"),
#                           lfcThreshold = 0.6, altHypothesis = "greaterAbs", alpha = 0.05)
#   res_PL23vR848_lfc <- results(dds, contrast = c("Treatment", "PL2-3", "R848"),
#                                lfcThreshold = 0.6, altHypothesis = "greaterAbs", alpha = 0.05)
#   res_2DG_PL23vR848_lfc <- results(dds, contrast = c("Treatment", "PL2-3+2DG", "R848+2DG"),
#                                    lfcThreshold = 1, altHypothesis = "greaterAbs", alpha = 0.05)
# 
#   summary(res_PL23_lfc)
#   summary(res_R848_lfc)
#   summary(res_PL23vR848_lfc)
#   summary(res_2DG_PL23vR848_lfc)

# Save count data to Excel file ------------------------------------------------
  # Add columns with alternative gene identifiers to the DESeq2 results --------
    res_PL23_full <- tibble::rownames_to_column(res_PL23, var = "ensembl_gene_id")
    res_PL23_full <- right_join(gene_IDs, res_PL23_full)
    head(res_PL23_full)
    
    res_R848_full <- tibble::rownames_to_column(res_R848, var = "ensembl_gene_id")
    res_R848_full <- right_join(gene_IDs, res_R848_full)
    head(res_R848_full)
    
    res_PL23vR848_full <- tibble::rownames_to_column(res_PL23vR848, var = "ensembl_gene_id")
    res_PL23vR848_full <- right_join(gene_IDs, res_PL23vR848_full)
    head(res_PL23vR848_full)
    
    res_2DG_PL23vR848_full <- tibble::rownames_to_column(res_2DG_PL23vR848, var = "ensembl_gene_id")
    res_2DG_PL23vR848_full <- right_join(gene_IDs, res_2DG_PL23vR848_full)
    head(res_2DG_PL23vR848_full)
    
  # Make a data.frame that summarizes the numbers of DEGs detected -------------
    full_summary <- full_join(summary_v2(res_PL23, "PL2-3+2DG vs PL2-3"), 
                                summary_v2(res_PL23, "PL2-3+2DG vs PL2-3", lfc_cutoff = 1))
    
    full_summary_2 <- full_join(summary_v2(res_R848, "R848+2DG vs R848"), 
                                summary_v2(res_R848, "R848+2DG vs R848", lfc_cutoff = 1))
    
    full_summary_3 <- full_join(summary_v2(res_PL23vR848, "PL2-3 vs R848"), 
                                summary_v2(res_PL23vR848, "PL2-3 vs R848", lfc_cutoff = 1))
    
    full_summary_4 <- full_join(summary_v2(res_2DG_PL23vR848, "PL2-3+2DG vs R848+2DG"), 
                                summary_v2(res_2DG_PL23vR848, "PL2-3+2DG vs R848+2DG", lfc_cutoff = 1))
    
    full_summary <- full_join(full_summary, full_summary_2)
    full_summary <- full_join(full_summary, full_summary_3)
    full_summary <- full_join(full_summary, full_summary_4)
    full_summary
    
    rm(full_summary_2, full_summary_3, full_summary_4)
    
    
  # Save results to an Excel file ----------------------------------------------
    wb <- createWorkbook("Output/DESeq2_results.xlsx")
    
    addWorksheet(wb, "DEG Summaries")
    addWorksheet(wb, "PL23_2DG_vs_Ctrl")
    addWorksheet(wb, "R848_2DG_vs_Ctrl")
    addWorksheet(wb, "PL23_Ctrl_vs_R848_Ctrl")
    addWorksheet(wb, "PL23_2DG_vs_R848_2DG")
    addWorksheet(wb, "Raw_Gene_Counts")
    addWorksheet(wb, "Normalized_Gene_Counts")
    
    writeData(wb, "DEG Summaries", full_summary)
    writeData(wb, "PL23_2DG_vs_Ctrl", res_PL23_full)
    writeData(wb, "R848_2DG_vs_Ctrl", res_R848_full)
    writeData(wb, "PL23_Ctrl_vs_R848_Ctrl", res_PL23vR848_full)
    writeData(wb, "PL23_2DG_vs_R848_2DG", res_2DG_PL23vR848_full)
    writeData(wb, "Raw_Gene_Counts", 
              counts(dds, normalized = FALSE), rowNames = TRUE)
    writeData(wb, "Normalized_Gene_Counts", 
              counts(dds, normalized = TRUE), rowNames = TRUE)
    
    saveWorkbook(wb, "Output/DESeq2_results.xlsx", overwrite = TRUE)
    
# Make DEG Venn Diagram --------------------------------------------------------
  DEG_lists <- list(
    # PL23_2DGvsCtrl = write_DEG_CSV(res_PL23_full, 1, "PL2-3_2DG_vs_Control"),
    # R848_2DGvsCtrl = write_DEG_CSV(res_R848_full, 1, "R848_2DG_vs_Control"),
    # PL23vR848 = write_DEG_CSV(res_PL23vR848_full, 1, "PL2-3_Ctrl_vs_R848_Ctrl")
    PL23_2DGvsCtrl = write_sig_LFCs(res_PL23_full, 1, "external_gene_name", "PL2-3_2DG_vs_Control"),
    R848_2DGvsCtrl = write_sig_LFCs(res_R848_full, 1, "external_gene_name", "R848_2DG_vs_Control"),
    PL23vR848 = write_sig_LFCs(res_PL23vR848_full, 1, "external_gene_name", "PL2-3_Ctrl_vs_R848_Ctrl")
  )
  nrow(DEG_lists$PL23_2DGvsCtrl$down)
  nrow(DEG_lists$R848_2DGvsCtrl$down)
  nrow(DEG_lists$PL23vR848$down)
  
  head(DEG_lists$PL23_2DGvsCtrl)
  
  
  write_DEG_CSV(res_PL23, 1, "PL2-3_2DG_vs_Control")
  write_DEG_CSV(res_R848, 1, "R848_2DG_vs_Control")
  write_DEG_CSV(res_PL23vR848, 1, "PL2-3_Ctrl_vs_R848_Ctrl")
  write_DEG_CSV(res_2DG_PL23vR848, 1, "PL2-3_2DG_vs_R848_2DG")
  
    
  venn_up <- venndetail(list("PL2-3+2DG vs PL2-3" = DEG_lists$PL23_2DGvsCtrl$up$ensembl_gene_id,
                             "R848+2DG vs R848" = DEG_lists$R848_2DGvsCtrl$up$ensembl_gene_id,
                             "PL2-3 vs R848" = DEG_lists$PL23vR848$up$ensembl_gene_id))
  plot(venn_up, mycol = c("goldenrod1", "darkorange1", "red"), 
       filename = "Output/venn_diagram_upreg.png",
       margin = 0.1, cat.cex = 0.5)
  dev.off()
  detail(venn_up)
  

    
  venn_down <- venndetail(list("PL2-3+2DG vs PL2-3" = DEG_lists$PL23_2DGvsCtrl$down$ensembl_gene_id,
                               "R848+2DG vs R848" = DEG_lists$R848_2DGvsCtrl$down$ensembl_gene_id,
                               "PL2-3 vs R848" = DEG_lists$PL23vR848$down$ensembl_gene_id))
  plot(venn_down, mycol = c("darkseagreen1", "dodgerblue", "orchid"), 
       filename = "Output/venn_diagram_downreg.png",
       margin = 0.1, cat.cex = 0.5)
  dev.off()
  detail(venn_down)
