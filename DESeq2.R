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
    temp_cts <- temp_cts[!duplicated(temp_cts$ensembl_gene_id), ]
    rownames(temp_cts) <- temp_cts$ensembl_gene_id
    return(temp_cts)
  }

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
  coldata <- coldata[, c(1:3, 5:6)]
  coldata
  
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
  QC_heatmaps(dds_AM14trans, "AM14_Adoptive_Transfer", "AM14 Adoptive Transfer", c("Treatment", "Cohort"))
  QC_heatmaps(dds_B18trans, "B18_Adoptive_Transfer", "B1-8 Adoptive Transfer", c("Treatment", "Cohort"))
  QC_heatmaps(dds_AM14MRLlpr, "AM14_MRLlpr", "AM14 MRL/lpr +/- 2DG", c("Treatment", "Cohort"))
  
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
  resultsNames(dds_AM14trans)
  resultsNames(dds_B18trans)
  resultsNames(dds_AM14MRLlpr)
  
# Get DESeq2 results -----------------------------------------------------------
  get_res <- function(dds, cons, IDs){
    res <- results(dds, contrast = cons)
    res <- data.frame(subset(res, !is.na(padj)))
      # Note: There are a few reasons why a p value or padj value would be NA
      # According to the DESeq2 manual, these are the reasons:
      # If within a row, all samples have zero counts, the baseMean column will be zero, and the LFC estimates, p value and padj will all be NA.
      # If a row contains a sample with an extreme count outlier then the p value and padj will be set to NA.
      # If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only padj will be set to NA.
    tempcols <- colnames(res)
    res <- data.frame(rownames(res), res)
    colnames(res) <- c("ensembl_gene_id", tempcols)
    return(right_join(IDs, res))
  }
  
  deseq_res <- list(
    AM14transfer = list(
      PL23_2DG_v_Ctrl = get_res(dds_AM14trans, c("Treatment", "PL23_2DG", "PL23"), gene_IDs$AM14trans),
      R848_2DG_v_Ctrl = get_res(dds_AM14trans, c("Treatment", "R848_2DG", "R848"), gene_IDs$AM14trans),
      PL23_vs_R848 = get_res(dds_AM14trans, c("Treatment", "PL23", "R848"), gene_IDs$AM14trans)),
    B18transfer = get_res(dds_B18trans, c("Treatment", "NP_2DG", "NP"), gene_IDs$B18trans),
    AM14MRLlpr = get_res(dds_AM14MRLlpr, c("Treatment", "2DG", "Control"), gene_IDs$AM14MRLlpr)
  )
  
  head(deseq_res$AM14transfer$PL23_2DG_v_Ctrl)
  head(deseq_res$B18transfer)
  head(deseq_res$AM14MRLlpr)

# Make a data.frame that summarizes the numbers of DEGs detected -------------
  full_summary <- join_all(list(summary_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "PL2-3+2DG vs PL2-3", "AM14 Adoptive Transfer"),
                                summary_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl, "R848+2DG vs R848", "AM14 Adoptive Transfer"),
                                summary_wrapper(deseq_res$AM14transfer$PL23_vs_R848, "PL2-3 vs R848", "AM14 Adoptive Transfer"),
                                summary_wrapper(deseq_res$B18transfer, "NP+2DG vs NP", "B1-8 Adoptive Transfer"),
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
  
  writeData(wb, "Raw_Counts - AM14 transfer", counts(dds_AM14trans, normalized = FALSE), rowNames = TRUE)
  writeData(wb, "Norm_Counts - AM14 transfer", counts(dds_AM14trans, normalized = TRUE), rowNames = TRUE)
  writeData(wb, "Raw_Counts - B1-8 transfer", counts(dds_B18trans, normalized = FALSE), rowNames = TRUE)
  writeData(wb, "Norm_Counts - B1-8 transfer", counts(dds_B18trans, normalized = TRUE), rowNames = TRUE)
  writeData(wb, "Raw_Counts - AM14 MRLlpr", counts(dds_AM14MRLlpr, normalized = FALSE), rowNames = TRUE)
  writeData(wb, "Norm_Counts - AM14 MRLlpr", counts(dds_AM14MRLlpr, normalized = TRUE), rowNames = TRUE)
  
  saveWorkbook(wb, "Output/DESeq2_gene_counts.xlsx", overwrite = TRUE)
  rm(wb)

# # Make DEG lists and export to CSV files -------------------------------------
#   write_DEG_CSV(res_PL23, 1, "PL2-3_2DG_vs_Control")
#   write_DEG_CSV(res_R848, 1, "R848_2DG_vs_Control")
#   write_DEG_CSV(res_PL23vR848, 1, "PL2-3_Ctrl_vs_R848_Ctrl")
#   write_DEG_CSV(res_2DG_PL23vR848, 1, "PL2-3_2DG_vs_R848_2DG")
    
# Make DEG Venn Diagram --------------------------------------------------------
  DEG_lists <- list(
    PL23_2DG_vs_Ctrl = deg_list_for_venn(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, 1),
    R848_2DG_vs_Ctrl = deg_list_for_venn(deseq_res$AM14transfer$R848_2DG_v_Ctrl, 1),
    PL23_vs_R848 = deg_list_for_venn(deseq_res$AM14transfer$PL23_vs_R848, 1),
    B18 = deg_list_for_venn(deseq_res$B18transfer, 1),
    MRLlpr = deg_list_for_venn(deseq_res$AM14MRLlpr, 1)
  )
  head(DEG_lists$B18$up)
  
  venn_up <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG\nvs PL2-3\n" = DEG_lists$PL23_2DG_vs_Ctrl$up,
                             "AM14 Transfer:\nR848 + 2DG\nvs R848\n" = DEG_lists$R848_2DG_vs_Ctrl$up,
                             "AM14 MRL/lpr:\n2DG vs Control\n" = DEG_lists$MRLlpr$up,
                             "B1-8 Transfer:\nNP+2DG vs NP\n" = DEG_lists$B18$up))
  plot(venn_up, mycol = met.brewer("Johnson", n = 8, direction = 1),
       filename = "Output/venn_diagram_upreg.png",
       margin = 0.05, cat.cex = 1, cex = 1.5)
  dev.off()
  
  venn_down <- venndetail(list("AM14 Transfer:\nPL2-3 + 2DG\nvs PL2-3\n" = DEG_lists$PL23_2DG_vs_Ctrl$down,
                               "AM14 Transfer:\nR848 + 2DG\nvs R848\n" = DEG_lists$R848_2DG_vs_Ctrl$down,
                               "AM14 MRL/lpr:\n2DG vs Control\n" = DEG_lists$MRLlpr$down,
                               "B1-8 Transfer:\nNP + 2DG vs NP\n" = DEG_lists$B18$down))
  plot(venn_down,  mycol = met.brewer("Johnson", n = 8, direction = -1),
       filename = "Output/venn_diagram_downreg.png",
       margin = 0.05, cat.cex = 1, cex = 1.5)
  dev.off()

  detail(venn_up)
  detail(venn_down)

# Volcano plots ----------------------------------------------------------------
  vp_data = list(
      AM14transfer = list(
        PL23_2DG_v_Ctrl = deseq_res$AM14transfer$PL23_2DG_v_Ctrl[, c(2, 6, 10)],
        R848_2DG_v_Ctrl = deseq_res$AM14transfer$R848_2DG_v_Ctrl[, c(2, 6, 10)],
        PL23_vs_R848 = deseq_res$AM14transfer$PL23_vs_R848[, c(2, 6, 10)]),
      B18transfer = deseq_res$B18transfer[, c(2, 6, 10)],
      AM14MRLlpr = deseq_res$AM14MRLlpr[, c(2, 6, 10)])
  head(vp_data$AM14MRLlpr)

  volcano_wrapper <- function(dat, plot_title, x_limits, y_limits, file_name){
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
                          pCutoff = 0.05, FCcutoff = 1,
                          x = 'log2FoldChange', y = 'padj', 
                          xlim = x_limits, ylim = y_limits, 
                          title = plot_title, subtitle = "(Adjusted p-values)", 
                          legendLabels = c("NS", expression(Log[2] ~ FC > 1),
                                           expression(p - value < ~ 0.05),
                                           expression(p - value < ~ 0.05 ~ and ~ log[2] ~ FC > ~ 1)),
                          drawConnectors = TRUE, min.segment.length = 1, 
                          max.overlaps = 12, labSize = 4)
    
    if(!missing(file_name)){
      if(file_name != FALSE & !is.null(file_name)){
        ggsave(stri_join("Output/Volcano_plots/Volcano plot - ", file_name, ".png"),
               units = "px", width = 3000, height = 3000, dpi = 300)
      } 
    } 
    vp
  }  
  
  volcano_wrapper(vp_data$AM14transfer$PL23_2DG_v_Ctrl, "AM14 Transfer: PL2-3 + 2DG vs PL2-3", 
                  c(-10, 5), c(0, 8), "AM14_Transfer-PL2-3_2DG_vs_Control")
  
  volcano_wrapper(vp_data$AM14transfer$R848_2DG_v_Ctrl, "AM14 Transfer: R848 + 2DG vs R848", 
                  c(-10, 5), c(0, 5), "AM14_Transfer-R848_2DG_vs_Control")
  
  volcano_wrapper(vp_data$AM14transfer$PL23_vs_R848, "AM14 Transfer: PL2-3 vs R848", 
                  c(-11, 6), c(0, 6), "AM14_Transfer-PL2-3_vs_R848")
  
  volcano_wrapper(vp_data$B18transfer, "B1-8 Transfer: NP + 2DG vs NP", 
                  c(-5, 3), c(0, 8), "B1-8_Transfer-NP_2DG_vs_NP")
  
  volcano_wrapper(vp_data$AM14MRLlpr, "AM14 MRL/lpr: 2DG vs Control", 
                  c(-5, 2.5), c(0, 15), "AM14_MRLlpr")
  