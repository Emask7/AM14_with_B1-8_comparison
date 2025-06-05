# This code runs DESeq2 on AM14 transfer PL2-3 vs B1-8 transfer NP
## This code isn't necessary, and I don't recommend running it.
## I'm only keeping the file so I don't have to re-write parts later if I end up wanting to try this analysis again.
## Since the experiments were performed completely separately and samples were sequenced separately, it probably isn't statistically sound to perform this analysis.

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
  
# Set up experimental factors --------------------------------------------------
  coldata <- read.csv("raw_data/sample_info.csv")
  coldata <- coldata[, c(1:3, 5:6, 10)]
  coldata
  
  cd_PL23vNP <- coldata[c(1:5, 26:29), ]
  cd_PL23vNP$Cohort <- factor(cd_PL23vNP$Cohort)
  cd_PL23vNP$Treatment <- factor(cd_PL23vNP$Treatment)
  cd_PL23vNP
    
  # Make sure the columns of cts and rows of coldata are in the same order -----
    for (x in 1:ncol(cts_PL23vNP)) {
      if(colnames(cts_PL23vNP)[x] == cd_PL23vNP$Sample[x]) print(c(x, "true")) 
      else print(c(x, "false"))
    }
  

# Set up DESeq Data Set --------------------------------------------------------
  dds_PL23vNP <- DESeqDataSetFromMatrix(countData = cts_PL23vNP, 
                                        colData = cd_PL23vNP,
                                        design = ~Treatment)
  dds_PL23vNP
  dds_PL23vNP <- DESeq(dds_PL23vNP)
  resultsNames(dds_PL23vNP)
  
  dds_PL23vNP_tempRes <- results_wrapper(dds_PL23vNP, c("Treatment", "PL23", "NP"), full_join(gene_IDs$B18trans, gene_IDs$AM14trans))
  summary_wrapper(dds_PL23vNP_tempRes, "PL2-3 vs NP", "AM14 and B1-8 Adoptive Transfers")

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

# Run QC steps -----------------------------------------------------------------
  QC_heatmaps(dds_PL23vNP, "PL23_vs_NP", "PL2-3 vs NP")
  
  # QC_PCAplot(dds_PL23vNP, "PL23_vs_NP", "PL2-3 vs NP", batch_effect = NULL)
  # dev.off()
  
# Differential Expression Analysis ---------------------------------------------
  resultsNames(dds_PL23vNP)

# Get DESeq2 results -----------------------------------------------------------
  deseq_res_PL23_v_NP <- results_wrapper(dds_PL23vNP, c("Treatment", "PL23", "NP"), full_join(gene_IDs$B18trans, gene_IDs$AM14trans))
  

  