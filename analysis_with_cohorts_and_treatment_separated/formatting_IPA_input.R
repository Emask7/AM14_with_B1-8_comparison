setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_and_treatment_separated")
getwd()

create_ipa_input <- function(ID_cols, file_name){
  if(missing(file_name)) file_name <- NULL
  shared_cols <- colnames(deseq_res$AM14transfer_PL23[, ID_cols])

  PL23 <- deseq_res$AM14transfer_PL23[, c(ID_cols, 8, 11:12)]
  R848 <- deseq_res$AM14transfer_R848[, c(ID_cols, 8, 11:12)]
  NP <- deseq_res$B18transfer[, c(ID_cols, 8, 11:12)]
  MRLlpr <- deseq_res$AM14MRLlpr[, c(ID_cols, 8, 11:12)]

  colnames(PL23) <- c(shared_cols, "PL23_2DGvCtrl.log2FoldChange", "PL23_2DGvCtrl.pvalue", "PL23_2DGvCtrl.fdr")
  colnames(R848) <- c(shared_cols, "R848_2DGvCtrl.log2FoldChange", "R848_2DGvCtrl.pvalue", "R848_2DGvCtrl.fdr")
  colnames(NP) <- c(shared_cols, "NP_2DGvCtrl.log2FoldChange", "NP_2DGvCtrl.pvalue", "NP_2DGvCtrl.fdr")
  colnames(MRLlpr) <- c(shared_cols, "MRLlpr_2DGvCtrl.log2FoldChange", "MRLlpr_2DGvCtrl.pvalue", "MRLlpr_2DGvCtrl.fdr")

  res <- join_all(list(PL23, R848, NP, MRLlpr), type = "full")
  res <- filter(res, )

  if(!is.null(file_name)) {
    if(!file.exists("Output/")) dir.create("Output/")
    if(!file.exists("Output/IPA_input/")) dir.create("Output/IPA_input/")

    wb <- createWorkbook(stri_join(c("Output/IPA_input/", file_name, ".xlsx"), collapse = ""))
    addWorksheet(wb, "sheet1")
    writeData(wb, "sheet1", res)
    saveWorkbook(wb, stri_join(c("Output/IPA_input/", file_name, ".xlsx"), collapse = ""), overwrite = TRUE)
    rm(wb)
  }

  return(res)
}

mouse_all <- create_ipa_input(c(1:3), "mouse_all_IDs")
human_all <- create_ipa_input(c(4:6), "human_all_IDs")

# mouse_ensembl <- create_ipa_input(c(1), "mouse_ensembl")
# mouse_symbol <- create_ipa_input(c(2), "mouse_symbol")
# mouse_entrez <- create_ipa_input(c(3), "mouse_entrez")

# human_ensembl <- create_ipa_input(c(4), "human_ensembl")
# human_symbol <- create_ipa_input(c(5), "human_symbol")
# human_entrez <- create_ipa_input(c(6), "human_entrez")
