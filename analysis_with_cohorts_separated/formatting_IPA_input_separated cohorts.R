
create_ipa_input <- function(ID_cols, file_name){
  if(missing(file_name)) file_name <- NULL
  shared_cols <- colnames(deseq_res$AM14transfer$co1$PL23_2DG_v_Ctrl[, ID_cols])

  PL23_co1 <- deseq_res$AM14transfer$co1$PL23_2DG_v_Ctrl[, c(ID_cols, 8, 11:12)]
  R848_co1 <- deseq_res$AM14transfer$co1$R848_2DG_v_Ctrl[, c(ID_cols, 8, 11:12)]
  MRLlpr_co1 <- deseq_res$AM14MRLlpr$co1[, c(ID_cols, 8, 11:12)]
  
  PL23_co2 <- deseq_res$AM14transfer$co2$PL23_2DG_v_Ctrl[, c(ID_cols, 8, 11:12)]
  R848_co2 <- deseq_res$AM14transfer$co2$R848_2DG_v_Ctrl[, c(ID_cols, 8, 11:12)]
  MRLlpr_co2 <- deseq_res$AM14MRLlpr$co2[, c(ID_cols, 8, 11:12)]
  
  NP <- deseq_res$B18transfer[, c(ID_cols, 8, 11:12)]

  colnames(PL23_co1) <- c(shared_cols, "PL23_2DGvCtrl_co1.log2FoldChange", "PL23_2DGvCtrl_co1.pvalue", "PL23_2DGvCtrl_co1.fdr")
  colnames(R848_co1) <- c(shared_cols, "R848_2DGvCtrl_co1.log2FoldChange", "R848_2DGvCtrl_co1.pvalue", "R848_2DGvCtrl_co1.fdr")
  colnames(MRLlpr_co1) <- c(shared_cols, "MRLlpr_2DGvCtrl_co1.log2FoldChange", "MRLlpr_2DGvCtrl_co1.pvalue", "MRLlpr_2DGvCtrl_co1.fdr")
  
  colnames(PL23_co2) <- c(shared_cols, "PL23_2DGvCtrl_co2.log2FoldChange", "PL23_2DGvCtrl_co2.pvalue", "PL23_2DGvCtrl_co2.fdr")
  colnames(R848_co2) <- c(shared_cols, "R848_2DGvCtrl_co2.log2FoldChange", "R848_2DGvCtrl_co2.pvalue", "R848_2DGvCtrl_co2.fdr")
  colnames(MRLlpr_co2) <- c(shared_cols, "MRLlpr_2DGvCtrl_co2.log2FoldChange", "MRLlpr_2DGvCtrl_co2.pvalue", "MRLlpr_2DGvCtrl_co2.fdr")
  
  colnames(NP) <- c(shared_cols, "NP_2DGvCtrl.log2FoldChange", "NP_2DGvCtrl.pvalue", "NP_2DGvCtrl.fdr")

  res <- join_all(list(PL23_co1, PL23_co2, R848_co1, R848_co2, MRLlpr_co1, MRLlpr_co2, NP), type = "full")
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

# For analysis with cohorts separated ------------------------------------------
  setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_separated")
  getwd()
  
  mouse_all <- create_ipa_input(c(1:3), "IPA_input_mouse_IDs_analysis_with_cohorts_separated")
  # human_all <- create_ipa_input(c(4:6), "human_all_IDs")
  
head(mouse_all)
