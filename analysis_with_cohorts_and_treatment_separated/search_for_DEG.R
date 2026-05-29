setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_and_treatment_separated")
getwd()

head(deseq_res$AM14transfer_PL23)
head(deseq_res$AM14transfer_R848)
head(deseq_res$B18transfer)
head(deseq_res$AM14MRLlpr)

vp_data = list(
  PL23_2DG_v_Ctrl = deseq_res$AM14transfer_PL23[, c(2, 8, 12)],
  R848_2DG_v_Ctrl = deseq_res$AM14transfer_R848[, c(2, 8, 12)],
  AM14MRLlpr = deseq_res$AM14MRLlpr[, c(2, 8, 12)],
  B18 = deseq_res$B18transfer[, c(2, 8, 12)])
head(vp_data$AM14MRLlpr)

search_for_gene <- function(gene_names, lfc_cutoff, ignore_case){
  if(missing(lfc_cutoff)) lfc_cutoff <- 0.6
  if(missing(ignore_case)) ignore_case <- FALSE
  
  for (x in 1:length(gene_names)) {
    temp_res <- gene_IDs_full[grepl(gene_names[x],
                                    gene_IDs_full$external_gene_name,
                                    ignore.case = ignore_case), ]
    
    if(x == 1) genelist <- temp_res
    else genelist <- full_join(genelist, temp_res)
  }
  
  genelist <- genelist$external_gene_name
  
  PL23_res <- deseq_res$AM14transfer_PL23[deseq_res$AM14transfer_PL23$external_gene_name %in% genelist, ]
  PL23_res$comparison <- "PL2-3"
  PL23_res$order <- 1
  
  MRLlpr_res <- deseq_res$AM14MRLlpr[deseq_res$AM14MRLlpr$external_gene_name %in% genelist, ]
  MRLlpr_res$comparison <- "MRL/lpr"
  MRLlpr_res$order <- 2
  
  R848_res <- deseq_res$AM14transfer_R848[deseq_res$AM14transfer_R848$external_gene_name %in% genelist, ]
  R848_res$comparison <- "R848"
  R848_res$order <- 3
  
  B18_res <- deseq_res$B18transfer[deseq_res$B18transfer$external_gene_name %in% genelist, ]
  B18_res$comparison <- "B1-8"
  B18_res$order <- 4
  
  full_res <- full_join(PL23_res, MRLlpr_res)
  full_res <- full_join(full_res, R848_res)
  full_res <- full_join(full_res, B18_res)
  full_res <- full_res[order(full_res$external_gene_name, full_res$order), ]
  full_res$log2FoldChange <- round(full_res$log2FoldChange, digits = 3)
  
  full_res$signif <- ifelse(
    test = full_res$padj <= 0.05, 
    yes = ifelse(
      test = full_res$log2FoldChange >= lfc_cutoff, 
      yes = "Up", 
      no = ifelse(
        test = full_res$log2FoldChange <= -1*lfc_cutoff, 
        yes = "Down", 
        no = "-"
      )
    ),
    no = "-"
  )
  
  full_res <- full_res[, c(2, 13, 15, 8, 12)]
  return(full_res)
}

search_for_gene(c("Tbx21", "Cr2"))
search_for_gene("Tbx21")

search_for_gene("Il2")

