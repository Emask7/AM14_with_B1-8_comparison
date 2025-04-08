summary_v2 <- function(res, title, p_cutoff, lfc_cutoff){
  if(missing(p_cutoff)) p_cutoff <- 0.05
  if(missing(lfc_cutoff)) lfc_cutoff <- 0.6
  if(missing(title)) title <- "# Genes"
  
  up_string <- stri_join(c("Up (LFC >= ", lfc_cutoff, ")"), collapse = "")
  down_string <- stri_join(c("Down (LFC <= -", lfc_cutoff, ")"), collapse = "")
  
  DEG_summary <- data.frame(c(up_string, down_string),
                            c(nrow(subset(res, res$padj <= p_cutoff & res$log2FoldChange >= lfc_cutoff)),
                              nrow(subset(res, res$padj <= p_cutoff & res$log2FoldChange <= (-1*lfc_cutoff)))))
  
  col1_name <- stri_join(c("padj <= ", p_cutoff), collapse = "")
  colnames(DEG_summary) <- c(col1_name, title)
  DEG_summary
}

write_DEG_CSV <- function(res, lfc_cutoff, file_start){
  up <- subset(res, res$padj <= 0.05 & res$log2FoldChange >= lfc_cutoff)[, c(1, 2, 4)]
  down <- subset(res, res$padj <= 0.05 & res$log2FoldChange <= (-1*lfc_cutoff))[, c(1, 2, 4)]
  all <- subset(res, res$padj <= 0.05 & abs(res$log2FoldChange) >= lfc_cutoff)[, c(1, 2, 4)]
  
  if(!is.null(file_start)){
    write.table(up, 
                file = stri_join(c("Output/Gene_Lists/", file_start, "_Up.csv"), collapse = ""),
                sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(down, 
                file = stri_join(c("Output/Gene_Lists/", file_start, "_Down.csv"), collapse = ""),
                sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(all,
                file = stri_join(c("Output/Gene_Lists/", file_start, "_All.csv"), collapse = ""),
                sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
  
  return(list(up = up, down = down, all = all))
}

write_sig_LFCs <- function(res, lfc_cutoff, ID_type, file_start){
  dat <- subset(res, res$padj <= 0.05 & abs(res$log2FoldChange) >= lfc_cutoff)
  
  if(missing(ID_type)) print("specify ID type: ensembl_gene_id, external_gene_name, or entrezgene")
  else if(ID_type == "ensembl_gene_id") dat <- dat[, c(1, 6)]
  else if(ID_type == "external_gene_name") dat <- dat[, c(2, 6)]
  else if(ID_type == "entrezgene") dat <- dat[, c(4, 6)]
  else print("specify ID type: ensembl_gene_id, external_gene_name, or entrezgene")
  
  colnames(dat) <- c("#", "LogFoldChange")
  
  if(!is.null(file_start)){
    write.table(dat,
                file = stri_join(c("Output/Gene_Lists/", file_start, "_LFC_values.csv"), collapse = ""),
                sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
  
  dat
}

