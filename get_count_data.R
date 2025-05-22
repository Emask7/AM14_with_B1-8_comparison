
find_ensembl_ID <- function(search_name) {
  dat <- dplyr::filter(gene_IDs$AM14trans, external_gene_name == search_name)
  dat$ensembl_gene_id
}

get_counts <- function(gene_name){
  ID <- find_ensembl_ID(gene_name)
  
  res <- counts(dds_AM14trans, normalized = TRUE)
  
  res <- res[rownames(res) %in% ID, ]
  # res
  
  df <- data.frame(res[1:5], res[6:10], res[11:15], res[16:20])
  colnames(df) <- c("PL2-3", "PL2-3+2DG", "R848", "R848+2DG")
  rownames(df) <- c("a", "b", "c", "d", "e")
  df
}

# Complement genes
  get_counts("C1qa")
  get_counts("C1qb")
  get_counts("C1qc")
  get_counts("Fcna")
  get_counts("C6")
  get_counts("C4b")
  get_counts("C3")
  get_counts("Cfh")
  get_counts("Ighg2c")
  get_counts("Ighm")
  get_counts("Ighg2b")
  

# TLR signaling
  get_counts("Unc93b1")
  get_counts("Tlr9")



get_counts("Ptpn22")
get_counts("Treml4")
get_counts("Rab7b")
get_counts("Colec12")
get_counts("Cxcr4")
get_counts("Gpr55")
