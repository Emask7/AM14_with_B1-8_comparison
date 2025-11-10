setwd("./analysis_with_cohorts_and_treatment_separated/GSVA_analysis")
getwd()

genesets <- read.xlsx("Bcell_subset.human.xlsx", startRow = 2)
genesets <- as.list(genesets)
genesets <- lapply(genesets, function(x) x[!is.na(x)])
head(genesets)

format_gsva_counts <- function(dds){
  cts <- counts(dds, normalized = TRUE)
  temp_colnames <- colnames(cts)
  cts <- data.frame(rownames(cts), cts)
  colnames(cts) <- c("ensembl_gene_id", temp_colnames)
  cts <- right_join(gene_IDs_full, cts, by = "ensembl_gene_id")
  cts <- cts[cts$hgnc_symbol != "", ]
  cts <- cts %>%
    mutate(RowSum = apply(cts, 1, function(x) rowSums(cts[, 7:ncol(cts)]))[, 1]) %>%
    group_by(hgnc_symbol) %>%
    slice_max(order_by = RowSum, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  cts <- data.frame(cts)
  rownames(cts) <- cts$hgnc_symbol
  cts <- cts[ , 7:(ncol(cts)-1)]
  cts <- as.matrix(cts)
  return(cts)
}

pl23_counts <- format_gsva_counts(dds_AM14trans_PL23)
head(pl23_counts)
nrow(pl23_counts)

# normcounts <- list("PL23" = matrix(
#                    "MRLlpr" = counts(dds_AM14MRLlpr, normalized = TRUE),
#                    "R848" = counts(dds_AM14trans_R848, normalized = TRUE),
#                    "NP" = counts(dds_B18trans, normalized = TRUE))
# head(normcounts$PL23)

# GSVA parameter object
  pl23_gsvaPar <- gsvaParam(pl23_counts, genesets)
  pl23_gsvaPar
  
# Calculate GSVA enrichment scores
  pl23_gsvaPar.es <- gsva(pl23_gsvaPar, verbose=FALSE)
  dim(pl23_gsvaPar.es)
  pl23_gsvaPar.es[1:5, 1:5]
  
  res <- igsva()

# p <- 10000 ## number of genes
# n <- 30    ## number of samples
# ## simulate expression values from a standard Gaussian distribution
# X <- matrix(rnorm(p*n), nrow=p,
#             dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
# X[1:5, 1:5]
# ## sample gene set sizes
# gs <- as.list(sample(10:100, size=100, replace=TRUE))
# ## sample gene sets
# gs <- lapply(gs, function(n, p) paste0("g", sample(1:p, size=n, replace=FALSE)), p)
# names(gs) <- paste0("gs", 1:length(gs))
# head(gs)



