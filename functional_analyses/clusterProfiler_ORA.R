# using the packaged "ReactomePA"

# extract EntrezID for Differentially Expressed Genes
de_genes <- sig$ENTREZID
head(de_genes)

# perform pathway enrichment analysis
pathway2 <- enrichPathway(gene=de_genes, pvalueCutoff = 0.05, readable=TRUE)
head(pathway2)
