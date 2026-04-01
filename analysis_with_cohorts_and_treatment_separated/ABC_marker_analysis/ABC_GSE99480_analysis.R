setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_and_treatment_separated/ABC_marker_analysis/")
# setwd("./analysis_with_cohorts_and_treatment_separated/ABC_marker_analysis/")
getwd()

# Load the necessary libraries
library(fgsea)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(ggrepel)

gsea_format <- function(dds_res){
  res_filtered <- dds_res[dds_res$external_gene_name != "", ]
  res_filtered <- res_filtered[res_filtered$pvalue <= 0.05, ]
  genelist <- res_filtered$log2FoldChange
  names(genelist) <- res_filtered$external_gene_name
  genelist <- sort(genelist, decreasing = TRUE)
  genelist <- genelist[!duplicated(names(genelist))]
  return(genelist)
}

compile_gsea <- function(res_list, names){
  for (res in 1:length(res_list)) {
    temp_df <- data.frame(res_list[res])
    temp_df <- data.frame(Comparison = rep(names[res], nrow(temp_df)), temp_df)
    if(res == 1) compiled_df <- temp_df
    else compiled_df <- full_join(compiled_df, temp_df)
  }
  return(compiled_df)
}

# --- Load your custom gene sets ---
# read.table reads the single-column file into a data frame,
# we then extract that column as a vector of gene names.
abc_up_genes <- read.table("../../public_GSE_datasets/GSE_lists/GSE99480_ABC_up_regulated_genes.txt", header = FALSE)$V1
abc_down_genes <- read.table("../../public_GSE_datasets/GSE_lists/GSE99480_ABC_down_regulated_genes.txt", header = FALSE)$V1

abc_gmt <- data.frame(
  "term" = c(rep("ABC Down, FoB Up", length(abc_down_genes)), rep("ABC Up, FoB Down", length(abc_up_genes))),
  "gene" = c(abc_down_genes, abc_up_genes)
)
head(abc_gmt)


# --- Prepare your ranked gene list ---
pl23_genelist <- gsea_format(deseq_res$AM14transfer_PL23)
MRLlpr_genelist <- gsea_format(deseq_res$AM14MRLlpr)
r848_genelist <- gsea_format(deseq_res$AM14transfer_R848)
np_genelist <- gsea_format(deseq_res$B18transfer)

pl23_ABC_gsea <- GSEA(pl23_genelist, TERM2GENE = abc_gmt, pvalueCutoff = 1)
MRLlpr_ABC_gsea <- GSEA(MRLlpr_genelist, TERM2GENE = abc_gmt, pvalueCutoff = 1)
r848_ABC_gsea <- GSEA(r848_genelist, TERM2GENE = abc_gmt, pvalueCutoff = 1)
np_ABC_gsea <- GSEA(np_genelist, TERM2GENE = abc_gmt, pvalueCutoff = 1)

ABC_vs_FoB_summary <- compile_gsea(res_list = list(pl23_ABC_gsea[, 2:11],
                                                   MRLlpr_ABC_gsea[, 2:11],
                                                   r848_ABC_gsea[, 2:11],
                                                   np_ABC_gsea[, 2:11]), 
                                   names = c("PL2-3", "MRLlpr", "R848", "NP"))
ABC_vs_FoB_summary
write.table(ABC_vs_FoB_summary, file = "ABC_GSE99480_summary.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)


# p1 <- gseaplot2(pl23_ABC_gsea, geneSetID = 1:2, subplots = 1:3, title = "PL2-3 +/- 2DG", color = c("#dd4124", "#0f85a0"), rel_heights = c(2, 0.5, 1))
# p2 <- gseaplot2(MRLlpr_ABC_gsea, geneSetID = 1:2, subplots = 1:2, title = "AM14 MRL/lpr +/- 2DG", color = c("#dd4124", "#0f85a0"), rel_heights = c(5, 0.5))
# p3 <- gseaplot2(r848_ABC_gsea, geneSetID = 1:2, subplots = 1:2, title = "R848 +/- 2DG", color = c("#dd4124", "#0f85a0"), rel_heights = c(5, 0.5))
# p4 <- gseaplot2(np_ABC_gsea, geneSetID = 1:2, subplots = 1:2, title = "NP +/- 2DG", color = c("#dd4124", "#0f85a0"), rel_heights = c(5, 0.5))
# 
# plot_list(p1, p2, p3, p4, ncol=2, tag_levels = 'A') 


genes <- c("Tbx21", "Itgax", "Itgam", "Cr2", "Cxcr3", "Zeb2",
           "Irf1", "Irf4", "Irf8", "Irf7", "Tlr7", "Tlr9", 
           "C1qc", "C1qb", "Ly6c1", "Batf3", "Ccl4", "Lck", "Gcsam", "Zbtb32",
           overlapping_DEGs_list)
abc_up_labels <- genes[genes %in% abc_up_genes]
abc_up_labels
abc_down_labels <- genes[genes %in% abc_down_genes]
abc_down_labels


# p1 <- gseaplot2(pl23_ABC_gsea, geneSetID = 1, subplots = 1:2, title = "PL2-3 +/- 2DG (ABC Up, FoB Down)", color = c("#dd4124"), rel_heights = c(2, 0.5), base_size = 15)
# p1[[1]] <- p1[[1]] + geom_gsea_gene(gene = abc_up_labels, geom = ggrepel::geom_label_repel, size = 5)
# # p1
# 
# p2 <- gseaplot2(MRLlpr_ABC_gsea, geneSetID = 1, subplots = 1:2, title = "AM14 MRL/lpr +/- 2DG (ABC Up, FoB Down)", color = c("#dd4124"), rel_heights = c(2, 0.5), base_size = 15)
# p3 <- gseaplot2(r848_ABC_gsea, geneSetID = 1, subplots = 1:2, title = "R848 +/- 2DG (ABC Up, FoB Down)", color = c("#dd4124"), rel_heights = c(2, 0.5), base_size = 15)
# p4 <- gseaplot2(np_ABC_gsea, geneSetID = 1, subplots = 1:2, title = "NP +/- 2DG (ABC Up, FoB Down)", color = c("#dd4124"), rel_heights = c(2, 0.5), base_size = 15)
# 
# p5 <- gseaplot2(pl23_ABC_gsea, geneSetID = 2, subplots = 1:2, title = "PL2-3 +/- 2DG (ABC Down, FoB Up)", color = c("#0f85a0"), rel_heights = c(2, 0.5), base_size = 15)
# p6 <- gseaplot2(MRLlpr_ABC_gsea, geneSetID = 2, subplots = 1:2, title = "AM14 MRL/lpr +/- 2DG (ABC Down, FoB Up)", color = c("#0f85a0"), rel_heights = c(2, 0.5), base_size = 15)
# p7 <- gseaplot2(r848_ABC_gsea, geneSetID = 2, subplots = 1:2, title = "R848 +/- 2DG (ABC Down, FoB Up)", color = c("#0f85a0"), rel_heights = c(2, 0.5), base_size = 15)
# p8 <- gseaplot2(np_ABC_gsea, geneSetID = 2, subplots = 1:2, title = "NP +/- 2DG (ABC Down, FoB Up)", color = c("#0f85a0"), rel_heights = c(2, 0.5), base_size = 15)


p1 <- gseaplot2(pl23_ABC_gsea, geneSetID = 1, subplots = 1:2, title = "PL2-3 +/- 2DG\nNES: -3.040     p.adjust: 2.00E-10", color = c("#dd4124"), rel_heights = c(2.25, 0.25), base_size = 15)
p1[[1]] <- p1[[1]] + geom_gsea_gene(gene = abc_up_labels, geom = ggrepel::geom_label_repel, size = 5, max.overlaps = Inf, force = 5)
# p1

p2 <- gseaplot2(MRLlpr_ABC_gsea, geneSetID = 1, subplots = 1:2, title = "AM14 MRL/lpr +/- 2DG\nNES: -2.180     p.adjust: 6.36E-10", color = c("#dd4124"), rel_heights = c(2.25, 0.25), base_size = 15)
p3 <- gseaplot2(r848_ABC_gsea, geneSetID = 1, subplots = 1:2, title = "R848 +/- 2DG\nNES: -2.503     p.adjust: 2.00E-10", color = c("#dd4124"), rel_heights = c(2.25, 0.25), base_size = 15)
p4 <- gseaplot2(np_ABC_gsea, geneSetID = 1, subplots = 1:2, title = "NP +/- 2DG\nNES: -1.370     p.adjust: 9.14E-02", color = c("#dd4124"), rel_heights = c(2.25, 0.25), base_size = 15)

p5 <- gseaplot2(pl23_ABC_gsea, geneSetID = 2, subplots = 1:2, title = "PL2-3 +/- 2DG\nNES: -1.285     p.adjust: 5.16E-02", color = c("#0f85a0"), rel_heights = c(2.25, 0.25), base_size = 15)
p6 <- gseaplot2(MRLlpr_ABC_gsea, geneSetID = 2, subplots = 1:2, title = "AM14 MRL/lpr +/- 2DG\nNES: -1.224     p.adjust: 7.77E-02", color = c("#0f85a0"), rel_heights = c(2.25, 0.25), base_size = 15)
p7 <- gseaplot2(r848_ABC_gsea, geneSetID = 2, subplots = 1:2, title = "R848 +/- 2DG\nNES: -1.821     p.adjust: 3.34E-04", color = c("#0f85a0"), rel_heights = c(2.25, 0.25), base_size = 15)
p8 <- gseaplot2(np_ABC_gsea, geneSetID = 2, subplots = 1:2, title = "NP +/- 2DG\nNES: 0.986     p.adjust: 5.12E-01", color = c("#0f85a0"), rel_heights = c(2.25, 0.25), base_size = 15)



p2[[1]] <- p2[[1]] + geom_gsea_gene(gene = abc_up_labels, geom = ggrepel::geom_label_repel, size = 5, max.overlaps = Inf, force = 5)
p3[[1]] <- p3[[1]] + geom_gsea_gene(gene = abc_up_labels, geom = ggrepel::geom_label_repel, size = 5, max.overlaps = Inf, force = 5)
p4[[1]] <- p4[[1]] + geom_gsea_gene(gene = abc_up_labels, geom = ggrepel::geom_label_repel, size = 5, max.overlaps = Inf, force = 5)

p5[[1]] <- p5[[1]] + geom_gsea_gene(gene = abc_down_labels, geom = ggrepel::geom_label_repel, size = 5, max.overlaps = Inf, force = 5)
p6[[1]] <- p6[[1]] + geom_gsea_gene(gene = abc_down_labels, geom = ggrepel::geom_label_repel, size = 5, max.overlaps = Inf, force = 5)
p7[[1]] <- p7[[1]] + geom_gsea_gene(gene = abc_down_labels, geom = ggrepel::geom_label_repel, size = 5, max.overlaps = Inf, force = 5)
p8[[1]] <- p8[[1]] + geom_gsea_gene(gene = abc_down_labels, geom = ggrepel::geom_label_repel, size = 5, max.overlaps = Inf, force = 5)

plot_list(p1, p2, p3, p4, ncol=2, tag_levels = 'A', tag_size = 30) 
plot_list(p1, p2, ncol=1, tag_levels = 'B') 

plot_list(p1, p2, p5, p6, ncol=2, tag_levels = 'A', tag_size = 30) 
ggsave(filename = "PL2-3 and AM14 MRLlpr ABC up and down.png", width = 10, height = 8, units = "in", dpi = 600)

plot_list(p1, p2, p5, p6, ncol=2, tag_levels = 'A', tag_size = 30) 
ggsave(filename = "R848 and NP ABC up and down.png", width = 10, height = 8, units = "in", dpi = 600)

plot_list(p1, p5, p2, p6, p3, p7, p4, p8, ncol=2, tag_levels = NULL) 
ggsave(filename = "All Comparisons - ABC up and down - vertical2.png", width = 12, height = 15, units = "in", dpi = 600)

plot_list(p1, p2, p3, p4, p5, p6, p7, p8, ncol=4, tag_levels = 'A', tag_size = 30)
ggsave(filename = "All Comparisons - ABC up and down - horizontal.png", width = 22, height = 10, units = "in", dpi = 600)
ggsave(filename = "All Comparisons - ABC up and down - horizontal2.png", width = 22, height = 9, units = "in", dpi = 600)











gseaplot2(pl23_ABC_gsea, geneSetID = 1, title = pl23_ABC_gsea$Description[1])
gseaplot2(pl23_ABC_gsea, geneSetID = 1:2, pvalue_table = TRUE)
gseaplot2(pl23_ABC_gsea, geneSetID = 1:2)
gseaplot2(pl23_ABC_gsea, geneSetID = 1, pvalue_table = TRUE, pvalue_table_columns = c("ID", "NES", "p.adjust"))

p1 <- gseaplot2(pl23_ABC_gsea, geneSetID = 1, subplots = 1, title = "ABC up, FoB down")
p2 <- gseaplot2(pl23_ABC_gsea, geneSetID = 2, subplots = 1:2, title = "ABC down, FoB up")
plot_list(p1, p2, ncol=1, tag_levels = 'A')

p3 <- gseaplot2(pl23_ABC_gsea, geneSetID = 1, title = "ABC up, FoB down", pvalue_table = TRUE, pvalue_table_rownames = NULL)
p3[[1]] <- p3[[1]] + geom_gsea_gene(genes, geom = geom_label)
p3


# --- Run the fgsea analysis ---
# Create the named list required by fgsea
custom_gene_sets <- list(
  "ABC_up_vs_FoB" = abc_up_genes,
  "ABC_down_vs_FoB" = abc_down_genes
)
head(custom_gene_sets$ABC_up_vs_FoB)

set.seed(42) # for reproducibility

pl23_fgsea_results <- fgsea(
  pathways = custom_gene_sets, stats = pl23_genelist, 
  minSize = 15, maxSize = 500
)

# View the results table, ordered by significance
print("GSEA Results Table:")
pl23_fgsea_results %>%
  arrange(padj) %>%
  as.data.frame()

# Plot the enrichment score for the "ABC_up_regulated" gene set
plotEnrichment(custom_gene_sets[["ABC_up_vs_FoB"]], pl23_genelist) +
  labs(title = "Enrichment of ABC Up-regulated Genes")

# Plot the enrichment score for the "ABC_down_regulated" gene set
plotEnrichment(custom_gene_sets[["ABC_down_vs_FoB"]], pl23_genelist) +
  labs(title = "Enrichment of ABC Down-regulated Genes")


# Create a summary plot of the GSEA results
plotGseaTable(
  pathways = custom_gene_sets[pl23_fgsea_results$pathway],
  stats = pl23_genelist,
  fgseaRes = pl23_fgsea_results,
  gseaParam = 0.5
)
