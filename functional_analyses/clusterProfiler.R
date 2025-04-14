### clusterProfiler
## https://guangchuangyu.github.io/software/clusterProfiler/
## https://yulab-smu.top/biomedical-knowledge-mining-book/index.html

# GO Enrichment Analysis -------------------------------------------------------
## background info: http://geneontology.org/docs/ontology-documentation/
enrichGO_wrapper <- function(dds_res){
  # Create background dataset for hypergeometric testing using all genes tested for significance in the results
    background <- as.character(dds_res$ensembl_gene_id)
  
  # Extract significant results
    sig <- dplyr::filter(dds_res, padj < 0.05)
    sig <- as.character(sig$ensembl_gene_id)
  
  # Run GO enrichment analysis 
    eGO <- enrichGO(gene = sig, universe = background, keyType = "ENSEMBL",
                    OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, readable = TRUE)
    
  return(list(eGO = eGO, summary = data.frame(eGO)))
}

enrichGO_results <- list(
    AM14trans = list(
      PL23_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl),
      R848_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl),
      PL23_vs_R848 = enrichGO_wrapper(deseq_res$AM14transfer$PL23_vs_R848)),
    B18trans = enrichGO_wrapper(deseq_res$B18transfer),
    AM14MRLlpr = enrichGO_wrapper(deseq_res$AM14MRLlpr)
)

head(enrichGO_results$AM14trans$PL23_2DG_v_Ctrl)

wb <- createWorkbook("Output/Functional_analyses/GO_enrichment.xlsx")
addWorksheet(wb, "AM14 - PL23_2DG_vs_Ctrl")
addWorksheet(wb, "AM14 - R848_2DG_vs_Ctrl")
addWorksheet(wb, "AM14 - PL23_vs_R848")
addWorksheet(wb, "B1-8 - 2DG_vs_Ctrl")
addWorksheet(wb, "AM14 MRLlpr - 2DG_vs_Ctrl")
writeData(wb, "AM14 - PL23_2DG_vs_Ctrl", enrichGO_results$AM14trans$PL23_2DG_v_Ctrl$summary)
writeData(wb, "AM14 - R848_2DG_vs_Ctrl", enrichGO_results$AM14trans$R848_2DG_v_Ctrl$summary)
writeData(wb, "AM14 - PL23_vs_R848", enrichGO_results$AM14trans$PL23_vs_R848$summary)
writeData(wb, "B1-8 - 2DG_vs_Ctrl", enrichGO_results$B18trans$summary)
writeData(wb, "AM14 MRLlpr - 2DG_vs_Ctrl", enrichGO_results$AM14MRLlpr$summary)
saveWorkbook(wb, "Output/Functional_analyses/GO_enrichment.xlsx", overwrite = TRUE)
rm(wb)

### Visualizing clusterProfiler results
# plot(barplot(enrichGO_results$AM14trans$PL23_2DG_v_Ctrl$eGO, showCateGOry=25))
# goplot(enrichGO_results$AM14trans$PL23_2DG_v_Ctrl$eGO)

# Dotplot Wrapper Function -----------------------------------------------------
dotplot_wrapper <- function(eGO_res, n, plot_title, file_name, h){
  if(missing(n)) n <- 25
  if(missing(h)) h <- 2200
  if(nrow(eGO_res) < n) n <- nrow(eGO_res)
  
  eGO_top <- eGO_res[order(eGO_res$p.adjust, decreasing = FALSE), ]
  eGO_top <- eGO_top[1:n, ]$Description
  
  if(!missing(file_name)){
    if(file_name != FALSE & !is.null(file_name)){
      png(filename = stri_join(c("Output/Functional_analyses/", file_name, ".png"),
                               collapse = ""),
          width = 1200, height = h, units = "px", pointsize = 10, res = 200,
          bg = "white", family = "", symbolfamily="default")
      dotplot(eGO_res, showCategory = eGO_top, title = plot_title)
    } else dotplot(eGO_res, showCategory = eGO_top, title = plot_title)
  } else dotplot(eGO_res, showCategory = eGO_top, title = plot_title)
  # NOTE: GeneRatio = # genes in your input list associated with the given GO term/total # of input genes.
}

dotplot_wrapper(enrichGO_results$AM14trans$PL23_2DG_v_Ctrl$eGO, 25, 
                "AM14 Adoptive Transfer:\nPL2-3 + 2DG vs PL2-3",
                "enrichGO_plots/AM14 Transfer - PL2-3 + 2DG vs PL2-3")
dev.off()

dotplot_wrapper(enrichGO_results$AM14trans$R848_2DG_v_Ctrl$eGO, 25, 
                "AM14 Adoptive Transfer:\nR848 + 2DG vs R848",
                "enrichGO_plots/AM14 Transfer - R848 + 2DG vs R848", 500)
dev.off()

dotplot_wrapper(enrichGO_results$AM14trans$PL23_vs_R848$eGO, 25, 
                "AM14 Adoptive Transfer:\nPL2-3 vs R848",
                "enrichGO_plots/AM14 Transfer - PL2-3 vs R848")
dev.off()

dotplot_wrapper(enrichGO_results$B18trans$eGO, 25, 
                "B1-8 Adoptive Transfer:\nNP + 2DG vs NP",
                "enrichGO_plots/B1-8 Transfer - NP + 2DG vs NP", 600)
dev.off()

dotplot_wrapper(enrichGO_results$AM14MRLlpr$eGO, 25, 
                "AM14 MRLlpr:\n2DG vs Control",
                "enrichGO_plots/AM14 MRLlpr")
dev.off()

# ## Enrichmap clusters the 25 most significant (by padj) GO terms to visualize relationships between terms
# eGO2 <- pairwise_termsim(enrichGO_results$AM14trans$PL23_2DG_v_Ctrl$eGO) 
# emapplot(eGO2)





# KEGG pathway over-representation analysis ------------------------------------
# see https://www.genome.jp/kegg/ for information about KEGG
enrichKEGG_wrapper <- function(dds_res){
  sig <- dplyr::filter(dds_res, padj < 0.05)
  sig <- as.character(sig$entrezgene)
  enrichKEGG(gene = sig, organism = 'mmu', pvalueCutoff = 0.05)
}

enrichKEGG_results <- list(
  AM14trans = list(
    PL23_2DG_v_Ctrl = enrichKEGG_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl),
    R848_2DG_v_Ctrl = enrichKEGG_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl),
    PL23_vs_R848 = enrichKEGG_wrapper(deseq_res$AM14transfer$PL23_vs_R848)),
  B18trans = enrichKEGG_wrapper(deseq_res$B18transfer),
  AM14MRLlpr = enrichKEGG_wrapper(deseq_res$AM14MRLlpr)
)

nrow(enrichKEGG_results$AM14MRLlpr)

# ## visualize enriched KEGG pathways
# browseKEGG(kk, 'hsa04510')
# 
# lg2FC <- sig$log2FoldChange
# names(lg2FC) <- sig$ENTREZID
# head(lg2FC)
# 
# hsa04510 <- pathview(gene.data  = lg2FC,
#                      pathway.id = "hsa04510",
#                      species    = "hsa",
#                      limit      = list(gene=5, cpd=1))
# 
# # set limit to 3 for optimal color scale - change this as needed


# GO Gene set enrichment analysis ----------------------------------------------
gseGO_wrapper <- function(dds_res){
  sig <- dplyr::filter(dds_res, padj < 0.05)
  sig_ordered <- sig[sig$baseMean > 50,]
  sig_ordered <- sig_ordered[order(-sig_ordered$stat), ]
  gene_list <- sig_ordered$stat
  names(gene_list) <- sig_ordered$ensembl_gene_id
  
  if(length(gene_list) <= 0) return(NULL)
  else {
    gse <- gseGO(gene_list, ont = "BP", keyType = "ENSEMBL", 
                 OrgDb = "org.Mm.eg.db", eps = 1e-300)
    return(list(gse = gse, gse_summary = as.data.frame(gse)))
  }
}

gseGO_results <- list(
  AM14trans = list(
    PL23_2DG_v_Ctrl = gseGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl),
    R848_2DG_v_Ctrl = gseGO_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl),
    PL23_vs_R848 = gseGO_wrapper(deseq_res$AM14transfer$PL23_vs_R848)),
  B18trans = gseGO_wrapper(deseq_res$B18transfer),
  AM14MRLlpr = gseGO_wrapper(deseq_res$AM14MRLlpr)
)
nrow(gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse_summary)
nrow(gseGO_results$AM14trans$R848_2DG_v_Ctrl$gse_summary)
nrow(gseGO_results$AM14trans$PL23_vs_R848$gse_summary)
nrow(gseGO_results$B18trans$gse_summary)
nrow(gseGO_results$AM14MRLlpr$gse_summary)

# visualize results
gseaplot(gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse, geneSetID = 1, 
         title = gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse$Description[1])
gseaplot2(gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse, geneSetID = 1, 
          title = gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse$Description[1])
gseaplot2(gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse, geneSetID = 1:4)
head(gseGO_results$AM14trans$PL23_2DG_v_Ctrl$gse)




# Reactome pathway over-representation analysis --------------------------------
# using ClusterProfiler
gsePathway_wrapper <- function(dds_res){
  sig <- dplyr::filter(dds_res, padj < 0.05)

  sig_up <- sig[sig$log2FoldChange > 0,]
  sig_up = sig_up[order(sig_up[,'entrezgene']),]
  sig_up = sig_up[!duplicated(sig_up$entrezgene),]
  sig_up <- sig_up[order(-sig_up$log2FoldChange), ]
  sig_up <- na.omit(sig_up)

  genes_up <- sig_up$log2FoldChange
  names(genes_up) <- sig_up$entrezgene

  pathway_up <- gsePathway(genes_up, organism = "mouse",
                           pvalueCutoff = 0.2, pAdjustMethod = "BH",
                           verbose = FALSE, scoreType="pos")

  sig_down <- sig[sig$log2FoldChange < 0,]
  sig_down = sig_down[order(sig_down[,'entrezgene']),]
  sig_down = sig_down[!duplicated(sig_down$entrezgene),]
  sig_down <- sig_down[order(-sig_down$log2FoldChange), ]
  sig_down <- na.omit(sig_down)

  genes_down <- sig_down$log2FoldChange
  names(genes_down) <- sig_down$entrezgene

  pathway_down <- gsePathway(genes_down, organism = "mouse",
                             pvalueCutoff = 0.2, pAdjustMethod = "BH",
                             verbose = FALSE, scoreType="pos")

  return(list(up = pathway_up, down = pathway_down))
}

gsePathway_res <- gsePathway_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl)
gsePathway_res
# head(pathway)
# 
# viewPathway("Signaling by GPCR", 
#             readable = TRUE, 
#             foldChange = genes_sorted)


# using the packaged "ReactomePA"

# extract EntrezID for Differentially Expressed Genes
de_genes <- sig$ENTREZID
head(de_genes)

# perform pathway enrichment analysis
pathway2 <- enrichPathway(gene=de_genes, pvalueCutoff = 0.05, readable=TRUE)
head(pathway2)



### Resources for functional analysis

# g:Profiler - http://biit.cs.ut.ee/gprofiler/index.cgi 
# DAVID - http://david.abcc.ncifcrf.gov/tools.jsp 
# clusterProfiler - http://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
# GeneMANIA - http://www.genemania.org/
# GenePattern -  http://www.broadinstitute.org/cancer/software/genepattern/ (need to register)
# WebGestalt - http://bioinfo.vanderbilt.edu/webgestalt/ (need to register)
# AmiGO - http://amigo.geneontology.org/amigo
# ReviGO (visualizing GO analysis, input is GO terms) - http://revigo.irb.hr/ 
# WGCNA - http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork
# GSEA - http://software.broadinstitute.org/gsea/index.jsp
# SPIA - https://www.bioconductor.org/packages/release/bioc/html/SPIA.html
# GAGE/Pathview - http://www.bioconductor.org/packages/release/bioc/html/gage.html
# Reactome Pathway Analysis - https://reactome.org/
# Metascape - https://metascape.org/gp/index.html#/main/step1
# Cytoscape - http://cytoscape.org
# Chea3 - https://maayanlab.cloud/chea3/ (transcription factor enrichment analysis)