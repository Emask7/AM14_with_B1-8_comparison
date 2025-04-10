QC_heatmaps <- function(dds, filename_start, plot_title, col_factors){
  # Transform data -------------------------------------------------------------
    vsd <- vst(dds)
    rld <- rlog(dds)
    ntd <- normTransform(dds)
  
  # Heatmap of count matrix ----------------------------------------------------
    select <- order(rowMeans(counts(dds,normalized=TRUE)), 
                    decreasing=TRUE)[1:20]
    df <- as.data.frame(colData(dds)[, c("Treatment", "Cohort")])
    
    png(filename = stri_join(c("QC_results/Heatmaps/", filename_start,
                               " - Normalized Counts Transformation.png"),
                             collapse = ""),
        width = 1200, height = 1200, units = "px", pointsize = 10, res = 200,
        bg = "white", family = "", # type = "windows", 
        symbolfamily="default")
    pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, 
             cluster_cols=TRUE, annotation_col=df,
             labels_col = colData(dds)$Label_Name,
             main = stri_join(c(plot_title, "(Normalized Counts Transformation)"), collapse = "\n"))
    dev.off()
    
    png(filename = stri_join(c("QC_results/Heatmaps/", filename_start,
                               " - Variance Stabilizing Transformation.png"),
                             collapse = ""),
        width = 1200, height = 1200, units = "px", pointsize = 10, res = 200,
        bg = "white", family = "", # type = "windows", 
        symbolfamily="default")
    pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
             cluster_cols=TRUE, annotation_col=df,
             labels_col = colData(dds)$Label_Name,
             main = stri_join(c(plot_title, "(Variance Stabilizing Transformation)"), collapse = "\n"))
    dev.off()

    png(filename = stri_join(c("QC_results/Heatmaps/", filename_start,
                               " - Regularized Log Transformation.png"),
                             collapse = ""),
        width = 1200, height = 1200, units = "px", pointsize = 10, res = 200,
        bg = "white", family = "", # type = "windows", 
        symbolfamily="default")
    pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
             cluster_cols=TRUE, annotation_col=df,
             labels_col = colData(dds)$Label_Name,
             main = stri_join(c(plot_title, "(Regularized Log Transformation)"), collapse = "\n"))

    dev.off()
    
  # Heatmap of sample-to-sample distances --------------------------------------
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$Label_Name)
    colnames(sampleDistMatrix) <- NULL
    
    png(filename = stri_join(c("QC_results/Heatmaps/", filename_start,
                               " - Sample-to-Sample Distances.png"),
                             collapse = ""),
        width = 1200, height = 1200, units = "px", pointsize = 10, res = 200, 
        bg = "white", family = "", # type = "windows", 
        symbolfamily="default")
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
             main = stri_join(c(plot_title, "(Sample-to-Sample Distances)"), collapse = "\n"))
    dev.off()
}

QC_PCAplot <- function(dds, filename_start, plot_title, batch_effect){
  vsd <- vst(dds)
  if(batch_effect == FALSE){
    pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Cohort"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    png(filename = stri_join(c("QC_results/PCA_plots/", filename_start, ".png"),
                             collapse = ""),
        width = 1500, height = 1500, units = "px", pointsize = 10, res = 200,
        bg = "white", family = "", symbolfamily="default")
    ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Cohort)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() +
      labs(title = plot_title)
  } else if(batch_effect == TRUE){
    mat <- assay(vsd)
    mm <- model.matrix(~Treatment, colData(vsd))
    mat <- removeBatchEffect(mat, batch=vsd$Cohort, design=mm)
    assay(vsd) <- mat
    pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Cohort"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))

    png(filename = stri_join(c("QC_results/PCA_plots/", filename_start, ".png"),
                             collapse = ""),
        width = 1500, height = 1500, units = "px", pointsize = 10, res = 200,
        bg = "white", family = "", symbolfamily="default")
    ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Cohort)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() +
      labs(title = plot_title)
  } else print("batch_effect must be either TRUE or FALSE")
}

# summary_v2 <- function(res, comparison, experiment, p_cutoff, lfc_cutoff){
#   if(missing(p_cutoff)) p_cutoff <- 0.05
#   if(missing(lfc_cutoff)) lfc_cutoff <- 1
#   # if(missing(title)) title <- "# Genes"
#   
#   up_string <- stri_join(c("Up (LFC >= ", lfc_cutoff, "p-adj <= ", p_cutoff, ")"), collapse = "")
#   down_string <- stri_join(c("Down (LFC <= -", lfc_cutoff, "p-adj <= ", p_cutoff,  ")"), collapse = "")
#   
#   DEG_summary <- data.frame(c("Experiment", up_string, down_string),
#                             c(experiment,
#                               nrow(subset(res, res$padj <= p_cutoff & res$log2FoldChange >= lfc_cutoff)),
#                               nrow(subset(res, res$padj <= p_cutoff & res$log2FoldChange <= (-1*lfc_cutoff)))))
#   
#   # col1_name <- stri_join(c("padj <= ", p_cutoff), collapse = "")
#   colnames(DEG_summary) <- c("Comparison", comparison)
#   DEG_summary
# }

summary_wrapper <- function(res, comparison, experiment, p_cutoff, lfc_cutoff){
  if(missing(p_cutoff)) p_cutoff <- 0.05
  if(missing(lfc_cutoff)) lfc_cutoff <- 1

  up_string <- stri_join(c("Up (LFC >= ", lfc_cutoff, "p-adj <= ", p_cutoff, ")"), collapse = "")
  down_string <- stri_join(c("Down (LFC <= -", lfc_cutoff, "p-adj <= ", p_cutoff,  ")"), collapse = "")
  
  degs_up <- nrow(subset(res, res$padj <= p_cutoff & res$log2FoldChange >= lfc_cutoff))
  degs_down <- nrow(subset(res, res$padj <= p_cutoff & res$log2FoldChange <= (-1*lfc_cutoff)))
  
  DEG_summary <- data.frame(experiment, comparison, lfc_cutoff, p_cutoff, degs_up, degs_down)
  colnames(DEG_summary) <- c("Experiment", "Comparison", "LFC_cutoff", "P_adj_cutoff", "Up", "Down")
  DEG_summary
}

deg_list_for_venn <- function(res, lfc_cutoff){
  list(up = subset(res, res$padj <= 0.05 & res$log2FoldChange >= lfc_cutoff)[, 1], 
       down = subset(res, res$padj <= 0.05 & res$log2FoldChange <= (-1*lfc_cutoff))[, 1], 
       all = subset(res, res$padj <= 0.05 & abs(res$log2FoldChange) >= lfc_cutoff)[, 1])
}



# write_DEG_CSV <- function(res, lfc_cutoff, file_start){
#   up <- subset(res, res$padj <= 0.05 & res$log2FoldChange >= lfc_cutoff)[, c(1, 2, 4)]
#   down <- subset(res, res$padj <= 0.05 & res$log2FoldChange <= (-1*lfc_cutoff))[, c(1, 2, 4)]
#   all <- subset(res, res$padj <= 0.05 & abs(res$log2FoldChange) >= lfc_cutoff)[, c(1, 2, 4)]
#   
#   if(!is.null(file_start)){
#     write.table(up, 
#                 file = stri_join(c("Output/Gene_Lists/", file_start, "_Up.csv"), collapse = ""),
#                 sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
#     write.table(down, 
#                 file = stri_join(c("Output/Gene_Lists/", file_start, "_Down.csv"), collapse = ""),
#                 sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
#     write.table(all,
#                 file = stri_join(c("Output/Gene_Lists/", file_start, "_All.csv"), collapse = ""),
#                 sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
#   }
#   
#   return(list(up = up, down = down, all = all))
# }

# write_sig_LFCs <- function(res, lfc_cutoff, ID_type, file_start){
#   dat <- subset(res, res$padj <= 0.05 & abs(res$log2FoldChange) >= lfc_cutoff)
#   
#   if(missing(ID_type)) print("specify ID type: ensembl_gene_id, external_gene_name, or entrezgene")
#   else if(ID_type == "ensembl_gene_id") dat <- dat[, c(1, 6)]
#   else if(ID_type == "external_gene_name") dat <- dat[, c(2, 6)]
#   else if(ID_type == "entrezgene") dat <- dat[, c(4, 6)]
#   else print("specify ID type: ensembl_gene_id, external_gene_name, or entrezgene")
#   
#   colnames(dat) <- c("#", "LogFoldChange")
#   
#   if(!is.null(file_start)){
#     write.table(dat,
#                 file = stri_join(c("Output/Gene_Lists/", file_start, "_LFC_values.csv"), collapse = ""),
#                 sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
#   }
#   
#   dat
# }

