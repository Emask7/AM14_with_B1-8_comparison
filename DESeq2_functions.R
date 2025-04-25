import_Rosalind_data <- function(res_file){
  temp_cts <- read.delim(res_file)
  temp_cts <- temp_cts %>%
    filter(!grepl("Igh", external_gene_name)) %>%
    filter(!grepl("Igk", external_gene_name)) # %>%
    # filter(!grepl("7SK", external_gene_name)) %>%
    # filter(!grepl("5_8S_rRNA", external_gene_name)) %>%
    # filter(!grepl("5S_rRNA", external_gene_name))
  temp_cts <- temp_cts[, c(1:4, 12:ncol(temp_cts))]
  temp_cts <- temp_cts[!duplicated(temp_cts), ]
  temp_cts <- temp_cts[!duplicated(temp_cts$ensembl_gene_id), ]
  rownames(temp_cts) <- temp_cts$ensembl_gene_id
  return(temp_cts)
}

QC_heatmaps <- function(dds, filename_start, plot_title){
  # Transform data -------------------------------------------------------------
    vsd <- vst(dds)
    rld <- rlog(dds)
    ntd <- normTransform(dds)
  
  # Heatmap of count matrix ----------------------------------------------------
    select <- order(rowMeans(counts(dds,normalized=TRUE)), 
                    decreasing=TRUE)[1:20]
    df <- as.data.frame(colData(dds)[, c("Treatment", "Cohort")])
    
    if(!is.null(filename_start)) {
      png(filename = stri_join(c("QC_results/Heatmaps/", filename_start,
                                 " - Normalized Counts Transformation.png"),
                               collapse = ""),
          width = 1200, height = 1200, units = "px", pointsize = 10, res = 200,
          bg = "white", family = "", symbolfamily="default")
      pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, 
               cluster_cols=TRUE, annotation_col=df, labels_col = colData(dds)$Label_Name,
               main = stri_join(c(plot_title, "(Normalized Counts Transformation)"), collapse = "\n"))
      dev.off()
      
      png(filename = stri_join(c("QC_results/Heatmaps/", filename_start,
                                 " - Variance Stabilizing Transformation.png"),
                               collapse = ""),
          width = 1200, height = 1200, units = "px", pointsize = 10, res = 200,
          bg = "white", family = "", symbolfamily="default")
      pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
               cluster_cols=TRUE, annotation_col=df, labels_col = colData(dds)$Label_Name,
               main = stri_join(c(plot_title, "(Variance Stabilizing Transformation)"), collapse = "\n"))
      dev.off()
      
      png(filename = stri_join(c("QC_results/Heatmaps/", filename_start,
                                 " - Regularized Log Transformation.png"),
                               collapse = ""),
          width = 1200, height = 1200, units = "px", pointsize = 10, res = 200,
          bg = "white", family = "", symbolfamily="default")
      pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
               cluster_cols=TRUE, annotation_col=df, labels_col = colData(dds)$Label_Name,
               main = stri_join(c(plot_title, "(Regularized Log Transformation)"), collapse = "\n"))
      dev.off()
    } else {
      pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, 
               cluster_cols=TRUE, annotation_col=df, labels_col = colData(dds)$Label_Name,
               main = stri_join(c(plot_title, "(Normalized Counts Transformation)"), collapse = "\n"))

      pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
               cluster_cols=TRUE, annotation_col=df, labels_col = colData(dds)$Label_Name,
               main = stri_join(c(plot_title, "(Variance Stabilizing Transformation)"), collapse = "\n"))

      pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
               cluster_cols=TRUE, annotation_col=df, labels_col = colData(dds)$Label_Name,
               main = stri_join(c(plot_title, "(Regularized Log Transformation)"), collapse = "\n"))
    }
    
    
  # Heatmap of sample-to-sample distances --------------------------------------
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$Label_Name)
    colnames(sampleDistMatrix) <- NULL
    
    if(!is.null(filename_start)) {
      png(filename = stri_join(c("QC_results/Heatmaps/", filename_start,
                                 " - Sample-to-Sample Distances.png"),
                               collapse = ""),
          width = 1200, height = 1200, units = "px", pointsize = 10, res = 200, 
          bg = "white", family = "", symbolfamily="default")
      pheatmap(sampleDistMatrix,
               clustering_distance_rows=sampleDists,
               clustering_distance_cols=sampleDists,
               col=colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
               main = stri_join(c(plot_title, "(Sample-to-Sample Distances)"), collapse = "\n"))
      dev.off()
    } else {
      pheatmap(sampleDistMatrix,
               clustering_distance_rows=sampleDists,
               clustering_distance_cols=sampleDists,
               col=colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
               main = stri_join(c(plot_title, "(Sample-to-Sample Distances)"), collapse = "\n"))
    }
}

QC_heatmaps_batch_corrected <- function(dds, n, filename_start, plot_title, remove_batch, nsplit){
  # Transform data -------------------------------------------------------------
    vsd <- vst(dds)
    mat <- assay(vsd)
    mm <- model.matrix(~Treatment, colData(vsd))
    
    if(remove_batch == "TRUE") {
      mat <- removeBatchEffect(mat, batch=vsd$Cohort, design=mm)
      assay(vsd) <- mat
    }

  # Heatmap of count matrix ----------------------------------------------------
    select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:n]
    # df <- as.data.frame(colData(dds)[, c("Treatment", "Cohort")])
    df <- as.data.frame(colData(dds)[, c("Treatment")])
    rownames(df) <- rownames(as.data.frame(colData(dds)))
    colnames(df) <- c("Treatment")
    
    if(!is.null(filename_start)) {
      png(filename = stri_join(c("QC_results/Heatmaps/", filename_start,".png"),
                               collapse = ""),
          width = 2000, height = 2000, units = "px", pointsize = 8, res = 250,
          bg = "white", family = "", symbolfamily="default")
      pheatmap(assay(vsd)[select,],
               # pheatmap(counts(dds,normalized=TRUE)[select,],
               cluster_rows=TRUE, show_rownames=FALSE, 
               # color = pnw_palette("Moth", 20, type = "continuous"),
               color = wes_palette("Zissou1", n = 100, type = "continuous"),
               # annotation_colors = list(Treatment=c(PL23="#1d457f", PL23_2DG="#d8aedd", R848="#cc5c76", R848_2DG = "#ffc3a3")),
               cutree_cols=nsplit, cluster_cols=TRUE, annotation_col=df, show_colnames = FALSE,
               # labels_col = colData(dds)$Label_Name,
               main = stri_join(c(plot_title, "(Variance Stabilizing Transformation)"), collapse = "\n"))
      dev.off()
    } else {
      pheatmap(assay(vsd)[select,],
      # pheatmap(counts(dds,normalized=TRUE)[select,],
               cluster_rows=TRUE, show_rownames=FALSE, 
               # color = pnw_palette("Moth", 20, type = "continuous"),
               color = wes_palette("Zissou1", n = 100, type = "continuous"),
               # annotation_colors = list(Treatment=c(PL23="#1d457f", PL23_2DG="#d8aedd", R848="#cc5c76", R848_2DG = "#ffc3a3")),
               cutree_cols=nsplit, cluster_cols=TRUE, annotation_col=df, show_colnames = FALSE,
               # labels_col = colData(dds)$Label_Name,
               main = stri_join(c(plot_title, "(Variance Stabilizing Transformation)"), collapse = "\n"))
    }
}


QC_PCAplot <- function(dds, filename_start, plot_title, batch_effect){
  vsd <- vst(dds)
  
  if(is.null(batch_effect)){
    pcaData <- plotPCA(vsd, intgroup=c("Treatment"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    if(!is.null(filename_start)){
      png(filename = stri_join(c("QC_results/PCA_plots/", filename_start, ".png"),
                               collapse = ""),
          width = 1500, height = 1500, units = "px", pointsize = 10, res = 200,
          bg = "white", family = "", symbolfamily="default")
    }
    ggplot(pcaData, aes(PC1, PC2, color=Treatment)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() +
      labs(title = plot_title)
  } else if(batch_effect == FALSE){
    pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Cohort"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    if(!is.null(filename_start)){
      png(filename = stri_join(c("QC_results/PCA_plots/", filename_start, ".png"),
                               collapse = ""),
          width = 1500, height = 1500, units = "px", pointsize = 10, res = 300,
          bg = "white", family = "", symbolfamily="default")
    }
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

    if(!is.null(filename_start)){
      png(filename = stri_join(c("QC_results/PCA_plots/", filename_start, ".png"),
                             collapse = ""),
        width = 1500, height = 1500, units = "px", pointsize = 10, res = 200,
        bg = "white", family = "", symbolfamily="default")
    }
    ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Cohort)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() +
      labs(title = plot_title)
  } else print("batch_effect must be either TRUE or FALSE")
}

results_wrapper <- function(dds, cons, IDs){
  res <- results(dds, contrast = cons)
  res <- data.frame(subset(res, !is.na(padj)))
  # Note: There are a few reasons why a p value or padj value would be NA
  # According to the DESeq2 manual, these are the reasons:
  # If within a row, all samples have zero counts, the baseMean column will be zero, and the LFC estimates, p value and padj will all be NA.
  # If a row contains a sample with an extreme count outlier then the p value and padj will be set to NA.
  # If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only padj will be set to NA.
  tempcols <- colnames(res)
  res <- data.frame(rownames(res), res)
  colnames(res) <- c("ensembl_gene_id", tempcols)
  return(right_join(IDs, res))
}

summary_wrapper <- function(res, comparison, experiment, p_cutoff, lfc_cutoff){
  if(missing(p_cutoff)) p_cutoff <- 0.05
  if(missing(lfc_cutoff)) lfc_cutoff <- 1

  degs_up <- nrow(subset(res, res$padj <= p_cutoff & res$log2FoldChange >= lfc_cutoff))
  degs_down <- nrow(subset(res, res$padj <= p_cutoff & res$log2FoldChange <= (-1*lfc_cutoff)))
  
  DEG_summary <- data.frame(experiment, comparison, lfc_cutoff, p_cutoff, degs_up, degs_down)
  colnames(DEG_summary) <- c("Experiment", "Comparison", "LFC_cutoff", "P_adj_cutoff", "Up", "Down")
  DEG_summary
}


sig_DEG_table <- function(res, ID_type, lfc_cutoff, padj_cutoff, file_start, file_format){
  if(missing(ID_type)) ID_type <- "all"
  if(missing(lfc_cutoff)) lfc_cutoff <- 1
  if(missing(padj_cutoff)) padj_cutoff <- 0.05
  if(missing(file_start)) file_start <- NULL
  if(missing(file_format)) file_format <- "xlsx"
  
  if(ID_type == "ensembl_gene_id") res <- res[, c(1, 6, 10)]
  else if(ID_type == "external_gene_name") res <- res[, c(2, 6, 10)]
  else if(ID_type == "entrezgene") res <- res[, c(4, 6, 10)]
  else if (ID_type == "all") res <- res[, c(1, 2, 4, 6, 10)]
  else print("specify ID type: ensembl_gene_id, external_gene_name, entrezgene, or all")
  
  up <- subset(res, res$padj <= padj_cutoff & res$log2FoldChange >= lfc_cutoff)
  down <- subset(res, res$padj <= padj_cutoff & res$log2FoldChange <= (-1*lfc_cutoff))
  all <- subset(res, res$padj <= padj_cutoff & abs(res$log2FoldChange) >= lfc_cutoff)
  
  if(!is.null(file_start)){
    if(file_format == "xlsx"){
      wb <- createWorkbook(stri_join(c("Output/Significant DEGs/", file_start, ".xlsx"), collapse = ""))
      addWorksheet(wb, "Upregulated DEGs")
      addWorksheet(wb, "Downregulated DEGs")
      addWorksheet(wb, "All DEGs")
      writeData(wb, "Upregulated DEGs", up, rowNames = FALSE)
      writeData(wb, "Downregulated DEGs", down, rowNames = FALSE)
      writeData(wb, "All DEGs", all, rowNames = FALSE)
      saveWorkbook(wb, stri_join(c("Output/Significant DEGs/", file_start, ".xlsx"), collapse = ""), overwrite = TRUE)
      rm(wb)
    } else if (file_format == "csv"){
      write.table(up,
                  file = stri_join(c("Output/", file_start, "_Up.csv"), collapse = ""),
                  sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
      write.table(down,
                  file = stri_join(c("Output/", file_start, "_Down.csv"), collapse = ""),
                  sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
      write.table(all,
                  file = stri_join(c("Output/", file_start, "_All.csv"), collapse = ""),
                  sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
    } else print("specify file type: xlsx or csv")
  }
  return(list(up = up, down = down, all = all))
}
