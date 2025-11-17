# Import Rosalind data ---------------------------------------------------------
## This function reads in a .txt file containing raw count data exported from ROSALIND.
## You have to use the "_rawCountsWithAnnotations.txt" file for the columns to be correct.
## It then filters the count data to exclude Ig gene fragments (using RegEx format),
## then removes duplicated lines, and sets the rownames to be the ensembl gene ID.
## There are a few lines (that are currently commented out) that remove ribosomal RNA sequences.
  import_Rosalind_data <- function(res_file, ipa_format){
    if(missing(ipa_format)) ipa_format <- FALSE
    
    temp_cts <- read.delim(res_file)
    temp_cts <- temp_cts %>%
      filter(!grepl("Ighd\\d", external_gene_name)) %>%
      filter(!grepl("Ighj\\d", external_gene_name)) %>%
      filter(!grepl("Ighv\\d", external_gene_name)) %>%
      filter(!grepl("Igkj\\d", external_gene_name)) %>%
      filter(!grepl("Igkv\\d", external_gene_name)) %>%
      filter(!grepl("Iglj\\d", external_gene_name)) %>%
      filter(!grepl("Iglv\\d", external_gene_name))
      # filter(!grepl("7SK", external_gene_name)) %>%
      # filter(!grepl("5_8S_rRNA", external_gene_name)) %>%
      # filter(!grepl("5S_rRNA", external_gene_name))
    
    if(ipa_format == TRUE) temp_cts <- temp_cts[, c(1:2, 4:5, 7, 9, 12:ncol(temp_cts))]
    else temp_cts <- temp_cts[, c(1:4, 12:ncol(temp_cts))]

    temp_cts <- temp_cts[!duplicated(temp_cts), ]
    temp_cts <- temp_cts[!duplicated(temp_cts$ensembl_gene_id), ]
    rownames(temp_cts) <- temp_cts$ensembl_gene_id
    return(temp_cts)
  }

# Generate heatmaps for quality control ----------------------------------------
## Function inputs:
### dds = a DESeqDataSet object
### filename_start = the beginning of whatever you want the resulting images to be saved as
### plot_title = the title of the plots
### plot_cohorts = a TRUE or FALSE value indicating whether to add cohorts to the column data. Default is TRUE
### max_genes = maximum number of genes to include in the heatmap (this corresponds to the number of rows)
## if you use NULL for filename_start, the image won't save but it'll appear in RStudio (this is good if you're running it for the first time and want to try things out)
  QC_heatmaps <- function(dds, filename_start, plot_title, plot_cohorts, max_genes){
    if(missing(plot_cohorts)) plot_cohorts <- TRUE
    if(missing(max_genes)) max_genes <- 20
    
    # Transform data -----------------------------------------------------------
      rld <- rlog(dds)

    # Heatmap of count matrix --------------------------------------------------
      select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:max_genes]
      
      if(plot_cohorts == TRUE) df <- as.data.frame(colData(dds)[, c("Treatment", "Cohort")])
      else {
        df <- as.data.frame(colData(dds)[, c("Treatment")])
        colnames(df) <- c("Treatment")
        rownames(df) <- colData(dds)[, c("Sample")]
      }

      if(!is.null(filename_start)) {
        # if filename_start is not NULL, then save the image to the "QC_results/Heatmaps/" folder
        
        if(!file.exists("QC_results/")) dir.create("QC_results/")
        if(!file.exists("QC_results/Heatmaps/")) dir.create("QC_results/Heatmaps/")
          # if there is not already a QC_results/Heatmaps/" folder, create one

        png(filename = stri_join(c("QC_results/Heatmaps/", filename_start, ".png"), collapse = ""),
            width = 1200, height = 1200, units = "px", pointsize = 10, res = 200,
            bg = "white", family = "", symbolfamily="default")
        pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
                 cluster_cols=TRUE, annotation_col=df, labels_col = colData(dds)$Label_Name,
                 main = plot_title)
        dev.off()
      } else {
        pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
                 cluster_cols=TRUE, annotation_col=df, labels_col = colData(dds)$Label_Name,
                 main = plot_title)
      }
  }
  
# Generate sample-to-sample distance plots -------------------------------------
## Function inputs:
### dds = a DESeqDataSet object
### filename_start = the beginning of whatever you want the resulting images to be saved as
### plot_title = the title of the plots
## if you use NULL for filename_start, the image won't save but it'll appear in RStudio (this is good if you're running it for the first time and want to try things out)
  sample_to_sample_plot <- function(dds, filename_start, plot_title){
    vsd <- vst(dds)
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$Label_Name)
    colnames(sampleDistMatrix) <- NULL
    
    if(!is.null(filename_start)) {
      if(!file.exists("QC_results/")) dir.create("QC_results/")
      if(!file.exists("QC_results/Sample-to-Sample_Distances/")) dir.create("QC_results/Sample-to-Sample_Distances/")
        
      png(filename = stri_join(c("QC_results/Sample-to-Sample_Distances/", filename_start,
                                 " - Sample Distances.png"),
                               collapse = ""),
          width = 1200, height = 1200, units = "px", pointsize = 10, res = 200, 
          bg = "white", family = "", symbolfamily="default")
      pheatmap(sampleDistMatrix,
               clustering_distance_rows=sampleDists,
               clustering_distance_cols=sampleDists,
               col=colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
               main = stri_join(c("Sample-to-Sample Distances:", plot_title), collapse = "\n"))
      dev.off()
    } else {
      pheatmap(sampleDistMatrix,
               clustering_distance_rows=sampleDists,
               clustering_distance_cols=sampleDists,
               col=colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
               main = stri_join(c("Sample-to-Sample Distances:", plot_title), collapse = "\n"))
    }
  }
  

# Generate PCA plots for quality control ---------------------------------------
## Function inputs:
### dds = a DESeqDataSet object
### filename_start = the beginning of whatever you want the resulting images to be saved as
### plot_title = the title of the plots
### batch_effect = a TRUE, FALSE, or NULL value indicating whether to run batch effect correction
#### (use batch_effects = NULL if the experiment does not contain batches)
### draw_ellipse = a TRUE or FALSE value indicating whether to draw ellipses around the groups
#### (if missing, this value will be set to FALSE)
## if you use NULL for filename_start, the image won't save but it'll appear in RStudio (this is good if you're running it for the first time and want to try things out)
  QC_PCAplot <- function(dds, filename_start, plot_title, batch_effect, draw_ellipse){
    if(missing(batch_effect)) batch_effect <- NULL
    if(missing(draw_ellipse)) draw_ellipse <- FALSE
    vsd <- vst(dds)
    
    if(!is.null(filename_start)){
      if(!file.exists("QC_results/")) dir.create("QC_results/")
      if(!file.exists("QC_results/PCA_plots/")) dir.create("QC_results/PCA_plots/")
      # if there is not already a QC_results/PCA_plots/" folder, create one

      png(filename = stri_join(c("./QC_results/PCA_plots/", filename_start, ".png"), collapse = ""),
          width = 1800, height = 1800, units = "px", pointsize = 10, res = 350,
          bg = "white", family = "", symbolfamily="default")
    }
    
    if(is.null(batch_effect)){
      pcaData <- plotPCA(vsd, intgroup=c("Treatment"), returnData=TRUE)
      percentVar <- round(100 * attr(pcaData, "percentVar"))
      
      ggplot(pcaData, aes(PC1, PC2, color=Treatment)) +
        geom_point(size=3) +
        geom_text(aes(label = dds$Mouse_Number), vjust = -0.5) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        coord_fixed() +
        labs(title = plot_title) +
        if(draw_ellipse == TRUE) stat_ellipse(aes(group = Treatment)) else {}
    } else {
      if(batch_effect == TRUE){
        mat <- assay(vsd)
        mm <- model.matrix(~Treatment, colData(vsd))
        mat <- removeBatchEffect(mat, batch=vsd$Cohort, design=mm)
        assay(vsd) <- mat
      } else if(batch_effect != FALSE){
        print("batch_effect must be either TRUE or FALSE")
        return(NULL)
      }
      
      pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Cohort"), returnData=TRUE)
      percentVar <- round(100 * attr(pcaData, "percentVar"))
      
      ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Cohort)) +
        geom_point(size=3) +
        geom_text(aes(label = dds$Mouse_Number), vjust = -0.5) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        coord_fixed() +
        theme(legend.position = 'bottom') +
        labs(title = plot_title) +
        labs(title = plot_title) +
          if(draw_ellipse == TRUE) stat_ellipse(aes(group = Treatment)) else {}
    } 
  }

# Wrapper for the DESeq2 results() function ------------------------------------
## Function inputs:
### dds = a DESeqDataSet object
### cons = which contrasts you want to perform
### IDs = a table containing additional data to add to the results, such as alternative gene IDs. Must contain a "ensembl_gene_id" column.
#### (if IDs is missing or NULL, this step will be skipped)
results_wrapper <- function(dds, cons, IDs){
  if(missing(IDs)) IDs <- NULL
  
  res <- results(dds, contrast = cons)
  res <- data.frame(subset(res, !is.na(padj)))
  #  Note: There are a few reasons why a p value or padj value would be NA
  # # According to the DESeq2 manual, these are the reasons:
  # # # If within a row, all samples have zero counts, the baseMean column will be zero, and the LFC estimates, p value and padj will all be NA.
  # # # If a row contains a sample with an extreme count outlier then the p value and padj will be set to NA.
  # # # If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only padj will be set to NA.
  tempcols <- colnames(res)
  res <- data.frame(rownames(res), res)
  colnames(res) <- c("ensembl_gene_id", tempcols)
  
  if(is.null(IDs)) return(res)
  else return(right_join(IDs, res))
}

# Make a heatmap of DEGs -------------------------------------------------------
DEG_heatmap <- function(dds, dds_results, factors, color_list, lfc_cutoff, padj_cutoff, count_type, dds_cols, plot_title, filename, h, w){
  if(missing(color_list)) {
    color_list <- list(Treatment = c(PL23 = "#0f85a0", PL23_2DG = "#dd4124",
                                      R848 = "#0f85a0", R848_2DG = "#dd4124",
                                     NP = "#0f85a0", NP_2DG = "#dd4124",
                                     Control = "#0f85a0", "2DG" = "#dd4124"),
                        Cohort = c(A1 = "#00496f", A2 = "#edd746", B1 = "#00496f",
                                   M1 = "#00496f", M2 = "#edd746"))
    
    # color_list <- list(Treatment = c(PL23 = "#00496f", PL23_2DG = "#0f85a0", 
    #                                  R848 = "#ed8b00", R848_2DG = "#dd4124",
    #                                  NP = "#0f85a0", NP_2DG = "#dd4124",
    #                                  Control = "#0f85a0", "2DG" = "#dd4124"),
    #                    Cohort = c(A1 = "#00496f", A2 = "#edd746", B1 = "#00496f", 
    #                               M1 = "#00496f", M2 = "#edd746"))
  }
  if(missing(lfc_cutoff)) lfc_cutoff <- 1
  if(missing(padj_cutoff)) padj_cutoff <- 0.05
  if(missing(count_type)) count_type <- "norm_counts"
  if(missing(dds_cols)) dds_cols <- NULL
  if(missing(plot_title)) plot_title <- "Plot Title"
  if(missing(filename)) filename <- NULL
  if(missing(h)) h <- 2000
  if(missing(w)) w <- 1500

  # make a list of DEGs from the DEG results -----------------------------------
    DEGs <- c()
    for (df in list(dds_results)) {
      res <- dplyr::filter(df, padj < padj_cutoff & abs(log2FoldChange) > abs(lfc_cutoff))
      res <- res$ensembl_gene_id
      DEGs <- append(DEGs, res)
    }
    DEGs <- DEGs[!duplicated(DEGs)]
  
  # Make a vector of count matrix row numbers in decreasing order of row means -
    select <- counts(dds, normalized=TRUE)
    if(!is.null(dds_cols)) select <- select[, dds_cols]
    select <- select[rownames(select) %in% DEGs, ]
    select <- order(rowMeans(select), decreasing=TRUE)
    head(select)
  
  # Make a dataframe with sample info ------------------------------------------
    # If dds_col is NULL, use all sample info from the dds object
    # Else, only get sample info from rows corresponding to the count columns of interest
    if(is.null(dds_cols)) {
      df <- as.data.frame(colData(dds)[, factors])
      rownames(df) <- rownames(as.data.frame(colData(dds)))
    } else {
      df <- as.data.frame(colData(dds)[dds_cols, factors])
      rownames(df) <- rownames(as.data.frame(colData(dds)))[dds_cols]
    }
    colnames(df) <- factors

  # Sort the normalized counts by decreasing mean counts -----------------------
    if(count_type == "norm_counts") {
      if(!is.null(dds_cols)) {
        htmp_cts <- counts(dds, normalized = TRUE)[select, dds_cols]
      } else htmp_cts <- counts(dds, normalized = TRUE)[select, ]
    } else if(count_type == "rlog") {
      rld <- rlog(dds)
      if(!is.null(dds_cols)) {
        htmp_cts <- assay(rld)[select, dds_cols]
      } else htmp_cts <- assay(rld)[select, ]
    } else {
      print("Error: count_type must be either norm_counts or rlog")
      return(NULL)
    }
    
  # Plot heatmap with unsupervised column clustering ---------------------------
    if(!is.null(filename)) {
      if(!file.exists("Output/")) dir.create("Output/")
      if(!file.exists("Output/DEG_heatmaps/")) dir.create("Output/DEG_heatmaps/")
      # if there is not already a QC_results/DEG_heatmaps/" folder, create one
      
      png(filename = stri_join(c("Output/DEG_heatmaps/Unsupervised - ", filename,".png"), collapse = ""),
          width = w, height = h, units = "px", pointsize = 8, res = 250,
          bg = "white", family = "", symbolfamily="default")
    }
    
    htmp <- ComplexHeatmap::pheatmap(htmp_cts,
                                     cluster_rows = TRUE, show_rownames = FALSE,
                                     cluster_cols = TRUE, show_colnames = FALSE,
                                     annotation_col = df, scale = "row",
                                     annotation_colors = color_list,
                                     heatmap_legend_param = list(title = "Z-score"),
                                     main = plot_title)
    draw(htmp, legend_grouping = "original", merge_legends = TRUE)
    
    if(!is.null(filename))  dev.off()
    
  # Plot heatmap with supervised column clustering -----------------------------
    if(!is.null(filename)) {
      png(filename = stri_join(c("Output/DEG_heatmaps/Supervised - ", filename,".png"), collapse = ""),
          width = w, height = h, units = "px", pointsize = 8, res = 250,
          bg = "white", family = "", symbolfamily="default")
    }
  
    htmp <- ComplexHeatmap::pheatmap(htmp_cts,
                                     cluster_rows = TRUE, show_rownames = FALSE,
                                     cluster_cols = FALSE, show_colnames = FALSE,
                                     annotation_col = df, scale = "row",
                                     annotation_colors = color_list,
                                     heatmap_legend_param = list(title = "Z-score"),
                                     main = plot_title)
    draw(htmp, legend_grouping = "original", merge_legends = TRUE)
    if(!is.null(filename)) dev.off()
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
