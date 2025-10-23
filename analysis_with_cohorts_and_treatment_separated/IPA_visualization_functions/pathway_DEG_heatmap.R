setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_and_treatment_separated/Output/chen_license_pval0.05_lfc0.6")
getwd()

get_pathway_molecules <- function(filename, pwname){
  res <- read.delim(stri_join("IPA_output_individual_files/", filename, sep = ""), 
                    row.names = NULL, col.names = c("Pathway", "pval", "Zscore", "Ratio", "Molecules"))
  res <- subset(res, Pathway == pwname)
  if(nrow(res) > 0) res <- stri_split_regex(res$Molecules, ",")[[1]]
  else return(NULL)
  return(res)
}

get_full_mols <- function(pwname){
  pl23_pw_res <- get_pathway_molecules("PL2-3_Canonical_Pathways.txt", pwname)
  mrllpr_pw_res <- get_pathway_molecules("AM14_MRLlpr_Canonical_Pathways.txt", pwname)
  r848_pw_res <- get_pathway_molecules("R848_Canonical_Pathways.txt", pwname)
  np_pw_res <- get_pathway_molecules("NP_Canonical_Pathways.txt", pwname)
  
  mols <- c()
  if(!is.null(pl23_pw_res)) mols <- c(mols, pl23_pw_res)
  if(!is.null(mrllpr_pw_res)) mols <- c(mols, mrllpr_pw_res)
  if(!is.null(r848_pw_res)) mols <- c(mols, r848_pw_res)
  if(!is.null(np_pw_res)) mols <- c(mols, np_pw_res)
  mols <- mols[!duplicated(mols)]
  return(mols)
}

pw_expression_single_dataset <- function(sheetname, mols, comparison_name){
  if(missing(comparison_name)) comparison_name <- NULL
  
  expr_res <- read.xlsx("../DESeq2_results.xlsx", sheet = sheetname)
  expr_res <- filter(expr_res, hgnc_symbol %in% mols)
  expr_res <- expr_res[, c(2, 5, 8, 11, 12)]
  
  expr_res <- data.frame(expr_res$external_gene_name, 
                         c(rep(comparison_name, nrow(expr_res))), 
                         -1*log10(expr_res$pvalue), 
                         expr_res$log2FoldChange)
  colnames(expr_res) <- c("ID", "comparison", "pvalue", "LFC")
  
  return(expr_res)
}

get_pathway_expression <- function(pwname){
  # pl23_pw_res <- get_pathway_molecules("PL2-3_Canonical_Pathways.txt", pwname)
  # mrllpr_pw_res <- get_pathway_molecules("AM14_MRLlpr_Canonical_Pathways.txt", pwname)
  # r848_pw_res <- get_pathway_molecules("R848_Canonical_Pathways.txt", pwname)
  # np_pw_res <- get_pathway_molecules("NP_Canonical_Pathways.txt", pwname)
  # 
  # mols <- c()
  # if(!is.null(pl23_pw_res)) mols <- c(mols, pl23_pw_res)
  # if(!is.null(mrllpr_pw_res)) mols <- c(mols, mrllpr_pw_res)
  # if(!is.null(r848_pw_res)) mols <- c(mols, r848_pw_res)
  # if(!is.null(np_pw_res)) mols <- c(mols, np_pw_res)
  # mols <- mols[!duplicated(mols)]
  
  mols <- get_full_mols(pwname)
  
  pl23_expr <- pw_expression_single_dataset("AM14 Co1 - PL23_2DG_vs_Ctrl", mols, "PL2-3")
  mrllpr_expr <- pw_expression_single_dataset("AM14 MRLlpr Co2 - 2DG_vs_Ctrl", mols, "AM14 MRLlpr")
  r848_expr <- pw_expression_single_dataset("AM14 Co2 - R848_2DG_vs_Ctrl", mols, "R848")
  np_expr <- pw_expression_single_dataset("B1-8 - 2DG_vs_Ctrl", mols, "NP")
  
  res <- full_join(pl23_expr, mrllpr_expr)
  res <- full_join(res, r848_expr)
  res <- full_join(res, np_expr)
  
  return(res)
}

pathway_molecules_dotplot <- function(pwname){
  expr_res <- get_pathway_expression(pwname)
  
  ggplot(data = expr_res, 
         aes(x = factor(comparison), y = factor(ID))) + 
    geom_point(aes(size = pvalue, color = LFC)) + 
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12)) +
    scale_colour_gradientn(colors = wes_palette("Zissou1")) +
    ylab("") + 
    xlab("") + 
    scale_y_discrete(position = "right") +
    labs(title = pwname) 
}

complement_system_mols <- get_full_mols("Complement System")
complement_system_mols

complement_system <- get_pathway_expression("Complement System")
complement_system

pathway_molecules_dotplot("Complement System")


# pathway_DEG_zScores <- function(dds, mols){
#   zscores <- counts(dds, normalized=TRUE)
#   zscores <- zscores[rownames(zscores) %in% mols, ]
#   res_colnames <- colnames(zscores)
#   zscores <- t(zscores)
#   zscores <- scale(zscores)
#   zscores <- t(zscores)
#   zscores <- data.frame(rownames(zscores), zscores)
#   colnames(zscores) <- c("ensembl_gene_id", res_colnames)
#   rm(res_colnames)
#   return(zscores)
# }
# 
# pl23_zscores <- pathway_DEG_zScores(dds_AM14trans_PL23, complement_system_mols)
# pl23_zscores









pl23_DEG_df <- deseq_res$AM14transfer_PL23[, c(1, 2, 8, 11, 12)]
pl23_DEG_df <- filter(pl23_DEG_df, abs(log2FoldChange) > 0.6 & padj <= 0.05)
nrow(pl23_DEG_df)
head(pl23_DEG_df)

MRLlpr_DEG_df <- deseq_res$AM14MRLlpr[, c(1, 2, 8, 11, 12)]
MRLlpr_DEG_df <- filter(MRLlpr_DEG_df, abs(log2FoldChange) > 0.6 & padj <= 0.05)
nrow(MRLlpr_DEG_df)
head(MRLlpr_DEG_df)

merged_DEG_df <- pl23_DEG_df
colnames(merged_DEG_df)  <- c("ensembl_gene_id", "external_gene_name", "PL23_LFC", "PL23_pvalue", "PL23_padj")
merged_DEG_df <- inner_join(merged_DEG_df, MRLlpr_DEG_df)
colnames(merged_DEG_df)  <- c("ensembl_gene_id", "external_gene_name", "PL23_LFC", "PL23_pvalue", "PL23_padj", "MRLlpr_LFC", "MRLlpr_pvalue", "MRLlpr_padj")
head(merged_DEG_df)
nrow(merged_DEG_df)

ID_conversion <- merged_DEG_df[, 1:2]
ID_conversion

DEGs <- merged_DEG_df$ensembl_gene_id
head(DEGs)


zscore_matrix <- function(dds){
  zscores <- counts(dds, normalized=TRUE)
  zscores <- zscores[rownames(zscores) %in% DEGs, ]
  res_colnames <- colnames(zscores)
  zscores <- t(zscores)
  zscores <- scale(zscores)
  zscores <- t(zscores)
  zscores <- data.frame(rownames(zscores), zscores)
  colnames(zscores) <- c("ensembl_gene_id", res_colnames)
  rm(res_colnames)
  return(zscores)
}

pl23_zscores <- zscore_matrix(dds_AM14trans_PL23)
head(pl23_zscores)

MRLlpr_zscores <- zscore_matrix(dds_AM14MRLlpr)
head(MRLlpr_zscores)

DEG_zscores <- inner_join(pl23_zscores, MRLlpr_zscores)
DEG_zscores <- inner_join(ID_conversion, DEG_zscores)
rownames(DEG_zscores) <- DEG_zscores$external_gene_name
DEG_zscores <- as.matrix(DEG_zscores[, c(3:ncol(DEG_zscores))])
head(DEG_zscores)
DEG_zscores

coldata_df <- data.frame(Treatment = c(rep("Control", 3), rep("2DG", 3), rep("Control", 6), rep("2DG", 7)),
                         Experiment = c(rep("AM14 + PL2-3", 6), rep("AM14 MRLlpr", 13)))
coldata_df

png(filename = ("Output/DEG_heatmaps/DEGs common in AM14 PL2-3 and AM14 MRLlpr.png"),
    width = 2000, height = 3500, units = "px", pointsize = 8, res = 400,
    bg = "white", family = "", symbolfamily="default")

htmp <- ComplexHeatmap::pheatmap(DEG_zscores,
                                 cluster_rows = TRUE, show_rownames = TRUE,
                                 cluster_cols = TRUE, show_colnames = FALSE, 
                                 annotation_col = coldata_df, 
                                 scale = "none",
                                 cutree_rows = 3, 
                                 annotation_colors = list(Experiment = c("AM14 + PL2-3" = "#0f85a0", "AM14 MRLlpr" = "#00496f"),
                                                          Treatment = c(Control = "#edd746", "2DG" = "#dd4124")),
                                 # main = "DEGs Common to AM14 Adoptive Transfer PL2-3 +/- 2DG and AM14 MRL/lpr +/- 2DG",
                                 heatmap_legend_param = list(title = "Z-score", direction = "horizontal"))

draw(htmp, legend_grouping = "original", merge_legends = TRUE, heatmap_legend_side = "top")
dev.off()

