setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_and_treatment_separated/")
getwd()

find_ensembl_ID <- function(search_name) {
  dat <- dplyr::filter(gene_IDs_full, external_gene_name == search_name)
  dat$ensembl_gene_id
}

normcount_table <- function(gene_name, dds){
  ID <- find_ensembl_ID(gene_name)
  
  res <- counts(dds, normalized = TRUE)
  res <- res[rownames(res) %in% ID, ]
  res <- data.frame(Sample = names(res), counts = res)
  
  if(nrow(res) > 0){
    res <- left_join(res, coldata, by = "Sample")
    res <- res[ , c(1:2, 4, 6, 8, 9)] %>%
      group_by(Treatment, Drug, heatmap_col1)
    
    colnames(res) <- c("Sample", "counts", "Mouse_Number", "Treatment", "Drug", "Experiment")
    
    return(res)
  } else return(NULL)
  
  # summary_df <- res %>%
  #   summarise(
  #     mean_counts = mean(counts),
  #     sd_counts = sd(counts),
  #     se_counts = sd_counts / sqrt(n()), # Calculate standard error
  #     .groups = 'drop'
  #   )
  # 
  # anova_result <- aov(counts ~ Treatment, data = res)
  # tukey_result <- TukeyHSD(anova_result)
  # 
  # return(list(counts_table = res,
  #             summary = summary_df,
  #             ANOVA = anova_result,
  #             Tukey = tukey_result))
}

# is_outlier <- function(x){
#   return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
# }

count_boxplot <- function(gene_name, counts_table, alt_title){
  if(missing(alt_title)) alt_title <- NULL
  
  if(is.null(alt_title)) alt_title <- gene_name
  
  ggboxplot(counts_table, x = "Drug", y = "counts", 
            # outliers = FALSE,
            add = c("jitter"), label = "Mouse_Number", repel = TRUE,
            ggtheme = theme_bw(), ylab = "counts", xlab = "Treatment") +
    # geom_point(aes(group = Treatment), position = position_jitter(width = 0.3), size = 1.5) +
    # geom_text(aes(label = Mouse_Number), vjust = -0.5) +
    labs(title = alt_title) + 
    facet_wrap( ~ Experiment, scales = 'fixed', ncol = 4)
}

multi_boxplot <- function(gene_name, alt_title){
  if(missing(alt_title)) alt_title <- NULL
  
  pl23_normcounts <- normcount_table(gene_name, dds_AM14trans_PL23)
  mrl_normcounts <- normcount_table(gene_name, dds_AM14MRLlpr)
  r848_normcounts <- normcount_table(gene_name, dds_AM14trans_R848)
  np_normcounts <- normcount_table(gene_name, dds_B18trans)
  
  if(!is.null(pl23_normcounts)) {
    full_table <- pl23_normcounts
    if(!is.null(mrl_normcounts)) full_table <- full_join(full_table, mrl_normcounts)
    if(!is.null(r848_normcounts)) full_table <- full_join(full_table, r848_normcounts)
    if(!is.null(np_normcounts)) full_table <- full_join(full_table, np_normcounts)
  } else if(!is.null(mrl_normcounts)){
    full_table <- mrl_normcounts
    if(!is.null(r848_normcounts)) full_table <- full_join(full_table, r848_normcounts)
    if(!is.null(np_normcounts)) full_table <- full_join(full_table, np_normcounts)
  } else if(!is.null(r848_normcounts)){
    full_table <- r848_normcounts
    if(!is.null(np_normcounts)) full_table <- full_join(full_table, np_normcounts)
  } else if(!is.null(np_normcounts)) {
    full_table <- np_normcounts
  } else {
    print("gene not found")
    return(NULL)
  }
  
  count_boxplot(gene_name, full_table, alt_title)
}

# T and B cell markers ---------------------------------------------------------
  multi_boxplot("Cd3d")
  ggsave(filename = "QC_results/cell_marker_checks/Cd3d.png", width = 8, height = 4, units = "in", dpi = 300)
  
  multi_boxplot("Cd3e")
  ggsave(filename = "QC_results/cell_marker_checks/Cd3e.png", width = 4, height = 4, units = "in", dpi = 300)
  
  multi_boxplot("Cd3g")
  ggsave(filename = "QC_results/cell_marker_checks/Cd3g.png", width = 8, height = 4, units = "in", dpi = 300)
  
  multi_boxplot("Cd4")
  ggsave(filename = "QC_results/cell_marker_checks/Cd4.png", width = 4, height = 4, units = "in", dpi = 300)
  
  multi_boxplot("Cd8a")
  ggsave(filename = "QC_results/cell_marker_checks/Cd8a.png", width = 8, height = 4, units = "in", dpi = 300)
  
  multi_boxplot("Cd8b")
  # ggsave(filename = "QC_results/cell_marker_checks/Cd8b.png", width = 8, height = 4, units = "in", dpi = 300)
  
  multi_boxplot("Icos")
  ggsave(filename = "QC_results/cell_marker_checks/Icos.png", width = 4, height = 4, units = "in", dpi = 300)
  
  multi_boxplot("Cd40lg")
  ggsave(filename = "QC_results/cell_marker_checks/Cd40lg.png", width = 6, height = 4, units = "in", dpi = 300)
  
  multi_boxplot("Cd19")
  ggsave(filename = "QC_results/cell_marker_checks/Cd19.png", width = 8, height = 4, units = "in", dpi = 300)
  
  multi_boxplot("Ptprc", "Ptprc (B220)")
  ggsave(filename = "QC_results/cell_marker_checks/B220.png", width = 8, height = 4, units = "in", dpi = 300)
  

  
# Genes of interest ------------------------------------------------------------
  multi_boxplot("Cr2", "Cr2 (CD21)")
  ggsave(filename = "Output/DEG_barplots/Cr2.png", width = 9, height = 4.5, units = "in", dpi = 300)
  
  multi_boxplot("Tbx21", "Tbx21 (T-bet)")
  ggsave(filename = "Output/DEG_barplots/Tbx21.png", width = 9, height = 4.5, units = "in", dpi = 300)
  
  multi_boxplot("Itgax", "Itgax (CD11c)")
  ggsave(filename = "Output/DEG_barplots/Itgax.png", width = 9, height = 4.5, units = "in", dpi = 300)
  
  multi_boxplot("Itgam", "Itgax (CD11b)")
  ggsave(filename = "Output/DEG_barplots/Itgam.png", width = 9, height = 4.5, units = "in", dpi = 300)
  
  multi_boxplot("Cxcr3")
  ggsave(filename = "Output/DEG_barplots/Cxcr3.png", width = 9, height = 4.5, units = "in", dpi = 300)
  
  multi_boxplot("C1qc")
  ggsave(filename = "Output/DEG_barplots/C1qc.png", width = 9, height = 4.5, units = "in", dpi = 300)
  
  multi_boxplot("C1qb")
  ggsave(filename = "Output/DEG_barplots/C1qb.png", width = 9, height = 4.5, units = "in", dpi = 300)
  
  multi_boxplot("Ly6c1")
  ggsave(filename = "Output/DEG_barplots/Ly6c1.png", width = 9, height = 4.5, units = "in", dpi = 300)
  
  multi_boxplot("Batf3")
  ggsave(filename = "Output/DEG_barplots/Batf3.png", width = 9, height = 4.5, units = "in", dpi = 300)
  
  multi_boxplot("Ccl4")
  ggsave(filename = "Output/DEG_barplots/Ccl4.png", width = 9, height = 4.5, units = "in", dpi = 300)
  
  multi_boxplot("Lck")
  ggsave(filename = "Output/DEG_barplots/Lck.png", width = 9, height = 4.5, units = "in", dpi = 300)
  
  multi_boxplot("")
  ggsave(filename = "Output/DEG_barplots/.png", width = 9, height = 4.5, units = "in", dpi = 300)
  
  multi_boxplot("")
  ggsave(filename = "Output/DEG_barplots/.png", width = 9, height = 4.5, units = "in", dpi = 300)
  
  multi_boxplot("")
  ggsave(filename = "Output/DEG_barplots/.png", width = 9, height = 4.5, units = "in", dpi = 300)
  
  multi_boxplot("")
  ggsave(filename = "Output/DEG_barplots/.png", width = 9, height = 4.5, units = "in", dpi = 300)



# plot_counts <- function(counts_table, summary_df){
#   ggplot(counts_table, 
#          aes(x = Treatment, y = counts, fill = Treatment)) +
#     geom_boxplot() +
#     # geom_col(position = position_dodge(width = 0.9), color = "black") + # Dodge bars and add outline
#     # geom_errorbar(aes(ymin = mean_counts - se_counts, ymax = mean_counts + se_counts),
#     #               width = 0.25,
#     #               position = position_dodge(width = 0.9)) + # Add error bars
#     # labs(title = "Normalized Gene Counts by Treatment and Genotype",
#     #      x = "Treatment",
#     #      y = "Normalized Counts (Mean +/- SE)",
#     #      fill = "Genotype") +
#     # geom_jitter(aes(y = counts_table$counts), position = position_dodge(width = 0.9)) +
#     # # facet_wrap(~ gene, scales = "free_y") + # Create separate plots for each gene
#     # theme_bw() + # Use a clean theme
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis labels
# }


# plot_counts(mrl_tbx21$counts_table, mrl_tbx21$summary)

  # mrl_tbx21 <- normcount_table("Tbx21", dds_AM14MRLlpr)
  # mrl_tbx21$counts_table
  # mrl_tbx21$summary
  # mrl_tbx21$ANOVA
  # mrl_tbx21$Tukey
  # coldata
  # 
#   ggboxplot(mrl_tbx21$counts_table, x = "Treatment", y = "counts", 
#             add = c("mean", "jitter"), # Add mean and individual points
#             ylab = "counts", xlab = "Treatment") 
#     # stat_compare_means()
# 
#   
#   
#   summary_df <- normcount_table("Tbx21", dds_AM14MRLlpr) %>%
#     summarise(
#       mean_counts = mean(counts),
#       sd_counts = sd(counts),
#       se_counts = sd_counts / sqrt(n()), # Calculate standard error
#       .groups = 'drop'
#     )
# 
#     summary_df
#     
#     anova_result <- aov(counts ~ Treatment, data = normcount_table("Tbx21", dds_AM14MRLlpr))
#     tukey_result <- TukeyHSD(anova_result)
# 
#     anova_result    
# tukey_result    
