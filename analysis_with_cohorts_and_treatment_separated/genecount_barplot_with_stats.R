find_ensembl_ID <- function(search_name) {
  dat <- dplyr::filter(gene_IDs_full, external_gene_name == search_name)
  dat$ensembl_gene_id
}

normcount_table <- function(gene_name, dds){
  ID <- find_ensembl_ID(gene_name)
  
  res <- counts(dds, normalized = TRUE)
  res <- res[rownames(res) %in% ID, ]
  res <- data.frame(Sample = names(res), counts = res)
  res <- left_join(res, coldata, by = "Sample")
  res <- res[ , c(1:2, 6, 8, 9)] %>%
    group_by(Treatment, Drug, heatmap_col1)
  
  colnames(res) <- c("Sample", "counts", "Treatment", "Drug", "Experiment")
  
  return(res)
  
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

count_boxplot <- function(gene_name, counts_table){
  # normcounts_res <- normcount_table(gene_name, dds)
  
  ggboxplot(counts_table, x = "Drug", y = "counts", 
            add = c("mean", "jitter"), # Add mean and individual points
            ylab = "counts", xlab = "Treatment") +
    # geom_boxplot() +
    # geom_jitter(width = 0.2, alpha = 0.5) +
    theme_bw() +
    labs(title = gene_name) + 
    geom_text(aes(label = Sample), nudge_y = 50, size = 3) +
    # geom_text(aes(label = ifelse(is_outlier(counts), Sample, NA)), nudge_y = 50, size = 3) +
    facet_wrap( ~ Experiment, scales = 'fixed', ncol = 4)
}

multi_boxplot <- function(gene_name){
  pl23_normcounts <- normcount_table(gene_name, dds_AM14trans_PL23)
  mrl_normcounts <- normcount_table(gene_name, dds_AM14MRLlpr)
  r848_normcounts <- normcount_table(gene_name, dds_AM14trans_R848)
  np_normcounts <- normcount_table(gene_name, dds_B18trans)
  
  full_table <- full_join(pl23_normcounts, mrl_normcounts)
  if(nrow(r848_normcounts) > 0) full_table <- full_join(full_table, r848_normcounts)
  if(nrow(np_normcounts) > 0) full_table <- full_join(full_table, np_normcounts)
  
  count_boxplot(gene_name, full_table)
  
  # count_boxplot(gene_name, pl23_normcounts) +
  #   count_boxplot(gene_name, mrl_normcounts) +
  #   count_boxplot(gene_name, r848_normcounts) +
  #   count_boxplot(gene_name, np_normcounts)
  
  # return(mrl_normcounts)
}

# cd3_mrl <- normcount_table("Cd3d", dds_AM14MRLlpr)
# count_boxplot("Cd3d", cd3_mrl)

# test <- normcount_table("Tbx21", dds_AM14MRLlpr)
# test

multi_boxplot("Tbx21")

multi_boxplot("Cd3d")
multi_boxplot("Cd8a")


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
