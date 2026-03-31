setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_and_treatment_separated")
getwd()

find_DEG_set_overlap <- function(deg1, colname1, deg2, colname2, p_value_cutoff, lfc_cutoff){
  if(missing(p_value_cutoff)) p_value_cutoff <- 0.05
  if(missing(lfc_cutoff)) lfc_cutoff <- 1
  
  filt1 <- deg1[deg1$padj < p_value_cutoff, ]
  # filt1 <- filt1[abs(filt1$log2FoldChange) > lfc_cutoff, ]
  filt1 <- filt1[ , 1:2]
  colnames(filt1) <- c("Gene", colname1)
  
  filt2 <- deg2[deg2$padj < p_value_cutoff, ]
  # filt2 <- filt2[abs(filt2$log2FoldChange) > lfc_cutoff, ]
  filt2 <- filt2[ , 1:2]
  colnames(filt2) <- c("Gene", colname2)
  
  merge_filt <- inner_join(filt1, filt2, by = "Gene")
  
  merge_filt$same_direction <- ifelse(
    test = merge_filt[, 2] < 0 & merge_filt[, 3] < 0, 
    yes = -1,
    no = ifelse(
      test = merge_filt[, 2] > 0 & merge_filt[, 3] > 0, 
      yes = 1, 
      no = 0
    )
  )
  
  return(merge_filt)
}




vp_data = list(
  PL23_2DG_v_Ctrl = deseq_res$AM14transfer_PL23[, c(2, 8, 12)],
  R848_2DG_v_Ctrl = deseq_res$AM14transfer_R848[, c(2, 8, 12)],
  AM14MRLlpr = deseq_res$AM14MRLlpr[, c(2, 8, 12)],
  B18 = deseq_res$B18transfer[, c(2, 8, 12)])
head(vp_data$AM14MRLlpr)


pl23_x_MRLlpr <- find_DEG_set_overlap(vp_data$PL23_2DG_v_Ctrl, "PL23",
                                      vp_data$AM14MRLlpr, "MRLlpr")
head(pl23_x_MRLlpr)

# ggplot(pl23_x_MRLlpr, aes(x = "AM14 PL2-3", y = "AM14 MRLlpr", fill = consistency)) +
ggplot(pl23_x_MRLlpr, aes(x = PL23, y = MRLlpr)) +
  geom_point() +
  # geom_point(aes(color = Direction_in_ABCs)) +
  # geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#444444") +
  # geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "#444444") +
  geom_hline(yintercept = c(0), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(0), linetype = "dashed", color = "gray") +
  geom_label_repel(
    data = pl23_x_MRLlpr[pl23_x_MRLlpr$same_direction != 0, ],
    aes(label = Gene),
    # size = 5,
    max.overlaps = Inf
  ) +
  # scale_fill_manual(
  #   values = c("direction matches" = "lightblue", "opposing directions" = "pink"),
  #   na.value = "white"
  # ) +
  # scale_color_manual(
  #   values = c("down" = "blue", "up" = "red", "neither" = "gray"),
  #   guide = guide_legend(title = "Direction in ABCs vs FoB\n(GSE99480)")
  # ) +
  labs(
    x = "LFC: PL2-3 2DG vs Ctrl",
    y = "LFC: AM14 MRL/lpr 2DG vs Ctrl",
    title = "Differentially Expressed Genes"
  ) +
  theme_bw()

