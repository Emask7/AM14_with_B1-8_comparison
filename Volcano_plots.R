# Format data for volcano plots ------------------------------------------------
  vpData_PL23 <- res_PL23_full[, c(2, 6, 10)]
  colnames(vpData_PL23) <- c("SYMBOL", "log2FoldChange", "padj")
  vpData_R848 <- res_R848_full[, c(2, 6, 10)]
  colnames(vpData_R848) <- c("SYMBOL", "log2FoldChange", "padj")
  vpData_PL23vR848 <- res_PL23vR848_full[, c(2, 6, 10)]
  colnames(vpData_PL23vR848) <- c("SYMBOL", "log2FoldChange", "padj")
  vpData_2DG_PL23vR848 <- res_2DG_PL23vR848_full[, c(2, 6, 10)]
  colnames(vpData_2DG_PL23vR848) <- c("SYMBOL", "log2FoldChange", "padj")
  
# Make plots -------------------------------------------------------------------
  EnhancedVolcano(vpData_PL23, 
                  lab = vpData_PL23$SYMBOL, 
                  pCutoff = 0.05, FCcutoff = 1,
                  x = 'log2FoldChange', y = 'padj', 
                  # xlim = c(-2.5, 1.5), 
                  ylim = c(0, 7),
                  title = "PL2-3 + 2DG vs PL2-3", 
                  subtitle = "(Adjusted p-values)", 
                  legendLabels = c("NS", expression(Log[2] ~ FC > 1), 
                                   expression(p - value < ~ 0.05),
                                   expression(p - value < ~ 0.05 ~ and ~ log[2] ~ FC > ~ 1)),
                  drawConnectors = TRUE, min.segment.length = 1, 
                  max.overlaps = 12, labSize = 4)
  ggsave("Output/Volcano_plots/Volcano plot - PL2-3+2DG vs PL2-3.png", 
         units = "px", width = 3000, height = 3000, dpi = 300)
  dev.off()
  
  EnhancedVolcano(vpData_R848, 
                  lab = vpData_R848$SYMBOL, 
                  pCutoff = 0.05, FCcutoff = 1,
                  x = 'log2FoldChange', y = 'padj', 
                  # xlim = c(-2.5, 1.5), 
                  ylim = c(0, 4),
                  title = "R848 + 2DG vs R848", 
                  subtitle = "(Adjusted p-values)", 
                  legendLabels = c("NS", expression(Log[2] ~ FC > 1), 
                                   expression(p - value < ~ 0.05),
                                   expression(p - value < ~ 0.05 ~ and ~ log[2] ~ FC > ~ 1)),
                  drawConnectors = TRUE, min.segment.length = 1, 
                  max.overlaps = 12, labSize = 4)
  ggsave("Output/Volcano_plots/Volcano plot - R848+2DG vs R848.png", 
         units = "px", width = 3000, height = 3000, dpi = 300)
  dev.off()

  EnhancedVolcano(vpData_PL23vR848, 
                  lab = vpData_PL23vR848$SYMBOL, 
                  pCutoff = 0.05, FCcutoff = 1,
                  x = 'log2FoldChange', y = 'padj', 
                  # xlim = c(-2.5, 1.5), 
                  ylim = c(0, 6),
                  title = "PL2-3 vs R848", 
                  subtitle = "(Adjusted p-values)", 
                  legendLabels = c("NS", expression(Log[2] ~ FC > 1), 
                                   expression(p - value < ~ 0.05),
                                   expression(p - value < ~ 0.05 ~ and ~ log[2] ~ FC > ~ 1)),
                  drawConnectors = TRUE, min.segment.length = 1, 
                  max.overlaps = 12, labSize = 4)
  ggsave("Output/Volcano_plots/Volcano plot - PL2-3 vs R848.png", 
         units = "px", width = 3000, height = 3000, dpi = 300)
  dev.off()
  
  EnhancedVolcano(vpData_2DG_PL23vR848, 
                  lab = vpData_2DG_PL23vR848$SYMBOL, 
                  pCutoff = 0.05, FCcutoff = 1,
                  x = 'log2FoldChange', y = 'padj', 
                  # xlim = c(-2.5, 1.5), 
                  ylim = c(0, 10),
                  title = "PL2-3 + 2DG vs R848 + 2DG", 
                  subtitle = "(Adjusted p-values)", 
                  legendLabels = c("NS", expression(Log[2] ~ FC > 1), 
                                   expression(p - value < ~ 0.05),
                                   expression(p - value < ~ 0.05 ~ and ~ log[2] ~ FC > ~ 1)),
                  drawConnectors = TRUE, min.segment.length = 1, 
                  max.overlaps = 12, labSize = 4)
  ggsave("Output/Volcano_plots/Volcano plot - PL2-3+2DG vs R848+2DG.png", 
         units = "px", width = 3000, height = 3000, dpi = 300)
  dev.off()
  