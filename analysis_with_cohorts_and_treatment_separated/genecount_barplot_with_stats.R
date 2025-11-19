setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_and_treatment_separated/")
getwd()

# Functions are in "data_visualization_functions.R" script 

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
