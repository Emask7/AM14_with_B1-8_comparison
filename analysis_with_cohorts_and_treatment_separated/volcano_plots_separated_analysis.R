setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_and_treatment_separated")
getwd()

vp_data = list(
  PL23_2DG_v_Ctrl = deseq_res$AM14transfer_PL23[, c(2, 8, 12)],
  R848_2DG_v_Ctrl = deseq_res$AM14transfer_R848[, c(2, 8, 12)],
  AM14MRLlpr = deseq_res$AM14MRLlpr[, c(2, 8, 12)],
  B18 = deseq_res$B18transfer[, c(2, 8, 12)])
head(vp_data$AM14MRLlpr)

vp_PL23 <- volcano_wrapper(vp_data$PL23_2DG_v_Ctrl, sub_title = "", lfc_cutoff = 0.6)
vp_PL23   <- ggdraw(vp_PL23) + 
  draw_label("A",  x = 0.02, y = 0.96, hjust = 0, fontface = "bold", size = 30) +
  draw_label("AM14 Transfer: PL2-3 + 2DG vs PL2-3",  x = 0.1, y = 0.93, hjust = 0, fontface = "bold", size = 22) +
  draw_label("Upregulated in PL2-3",  x = 0.1, y = 0.04, hjust = 0, fontface = "italic", size = 18) +
  draw_label("Upregulated in PL2-3 + 2DG",  x = 0.98, y = 0.04, hjust = 1, fontface = "italic", size = 18)

vp_R848 <- volcano_wrapper(vp_data$R848_2DG_v_Ctrl, sub_title = "", lfc_cutoff = 0.6)
vp_R848   <- ggdraw(vp_R848) + 
  draw_label("B",  x = 0.02, y = 0.96, hjust = 0, fontface = "bold", size = 30) +
  draw_label("AM14 Transfer: R848 + 2DG vs R848",  x = 0.1, y = 0.93, hjust = 0, fontface = "bold", size = 22) +
  draw_label("Upregulated in R848",  x = 0.1, y = 0.04, hjust = 0, fontface = "italic", size = 18) +
  draw_label("Upregulated in R848 + 2DG",  x = 0.98, y = 0.04, hjust = 1, fontface = "italic", size = 18)


vp_B18 <- volcano_wrapper(vp_data$B18, sub_title = "", lfc_cutoff = 0.6)
vp_B18    <- ggdraw(vp_B18) + 
  draw_label("C",  x = 0.02, y = 0.96, hjust = 0, fontface = "bold", size = 30) +
  draw_label("B1-8 Transfer: NP-OVA + 2DG vs NP-OVA",  x = 0.1, y = 0.93, hjust = 0, fontface = "bold", size = 22) +
  draw_label("Upregulated in NP-OVA",  x = 0.1, y = 0.04, hjust = 0, fontface = "italic", size = 18) +
  draw_label("Upregulated in NP-OVA + 2DG", x = 0.98, y = 0.04, hjust = 1, fontface = "italic", size = 18)


vp_MRLlpr <- volcano_wrapper(vp_data$AM14MRLlpr, sub_title = "", lfc_cutoff = 0.6)
vp_MRLlpr <- ggdraw(vp_MRLlpr) + 
  draw_label("D",  x = 0.02, y = 0.96, hjust = 0, fontface = "bold", size = 30) +
  draw_label("AM14 MRL/lpr: 2DG vs Control",  x = 0.1, y = 0.93, hjust = 0, fontface = "bold", size = 22) +
  draw_label("Upregulated in Control",  x = 0.1, y = 0.04, hjust = 0, fontface = "italic", size = 18) +
  draw_label("Upregulated in 2DG", x = 0.98, y = 0.04, hjust = 1, fontface = "italic", size = 18)

combined <- (vp_PL23 | vp_R848) / (vp_B18 | vp_MRLlpr) 

ggsave("Output/Volcano_plots/Combined_volcano.png",
       plot = combined,
       width = 20, height = 20, dpi = 600, units = "in",
       bg = "white")

rm(vp_PL23, vp_R848, vp_B18, vp_MRLlpr)
