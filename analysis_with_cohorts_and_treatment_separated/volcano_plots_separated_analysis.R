setwd("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_and_treatment_separated")
getwd()

vp_data = list(
  PL23_2DG_v_Ctrl = deseq_res$AM14transfer_PL23[, c(2, 8, 12)],
  R848_2DG_v_Ctrl = deseq_res$AM14transfer_R848[, c(2, 8, 12)],
  AM14MRLlpr = deseq_res$AM14MRLlpr[, c(2, 8, 12)],
  B18 = deseq_res$B18transfer[, c(2, 8, 12)])
head(vp_data$AM14MRLlpr)


# confirm_foldchange(dds_AM14trans, "Crisp1", gene_IDs$AM14trans)
# confirm_foldchange(dds_B18trans, "Dmrt2", gene_IDs$B18trans)
# confirm_foldchange(dds_PL23vNP, "Lars2", gene_IDs$AM14trans)

volcano_wrapper(vp_data$PL23_2DG_v_Ctrl, "AM14 Transfer: PL2-3 + 2DG vs PL2-3",
                "Upregulated in PL2-3                                               Upregulated in PL2-3 + 2DG",
                0.6, 0.05, "AM14_Transfer-PL2-3_2DG_vs_Control")

volcano_wrapper(vp_data$R848_2DG_v_Ctrl, "AM14 Transfer: R848 + 2DG vs R848", 
                "Upregulated in R848                                                Upregulated in R848 + 2DG",
                0.6, 0.05, "AM14_Transfer-R848_2DG_vs_Control")

volcano_wrapper(vp_data$AM14MRLlpr, "AM14 MRL/lpr: 2DG vs Control",
                "Upregulated in Control                                                              Upregulated in 2DG",
                0.6, 0.05, "AM14_MRLlpr")

volcano_wrapper(vp_data$B18, "B1-8 Transfer: NP + 2DG vs NP",
                "Upregulated in NP                                                         Upregulated in NP + 2DG",
                0.6, 0.05, "B18_Transfer")


# Scatter plots correlating two datasets ---------------------------------------
  # Lists of potentially interesting genes -------------------------------------
    # highlight_genes <- c("Lck", "Tbx21", "Irf1"     
    #                      "Ccl4"      "Slc25a19"  "Cyp39a1"   "Slc37a1"   "Aif1"      "Aldh3b1"   "Cr2"      
    #                      "Batf3"     "Mrc1"      "Zbp1"      "Pim2"      "Plscr1"    "Tmem106a"  "Capn5"    
    #                      "C1qc"      "C1qb"      "Pml", "Rftn1"     "Dusp10"    "Spns2"     "Serpina3g"
    #                      "Irf8"      "Trafd1"    "Clec4a3"   "Irgm1"     "Fcgr4"     "Tnip2"     "Lyz1"     
    #                      "Wfdc17"    "Hpcal1"    "Ass1"      "Igtp"      "Ly6c1"     "Psmb9"    )
    
    abc_up_genes <- read.table("../public_GSE_datasets/GSE_lists/GSE99480_ABC_up_regulated_genes.txt", header = FALSE)$V1
    abc_down_genes <- read.table("../public_GSE_datasets/GSE_lists/GSE99480_ABC_down_regulated_genes.txt", header = FALSE)$V1

    jack_gmt_sets <- read.gmt("../public_GSE_datasets/Bcell_subset.human.jack.gmt")
    jack_gmt_sets <- jack_gmt_sets[jack_gmt_sets$gene != "", ]
    head(jack_gmt_sets)
    
    ABC_signature_in_Lupus_Wu <- jack_gmt_sets[jack_gmt_sets$term == "ABC signature in Lupus (ref: Wu et al.)", ]
    ABC_signature_in_Lupus_Wu <- ABC_signature_in_Lupus_Wu$gene
    ABC_signature_in_Lupus_Wu
    
  # Function to find significant genes in two data sets and make one table -----  
    find_DEG_set_overlap <- function(deg1, colname1, deg2, colname2, p_value_cutoff, lfc_cutoff){
      if(missing(colname1)) colname1 <- "x"
      if(missing(colname2)) colname2 <- "y"
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
        yes = "down",
        no = ifelse(
          test = merge_filt[, 2] > 0 & merge_filt[, 3] > 0, 
          yes = "up", 
          no = NA
        )
      )
      
      merge_filt$ABC_direction <- ifelse(
        test = merge_filt$Gene %in% abc_up_genes,
        yes = "ABC_up_FoB_down",
        no = ifelse(
          test = merge_filt$Gene %in% abc_down_genes,
          yes = "ABC_down_FoB_up",
          no = NA
        )
      )
      
      return(merge_filt)
    }
    

  # AM14 PL2-3 vs AM14 MRLlpr --------------------------------------------------
    pl23_x_MRLlpr <- find_DEG_set_overlap(vp_data$PL23_2DG_v_Ctrl, "PL23",
                                          vp_data$AM14MRLlpr, "MRLlpr")
    head(pl23_x_MRLlpr)

    ggplot(pl23_x_MRLlpr, aes(x = PL23, y = MRLlpr, fill = ABC_direction)) +
      geom_point(aes(color = same_direction)) +
      scale_color_manual(
        values = c("down" = "#00496f", "up" = "#F21A00", "neither" = "darkgray"),
        na.value = "gray",
        guide = NULL # guide = guide_legend(title = "title")
      ) +
      geom_hline(yintercept = c(0), linetype = "dotted", color = "#555555") +
      geom_vline(xintercept = c(0), linetype = "dotted", color = "#555555") +
      # geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", color = "#555555") +
      # geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "#555555") +
      geom_text_repel(
        data = pl23_x_MRLlpr[is.na(pl23_x_MRLlpr$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = 20
      ) +
      geom_label_repel(
        data = pl23_x_MRLlpr[!is.na(pl23_x_MRLlpr$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = Inf
      ) +      
      scale_fill_manual(
        values = c("ABC_up_FoB_down" = "#ed8b00", "ABC_down_FoB_up" = "#78B7C5"),
        na.value = "white"
      ) +
      labs(
        x = "LFC: AM14 PL2-3 2DG vs Ctrl",
        y = "LFC: AM14 MRL/lpr 2DG vs Ctrl",
        title = "AM14 PL2-3 +/- 2DG vs AM14 MRLlpr +/- 2DG"
      ) +
      theme_classic()
    
    ggsave(filename = "Output/Correlation Scatterplots/AM14 PL2-3 vs AM14 MRLlpr.png", 
           width = 8, height = 7, units = "in", dpi = 600)
    

  # AM14 PL2-3 vs AM14 R848 ----------------------------------------------------
    pl23_x_r848 <- find_DEG_set_overlap(vp_data$PL23_2DG_v_Ctrl, "PL23",
                                        vp_data$R848_2DG_v_Ctrl, "R848")
    head(pl23_x_r848)
    
    ggplot(pl23_x_r848, aes(x = PL23, y = R848, fill = ABC_direction)) +
      geom_point(aes(color = same_direction)) +
      scale_color_manual(
        values = c("down" = "#00496f", "up" = "#F21A00", "neither" = "darkgray"),
        na.value = "gray",
        guide = NULL # guide = guide_legend(title = "title")
      ) +
      geom_hline(yintercept = c(0), linetype = "dotted", color = "#555555") +
      geom_vline(xintercept = c(0), linetype = "dotted", color = "#555555") +
      # geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", color = "#555555") +
      # geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "#555555") +
      geom_text_repel(
        data = pl23_x_r848[is.na(pl23_x_r848$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = 100
      ) +
      geom_label_repel(
        data = pl23_x_r848[!is.na(pl23_x_r848$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = Inf
      ) +
      scale_fill_manual(
        values = c("ABC_up_FoB_down" = "#ed8b00", "ABC_down_FoB_up" = "#78B7C5"),
        na.value = "white"
      ) +
      labs(
        x = "LFC: AM14 PL2-3 2DG vs Ctrl",
        y = "LFC: AM14 R848 2DG vs Ctrl",
        title = "AM14 PL2-3 +/- 2DG vs AM14 R848 +/- 2DG"
      ) +
      # theme_bw()
      theme_classic()
    
    ggsave(filename = "Output/Correlation Scatterplots/AM14 PL2-3 vs AM14 R848.png", 
           width = 6, height = 4.5, units = "in", dpi = 600)
    
    
  # AM14 PL2-3 vs B1-8 NP-OVA --------------------------------------------------
    pl23_x_np <- find_DEG_set_overlap(vp_data$PL23_2DG_v_Ctrl, "PL23",
                                      vp_data$B18, "NP")
    head(pl23_x_np)
    
    ggplot(pl23_x_np, aes(x = PL23, y = NP, fill = ABC_direction)) +
      geom_point(aes(color = same_direction)) +
      scale_color_manual(
        values = c("down" = "#00496f", "up" = "#F21A00", "neither" = "darkgray"),
        na.value = "gray",
        guide = NULL # guide = guide_legend(title = "title")
      ) +
      geom_hline(yintercept = c(0), linetype = "dotted", color = "#555555") +
      geom_vline(xintercept = c(0), linetype = "dotted", color = "#555555") +
      # geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", color = "#555555") +
      # geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "#555555") +
      geom_text_repel(
        data = pl23_x_np[is.na(pl23_x_np$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = 100
      ) +
      geom_label_repel(
        data = pl23_x_np[!is.na(pl23_x_np$same_direction), ],
        aes(label = Gene),
        size = 3,
        max.overlaps = Inf
      ) +
      scale_fill_manual(
        values = c("ABC_up_FoB_down" = "#ed8b00", "ABC_down_FoB_up" = "#78B7C5"),
        na.value = "white"
      ) +
      labs(
        x = "LFC: AM14 PL2-3 2DG vs Ctrl",
        y = "LFC: B1-8 NP-OVA 2DG vs Ctrl",
        title = "AM14 PL2-3 +/- 2DG vs B1-8 NP-OVA +/- 2DG"
      ) +
      # theme_bw()
      theme_classic()
    
    ggsave(filename = "Output/Correlation Scatterplots/AM14 PL2-3 vs B1-8 NP-OVA.png", 
           width = 6, height = 4.5, units = "in", dpi = 600)
    