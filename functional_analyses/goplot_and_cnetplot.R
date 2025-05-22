interesting_GO_PL23_2DG_v_Ctrl <- read.xlsx("Output/Functional_analyses/enrich_GOBP_results_allDEGs.xlsx",
                                            sheet = "AM14 PL23_2DG_vs_Ctrl - Simp")
interesting_GO_PL23_2DG_v_Ctrl <- interesting_GO_PL23_2DG_v_Ctrl$Description
head(interesting_GO_PL23_2DG_v_Ctrl)

lfc_data <- deseq_res$AM14transfer$PL23_2DG_v_Ctrl
lfc_data <- dplyr::filter(lfc_data, padj < 0.05 & abs(log2FoldChange) > 0.5)
lfc_names <- lfc_data$external_gene_name
lfc_data <- lfc_data[, 6]
names(lfc_data) <- lfc_names
head(lfc_data)

# test <- enrichGO_down$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified
# test <- sample(test$Description, 10)
# test

goplot(enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified)

cnetplot(enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified,
         showCategory = interesting_GO_PL23_2DG_v_Ctrl[1:5],
         # node_label = "share",
         # node_label = "exclusive",
         foldChange = lfc_data) #+
  # scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous"))

# cnetplot(enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified,
#          showCategory = interesting_GO_PL23_2DG_v_Ctrl[c(1:10)],
#          # node_label = "share",
#          # node_label = "exclusive",
#          foldChange = lfc_data)
# 
# cnetplot(enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified,
#          showCategory = interesting_GO_PL23_2DG_v_Ctrl[c(1:5)],
#          # node_label = "share",
#          # node_label = "exclusive",
#          foldChange = lfc_data)
# 
# cnetplot(enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified,
#          showCategory = interesting_GO_PL23_2DG_v_Ctrl[c(2, 5:10)],
#          # node_label = "share",
#          node_label = "exclusive",
#          foldChange = lfc_data)
# 
# cnetplot(enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified,
#          showCategory = interesting_GO_PL23_2DG_v_Ctrl[c(2, 5:7, 8)],
#          # node_label = "share",
#          # node_label = "exclusive",
#          foldChange = lfc_data)



