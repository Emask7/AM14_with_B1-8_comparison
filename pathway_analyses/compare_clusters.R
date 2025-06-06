# get_DEG_list <- function(dds_res, treatment, lfc_direction, lfc_cutoff){
#   if(missing(lfc_direction)) lfc_direction <- NULL
#   if(missing(lfc_cutoff)) lfc_cutoff <- 0.5
# 
#   if(is.null(lfc_direction)) {
#     sig <- dplyr::filter(dds_res, padj < 0.05 & abs(log2FoldChange) > abs(lfc_cutoff))
#   } else if(lfc_direction == "up"){
#     sig <- dplyr::filter(dds_res, padj < 0.05 & log2FoldChange > lfc_cutoff)
#   } else if(lfc_direction == "down"){
#     lfc_cutoff <- -1*abs(lfc_cutoff)
#     sig <- dplyr::filter(dds_res, padj < 0.05 & log2FoldChange < lfc_cutoff)
#   } else return(print("lfc_direction must be either up, down, or NULL"))
#   
#   if(nrow(sig) < 1) {
#     print("0 significant DEGs")
#     return(NULL)
#   }
#   sig <- as.character(sig$entrezgene)
#   treatment_df <- data.frame(rep(treatment, length(sig)))
#   
#   res <- data.frame(Genes = sig, Treatment = treatment_df)
#   
#   return(res)
# }
# 
# temp <- get_DEG_list(deseq_res$AM14MRLlpr, "AM14 MRLlpr")
# head(temp)
# nrow(temp)
# 
# ccdata <- rbind(get_DEG_list(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, "PL2-3 2DG vs PL2-3", NULL),
#             get_DEG_list(deseq_res$AM14transfer$R848_2DG_v_Ctrl, "R848 2DG vs R848", NULL),
#             # get_DEG_list(deseq_res$B18transfer, "NP + 2DG vs NP", NULL),
#             get_DEG_list(deseq_res$AM14MRLlpr, "AM14 MRLlpr", NULL))
# colnames(ccdata) <- c("Entrez", "group")
# head(ccdata)
# 
# cc_res <- compareCluster(Entrez~group, data=ccdata, fun='enrichGO', OrgDb='org.Mm.eg.db')
# head(cc_res)
# 
# write.csv(cc_res, "Output/Functional_analyses/CompareClusters_allDEGs.csv")
# 
# ggplot(cc_res, showCategory = 10,
#        aes(group, fct_reorder(Description, group))) +
#   # geom_segment(aes(xend=0, yend = Description)) +
#   geom_point(aes(color=p.adjust, size = RichFactor)) +
#   # geom_point(aes(color=RichFactor, size = Count)) +
#   scale_color_gradientn(colors = pnw_palette("Bay", 8, type = "continuous")) +
#   scale_size_continuous(range=c(3, 10)) +
#   theme_dose(12) +
#   # xlab("RichFactor") +
#   # xlab("FDR") +
#   ylab(NULL) +
#   # scale_y_discrete(labels = label_wrap(30)) +
#   ggtitle("CompareClusters")

top_terms <- function(go_res, group, n){
  if(nrow(go_res) < n) n <- nrow(go_res)
  
  # sig <- subset(go_res, go_res$p.adjust < 0.05)
  
  if(is.null(go_res)) {
    top_res <- data.frame(group, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
    colnames(top_res) <- c("group", "ID", "Description", "GeneRatio", "BgRatio", "RichFactor", "FoldEnrichment",
                            "zScore", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
  } else {
    top_res <- as.data.frame(go_res)
    top_res <- top_res[order(-top_res$RichFactor), ]
    top_res <- data.frame(group = rep(group, n), top_res[1:n, ])
  }
  
  # full_res <- data.frame(group = rep(group, nrow(go_res)), go_res)
  
  return(top_res)
}


test <- top_terms(enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified, 10, "PL2-3 + 2DG vs PL2-3")
test


# enrichGO_padj1.0_all <- list(
#   AM14trans = list(
#     PL23_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$PL23_2DG_v_Ctrl, NULL, 0.5, 1.0),
#     # R848_2DG_v_Ctrl = enrichGO_wrapper(deseq_res$AM14transfer$R848_2DG_v_Ctrl, NULL, 0.5, 1.0),
#     PL23_vs_R848 = enrichGO_wrapper(deseq_res$AM14transfer$PL23_vs_R848, NULL, 0.5, 1.0)),
#   # PL23_v_NP = enrichGO_wrapper(deseq_res$PL23_v_NP, NULL, 0.5, 1.0),
#   # B18trans = enrichGO_wrapper(deseq_res$B18transfer, NULL, 0.5, 1.0),
#   AM14MRLlpr = enrichGO_wrapper(deseq_res$AM14MRLlpr, NULL, 0.5, 1.0))
# 
# head(enrichGO_padj1.0_all$AM14trans$PL23_2DG_v_Ctrl$eGO)
# max(enrichGO_padj1.0_all$AM14trans$PL23_2DG_v_Ctrl$eGO$p.adjust)
# 
# top_GO <- list(
#   PL23_2DG_v_Ctrl = top_terms(enrichGO_padj1.0_all$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified, "PL2-3 2DG\nvs PL2-3", 10),
#   R848_2DG_v_Ctrl = top_terms(enrichGO_padj1.0_all$AM14trans$R848_2DG_v_Ctrl$eGO, "R848 2DG\nvs R848", 10),
#   AM14MRLlpr = top_terms(enrichGO_padj1.0_all$AM14MRLlpr$eGO_simplified, "AM14 MRLlpr\n2DG vs Ctrl", 10))


temp <- top_terms(enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified, "PL2-3 2DG\nvs PL2-3", 10)
head(temp$all)


ccdata <- rbind(top_terms(enrichGO_all$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified, "PL2-3 2DG\nvs PL2-3", 10),
                top_terms(enrichGO_all$AM14trans$R848_2DG_v_Ctrl$eGO, "R848 2DG\nvs R848", 10),
                top_terms(enrichGO_all$AM14MRLlpr$eGO_simplified, "AM14 MRLlpr\n2DG vs Ctrl", 10))
                # top_terms(enrichGO_all$B18trans, "NP + 2DG\nvs NP", 10))
tail(ccdata, 12)
colnames(ccdata)

ggplot(ccdata, 
       aes(group, fct_reorder(Description, group))) +
  geom_point(aes(color=RichFactor, size = -1*log10(p.adjust))) +
  scale_color_gradientn(colors = pnw_palette("Bay", 8, type = "continuous"),
                        guide = guide_colorbar(order = 1)) +
  scale_size_continuous(range=c(3, 12), 
                        guide = guide_legend(title = expression('-Log'[10]*'(FDR)'),
                                             order = 2)) +
  theme_dose(12) +
  xlab(NULL) +
  ylab(NULL) +
  scale_y_discrete(labels = label_wrap(35)) +
  ggtitle("CompareClusters")


ccdata_down <- rbind(top_terms(enrichGO_down$AM14trans$PL23_2DG_v_Ctrl$eGO_simplified, "PL2-3 2DG\nvs PL2-3", 10),
                     # top_terms(enrichGO_down$AM14trans$R848_2DG_v_Ctrl$eGO, "R848 2DG\nvs R848", 10),
                     top_terms(enrichGO_down$AM14MRLlpr$eGO_simplified, "AM14 MRLlpr\n2DG vs Ctrl", 10))

ggplot(ccdata_down, 
       aes(group, fct_reorder(Description, group))) +
  geom_point(aes(color=RichFactor, size = -1*log10(p.adjust))) +
  scale_color_gradientn(colors = pnw_palette("Bay", 8, type = "continuous"),
                        guide = guide_colorbar(order = 1)) +
  scale_size_continuous(range=c(3, 12), 
                        guide = guide_legend(title = expression('-Log'[10]*'(FDR)'),
                                             order = 2)) +
  theme_dose(12) +
  xlab(NULL) +
  ylab(NULL) +
  scale_y_discrete(labels = label_wrap(35)) +
  ggtitle("CompareClusters")

