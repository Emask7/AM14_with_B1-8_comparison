mk_DEG_list <- function(dat, lfc_cutoff) {
  dat <- subset(dat, padj < 0.05 & abs(log2FoldChange) >= lfc_cutoff)
  gene_list <- dat$log2FoldChange
  names(gene_list) <- dat$entrezgene
  gene_list <- subset(gene_list, !is.na(names(gene_list)))
  return(sort(gene_list, decreasing = TRUE))
}

pathview_wrap <- function(dat, pw_ID){
  pathview(gene.data = dat, gene.idtype = "entrez", pathway.id = pw_ID, 
           species = "mmu", limit = list(gene=max(abs(dat)), cpd=1))
}

KEGG_dotplot <- function(kegg_dat, plot_title, name_type, NES_min, NES_max, ignore_limit){
  if(missing(ignore_limit)) ignore_limit <- FALSE
  if(nrow(kegg_dat) > 35 & ignore_limit == FALSE){
    kegg_dat <- kegg_dat[order(kegg_dat$FDR.q.val, decreasing = FALSE), ]
    kegg_dat <- kegg_dat[1:35, ]
    plot_title <- stri_join(c(plot_title, "(Top 35)"), collapse = "\n")
  }
  
  if(max(kegg_dat$FDR.q.val) <= 0.01) x_max <- max(kegg_dat$FDR.q.val)
  else if(max(kegg_dat$FDR.q.val) <= 0.05) x_max <- 0.05
  else x_max <- max(kegg_dat$FDR.q.val)
  
  if(missing(name_type)) name_type <- "na"
  else if(name_type == "full") kegg_dat$NAME <- kegg_dat$FULL_NAME
  
  ggplot(kegg_dat,
         aes(x = FDR.q.val, y = reorder(NAME, FDR.q.val, decreasing = TRUE),
             size = SIZE, color = NES)) +
    geom_point() +
    scale_size_area(max_size = 5, limits = c(1, max(kegg_dat$SIZE))) +
    scale_colour_gradient2(name = "NES", 
                           low = "blue", mid = "white", high = "red", 
                           limits = c(NES_min, NES_max)) +
    xlab("FDR") +
    scale_x_continuous(limits = c(0.00, x_max)) +
    ylab("") +
    ggtitle(plot_title) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
}


# PL2-3 + 2DG vs PL2-3 ---------------------------------------------------------
  setwd("./KEGG_pathways/PL23_2DG_vs_PL23")
  PL23_kegg_data <- mk_DEG_list(res_PL23_full, 1)
  head(PL23_kegg_data)
  
  kegga_PL23 <- kegga(PL23_kegg_data$entrezgene, 
                      universe = res_PL23_full$entrezgene,
                      species = "Mm")
  head(kegga_PL23)
  topKEGG(kegga_PL23, n = 22)
  nrow(subset(kegga_PL23, P.DE < 0.05))
  
  
  # Complement and coagulation cascades
    pathview_wrap(PL23_kegg_data, 'mmu04610')
    
    
# Reset the working directory --------------------------------------------------
  setwd("./../..")
