
find_ensembl_ID <- function(search_name) {
  dat <- dplyr::filter(gene_IDs$AM14trans, external_gene_name == search_name)
  dat$ensembl_gene_id
}

get_counts_AM14trans <- function(gene_name){
  ID <- find_ensembl_ID(gene_name)
  
  res_PL23 <- counts(dds_AM14trans_PL23, normalized = TRUE)
  res_PL23 <- res_PL23[rownames(res_PL23) %in% ID, ]
  
  res_R848 <- counts(dds_AM14trans_R848, normalized = TRUE)
  res_R848 <- res_R848[rownames(res_R848) %in% ID, ]
  
  df <- data.frame(res_PL23[1:3], res_PL23[4:6], res_R848[1:3], res_R848[4:6])
  colnames(df) <- c("PL2-3", "PL2-3+2DG", "R848", "R848+2DG")
  rownames(df) <- c("a", "b", "c")
  df
}

get_counts_AM14MRLlpr <- function(gene_name){
  ID <- find_ensembl_ID(gene_name)
  
  res_AM14MRLlpr <- counts(dds_AM14MRLlpr, normalized = TRUE)
  res_AM14MRLlpr <- res_AM14MRLlpr[rownames(res_AM14MRLlpr) %in% ID, ]
  
  df <- data.frame(c(res_AM14MRLlpr[1:6], "NA"), res_AM14MRLlpr[7:13])
  colnames(df) <- c("Control", "2DG")
  rownames(df) <- c("a", "b", "c", "d", "e", "f", "g")
  df
}

get_counts_B18trans <- function(gene_name){
  ID <- find_ensembl_ID(gene_name)
  
  res_B18 <- counts(dds_B18trans, normalized = TRUE)
  res_B18 <- res_B18[rownames(res_B18) %in% ID, ]
  
  df <- data.frame(c(res_B18[6:9], "NA"), res_B18[1:5])
  colnames(df) <- c("NP", "NP+2DG")
  rownames(df) <- c("a", "b", "c", "d", "e")
  df
}



# Complement genes
  get_counts("C1qa")
  get_counts("C1qb")
  get_counts("C1qc")
  get_counts("Fcna")
  get_counts("C6")
  get_counts("C4b")
  get_counts("C3")
  get_counts("Cfh")
  get_counts("Ighg2c")
  get_counts("Ighm")
  get_counts("Ighg2b")
  
# TLR signaling
  get_counts("Unc93b1")
  get_counts("Tlr9")
  get_counts("Ptpn22")
  get_counts("Treml4")
  get_counts("Rab7b")
  get_counts("Colec12")
  get_counts("Cxcr4")
  get_counts("Gpr55")
  
# Th1 pathway
  get_counts_AM14trans("Tbx21")
  get_counts_AM14MRLlpr("Tbx21")
  get_counts_B18trans("Tbx21")
  
  get_counts_AM14trans("Il12a")
  get_counts_AM14MRLlpr("Il12a")
  get_counts_B18trans("Il12a")
  


normcounts_plot <- function(gene_list, dds, ID_df) {
  cts <- counts(dds, normalized = TRUE)
  cts <- data.frame(rownames(cts), cts)
  colnames(cts) <- c("ensembl_gene_id", dds$barplot_info)
  cts <- right_join(ID_df[, 1:2], cts, by = "ensembl_gene_id")
  cts <- cts[, c(2:ncol(cts))]
  
  plot_data <- cts %>%
    pivot_longer(
      cols = -(1), # by using a `-`, we can exclude the first two columns
      values_to = c("normCounts"),
      names_to = c("stim", "drug", "mouse"),
      names_sep  = "_"
    ) %>%
    group_by(external_gene_name, stim, drug) %>%
    summarise(n = n(),
              reads = normCounts,
              mean = mean(normCounts),
              sem = std.error(normCounts),
              sd = sd(normCounts)) %>%
    filter(external_gene_name %in% gene_list)

  # plot_data %>%
  #   ggplot(aes(x = stim, y = mean, fill = drug)) +
  #   geom_bar(position = position_dodge(width = 0.9), stat = "identity") +
  #   ggtitle("title") +
  #   geom_errorbar(
  #     aes(x = stim, ymin = mean - sem, ymax = mean + sem),
  #     width = 0.2,
  #     position = position_dodge(width = 0.9),
  #     stat = "identity"
  #   ) +
  #   geom_jitter(aes(y = reads), position = position_dodge(width = 0.9)) +
  #   facet_wrap( ~ external_gene_name, scales = 'free')
  
  return(plot_data)
}


normcounts_plot(c("Unc93b1"), dds_AM14MRLlpr, gene_IDs$AM14MRLlpr) 
p1 <- normcounts_plot(c("Unc93b1"), dds_AM14trans, gene_IDs$AM14trans) 
head(p1)
p1

df.aov <- p1 %>%
  group_by(external_gene_name) %>%
  nest() %>%
  mutate(aov = map(data, ~aov(reads ~ stim * drug, data = .x)))
head(df.aov)
coefficients(df.aov$aov[[1]])


normcounts_plot(c("Unc93b1", "C1qa", "C1qb", "C1qc"), dds_AM14trans, gene_IDs$AM14trans) 

p1 <- normcounts_plot(c("Unc93b1", "C1qa", "C1qb", "C1qc"), dds_AM14trans, gene_IDs$AM14trans) 
p1
p1 %>%
  ggplot(aes(x = treatment, y = mean, fill = treatment)) +
  geom_bar(position = position_dodge(width = 0.9), stat = "identity") +
  ggtitle("title") +
  geom_errorbar(
    aes(x = treatment, ymin = mean - sem, ymax = mean + sem),
    width = 0.2, 
    position = position_dodge(width = 0.9),
    stat = "identity"
  ) +
  geom_jitter(aes(y = reads), position = position_dodge(width = 0.9)) +
  facet_wrap( ~ external_gene_name, scales = 'free')


bp <- tidy_reads("Unc93b1", "AM14trans")
bp$summary
bp$df
head(bp$cts)
head(bp$cts_with_symbols)

test <- bp$cts_with_symbols %>%
  pivot_longer(
    cols = -(1:2), # by using a `-`, we can exclude the first two columns
    values_to = c("normCounts"),
    names_to = c("treatment", "mouse"),
    names_sep  = "_"
  )
head(test)
rm(test)

test %>%
  filter(external_gene_name == "Unc93b1") %>%
  ggplot(aes(x = treatment)) +
  geom_bar() +
  ggtitle("geom_bar")

# temp_cts <- counts(dds_AM14trans, normalized = TRUE)
# rm(temp_cts)


# summary1 <-
#   bp$df %>%
#   group_by(gene, treatment) %>%
#   summarise(n = n(),
#             mean = mean(normCounts),
#             sd = sd(normCounts))
# 
#   summary1 %>% knitr::kable()
  