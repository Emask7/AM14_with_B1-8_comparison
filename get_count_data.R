
find_ensembl_ID <- function(search_name) {
  dat <- dplyr::filter(gene_IDs$AM14trans, external_gene_name == search_name)
  dat$ensembl_gene_id
}

get_counts <- function(gene_name){
  ID <- find_ensembl_ID(gene_name)
  
  res <- counts(dds_AM14trans, normalized = TRUE)
  
  res <- res[rownames(res) %in% ID, ]
  # res
  
  df <- data.frame(res[1:5], res[6:10], res[11:15], res[16:20])
  colnames(df) <- c("PL2-3", "PL2-3+2DG", "R848", "R848+2DG")
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
  