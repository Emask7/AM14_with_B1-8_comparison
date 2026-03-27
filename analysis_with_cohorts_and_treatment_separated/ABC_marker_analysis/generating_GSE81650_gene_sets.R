# --- Step 1: Load the Required Libraries ---
library(limma)
library(dplyr)
library(tibble)

# --- Step 2: Load and Prepare Your Data ---
# Source paper: 
# Pubmed link: 

# Load the data from the CSV files
gse81650_data <- read.delim("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/public_GSE_datasets/GSE81650_series_matrix_edited.txt")
head(gse81650_data)

rownames(gse81650_data) <- gse81650_data$genes
gse81650_data <- gse81650_data[, 7:8]
colnames(gse81650_data) <- c("ABC_1_log2cpm", "ABC_2_log2cpm")

# --- Step 3: Filter for WT Samples and Combine ---
# First, ensure we only keep columns that contain numeric expression data.
gse81650_data <- gse81650_data %>%
  select_if(is.numeric)
head(gse81650_data)

fob_data <- fob_data %>%
  select_if(is.numeric)
head(fob_data)

# To safely combine them, we'll find the genes that are common to both datasets.
common_genes <- intersect(rownames(gse81650_data), rownames(fob_data))

# Now, create a single expression matrix with the common genes and all WT samples.
# The columns will be all gse81650_data samples followed by all fob_data samples.
expression_matrix <- cbind(gse81650_data[common_genes, ], fob_data[common_genes, ])

# --- Step 4: Create the Experimental Design ---
# We need to tell limma which sample belongs to which group.
# First, create a factor that defines the groups.
group <- factor(
  c(rep("ABC", ncol(gse81650_data)), rep("FoB", ncol(fob_data))),
  levels = c("FoB", "ABC") # This sets FoB as the baseline for comparison
)

# Next, create a design matrix. This is a numerical representation of your experiment.
design <- model.matrix(~group)
colnames(design) <- c("Intercept", "ABCvsFoB")

# --- Step 5: Run the Limma Analysis ---
# Fit the linear model for each gene based on the design matrix.
fit <- lmFit(expression_matrix, design)

# Apply empirical Bayes smoothing to the standard errors. This makes the
# analysis more powerful by borrowing information across all genes.
fit <- eBayes(fit)

# Extract a table of the top results. We specify the coefficient we are
# interested in, which is the comparison of ABC vs FoB.
# `n = Inf` ensures that results for all genes are returned.
results_table <- topTable(fit, coef = "ABCvsFoB", n = Inf)

# You can now view the top differential expressed genes
print("Top results from the differential expression analysis:")
head(results_table)

# --- Step 6: Create the Final Gene Sets ---
# Define the cutoffs for significance, as the original paper likely did.
# These are common values, but you can adjust them.
p_value_cutoff <- 0.05
logfc_cutoff <- 1 # A log2 fold-change of 1 means a 2-fold change in expression

# Filter the results table to get your "ABC up-regulated" genes
abc_up_regulated <- results_table %>%
  filter(adj.P.Val < p_value_cutoff & logFC > logfc_cutoff)

# Filter the results table to get your "ABC down-regulated" genes
abc_down_regulated <- results_table %>%
  filter(adj.P.Val < p_value_cutoff & logFC < -logfc_cutoff)

# Extract just the gene names from these tables
abc_up_genes <- rownames(abc_up_regulated)
abc_down_genes <- rownames(abc_down_regulated)

print("---")
print(paste("Found", length(abc_up_genes), "up-regulated genes."))
print(paste("Found", length(abc_down_genes), "down-regulated genes."))

# --- Step 7 (Optional): Save Gene Sets to Files ---
# This creates simple text files, perfect for input into GSEA software.

write.table(abc_up_genes, file = "Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_and_treatment_separated/ABC_marker_analysis/GSE_lists/GSE99480_ABC_up_regulated_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(abc_down_genes, file = "Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/analysis_with_cohorts_and_treatment_separated/ABC_marker_analysis/GSE_lists/GSE99480_ABC_down_regulated_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

print("Gene sets have been saved to .txt files in your current R working directory.")