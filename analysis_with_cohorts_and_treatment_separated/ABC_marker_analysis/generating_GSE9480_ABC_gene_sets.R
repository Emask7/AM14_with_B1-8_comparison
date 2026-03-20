if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
install.packages("dplyr")
install.packages("tibble")

# --- Step 1: Load the Required Libraries ---
library(limma)
library(dplyr)
library(tibble)

# --- Step 2: Load and Prepare Your Data ---

# !!! IMPORTANT: Replace with the actual file paths to your unzipped .csv files !!!
abc_file_path <- "path/to/your/GSE99480_ABC_WT_DKO_09_04_2017.csv"
fob_file_path <- "path/to/your/GSE99480_FoB_WT_DKO_ABC_22_07_2016.csv"

# Load the data from the CSV files
abc_data <- read.csv("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/GSE99480/GSE99480_ABC_WT_DKO_09_04_2017.csv")
abc_data <- abc_data[ , 1:8]
rownames(abc_data) <- abc_data$genes
head(abc_data)

fob_data <- read.csv("Z:/Emma Mask/R_projects/AM14_with_B1-8_comparison/GSE99480/GSE99480_FoB_WT_DKO_ABC_22_07_2016.csv")
fob_data <- fob_data[ , c(1:8)]
rownames(fob_data) <- fob_data$genes
head(fob_data)

# --- Step 3: Filter for WT Samples and Combine ---
# First, ensure we only keep columns that contain numeric expression data.
abc_numeric <- abc_data %>%
  select_if(is.numeric)

fob_numeric <- fob_data %>%
  select_if(is.numeric)

# Now, from the numeric-only data, select the columns that contain "WT" (Wild Type).
abc_wt <- abc_numeric %>%
  select(contains("WT"))

fob_wt <- fob_numeric %>%
  select(contains("WT"))

colnames(abc_wt) <- c("ABC_1_log2cpm", "ABC_2_log2cpm")
colnames(fob_wt) <- c("FoB_1_log2cpm", "FoB_2_log2cpm")

# To safely combine them, we'll find the genes that are common to both datasets.
common_genes <- intersect(rownames(abc_wt), rownames(fob_wt))

# Now, create a single expression matrix with the common genes and all WT samples.
# The columns will be all ABC_WT samples followed by all FoB_WT samples.
expression_matrix <- cbind(abc_wt[common_genes, ], fob_wt[common_genes, ])

# --- Step 4: Create the Experimental Design ---

# We need to tell limma which sample belongs to which group.
# First, create a factor that defines the groups.
group <- factor(
  c(rep("ABC", ncol(abc_wt)), rep("FoB", ncol(fob_wt))),
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

# You can now view the top differentially expressed genes
print("Top results from the differential expression analysis:")
head(results_table)

# --- Step 6: Create the Final Gene Sets ---
# Define the cutoffs for significance, as the original paper likely did.
# These are common values, but you can adjust them.
p_value_cutoff <- 0.05
logfc_cutoff <- 1.0 # A log2 fold-change of 1 means a 2-fold change in expression

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

write.table(abc_up_genes, file = "analysis_with_cohorts_and_treatment_separated/ABC_marker_analysis/GSE99480_ABC_up_regulated_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(abc_down_genes, file = "analysis_with_cohorts_and_treatment_separated/ABC_marker_analysis/GSE99480_ABC_down_regulated_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

print("Gene sets have been saved to .txt files in your current R working directory.")