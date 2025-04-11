# Bulk RNA-seq Differential Expression Analysis
# Main analysis script

# Clean workspace
rm(list = ls())

# Package installation helper
install_if_missing <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      tryCatch({
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
      }, error = function(e) {
        message(paste("Attempting to install", pkg, "from Bioconductor..."))
        BiocManager::install(pkg)
        library(pkg, character.only = TRUE)
      })
    }
  }
}

# Required packages
packages_list = c(
  "GEOquery",
  "Seurat",
  "dplyr",
  "tidyverse",
  "airway",
  "DESeq2",
  "ggplot2",
  "patchwork",
  "enrichR",
  "pheatmap",
  "org.Hs.eg.db"
)

# Uncomment to install packages
# install_if_missing(packages_list)
rm(install_if_missing)

# Load packages
for (pkg in packages_list) {
  library(pkg, character.only = TRUE)
}

rm(list = ls())

# Load data
data("airway")
airway

# Extract count data and metadata
count_data <- assay(airway)
metadata <- colData(airway)

head(count_data)
dim(count_data)

# Prepare metadata
sample_info <- as.data.frame(metadata)
sample_info = sample_info[, c(2, 3)]
sample_info$dex <- gsub("trt", "treated", sample_info$dex)
sample_info$dex <- gsub("untrt", "untreated", sample_info$dex)
names(sample_info) <- c("cell_line", "dexamethasone")

head(sample_info)
summary(sample_info)

# Check for missing values
sum(is.na(count_data))

# Validate data structure
all(colnames(count_data) %in% rownames(sample_info))
all(colnames(count_data) == rownames(sample_info))

message("Data loaded and validated successfully.")

# Create DESeq2 dataset
dex_init <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sample_info,
  design = ~dexamethasone
)

# Pre-filtering: remove genes with low counts
dex = dex_init[rowSums(counts(dex_init)) >= 3, ]

# Run DESeq2 analysis
dex = DESeq(dex)

# Get results
result = results(dex)
summary(result)

# Get results at different significance thresholds
result_0_05 <- results(dex, alpha = 0.05)
result_0_001 <- results(dex, alpha = 0.001)

# Sort by adjusted p-value
result_0_05_sorted <- result_0_05[order(result_0_05$padj), ]
result_0_001_sorted <- result_0_001[order(result_0_001$padj), ]

# Extract top 10 genes
top_10_genes_0_05 <- head(result_0_05_sorted, 10)
top_10_genes_0_001 <- head(result_0_001_sorted, 10)

top_10_genes_0_05
top_10_genes_0_001

# Generate MA plots
plotMA(result)
plotMA(result_0_05)
plotMA(result_0_001)

# Gene annotation: Convert Ensembl IDs to gene symbols
ens_to_gene <- function(ensembl_id) {
  return(mapIds(org.Hs.eg.db, keys = ensembl_id, column = "SYMBOL", keytype = "ENSEMBL"))
}

gene_names_0_05 = ens_to_gene(rownames(top_10_genes_0_05))
gene_names_0_001 = ens_to_gene(rownames(top_10_genes_0_001))

rownames(top_10_genes_0_05) <- gene_names_0_05
rownames(top_10_genes_0_001) <- gene_names_0_001

message("DESeq2 analysis complete. Gene symbols annotated.")

