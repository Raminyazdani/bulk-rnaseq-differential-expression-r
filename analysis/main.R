# Bulk RNA-seq Differential Expression Analysis
# Main analysis script

# Clean workspace
rm(list = ls())

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

# Load packages
for (pkg in packages_list) {
  library(pkg, character.only = TRUE)
}

# TODO: Add analysis code
message("Analysis script initialized. Ready to add analysis code.")
