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
