rm(list = ls())
install_if_missing <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      # First try installing from CRAN
      tryCatch({
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
      }, error = function(e) {
        # If CRAN installation fails, try Bioconductor
        message(paste("Attempting to install", pkg, "from Bioconductor..."))
        BiocManager::install(pkg)
        library(pkg, character.only = TRUE)
      })
    }
  }
}



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


#install_if_missing(packages_list)
#install_if_missing(c("org.Hs.eg.db"))
rm(install_if_missing)
# Load the packages

for (pkg in packages_list) {
  library(pkg, character.only = TRUE)
}

rm(list = ls())

data("airway")

airway

# Extract count data and metadata
count_data <- assay(airway)
metadata <- colData(airway)

head(count_data)
dim(count_data)

# Explore metadata
sample_info <- as.data.frame(metadata)
sample_info = sample_info[, c(2, 3)]
sample_info$dex <- gsub("trt", "treated", sample_info$dex)
sample_info$dex <- gsub("untrt", "untreated", sample_info$dex)
names(sample_info) <- c("cell_line", "dexamethasone")

# Check the modified metadata
head(sample_info)

# Summary of sample metadata
summary(sample_info)

# Check for missing values in count data
sum(is.na(count_data))


# check the count data matrix
all(colnames(count_data) %in% rownames(sample_info))


all(colnames(count_data) == rownames(sample_info))


dex_init <- DESeqDataSetFromMatrix(countData = count_data,colData = sample_info,design = ~dexamethasone)


# pre filtering rows that has low count of genes
dex = dex_init[rowSums(counts(dex_init))>=3,]

dex = DESeq(dex)

result = results(dex)


# View the summary of results
summary(result)

result_0_05 <- results(dex, alpha = 0.05)

result_0_001 <- results(dex, alpha = 0.001)

result_0_05_sorted <- result_0_05[order(result_0_05$padj), ]
result_0_001_sorted <- result_0_001[order(result_0_001$padj), ]


top_10_genes_0_05 <- head(result_0_05_sorted, 10)
top_10_genes_0_001 <- head(result_0_001_sorted, 10)

top_10_genes_0_05
top_10_genes_0_001

plotMA(result)
plotMA(result_0_05)
plotMA(result_0_001)


# genes enrich
ens_to_gene <- function(ensembl_id){
  return(mapIds(org.Hs.eg.db, keys = ensembl_id, column = "SYMBOL", keytype = "ENSEMBL"))
}
gene_names_0_05 = ens_to_gene(rownames(top_10_genes_0_05))
gene_names_0_001 = ens_to_gene(rownames(top_10_genes_0_001))

rownames(top_10_genes_0_05)<-gene_names_0_05
rownames(top_10_genes_0_001)<-gene_names_0_001


enrich_genes <- function(gene_list, databases){
  enr <- enrichr(gene_list, databases)
  for (i in seq_along(enr)){
    cat("\nDatabase:",names(enr)[i],"\n")
    print(head(enr[[i]]),10)
  }
  return(enr)
}

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
if (websiteLive) dbs <- listEnrichrDbs()

library_names <- c("KEGG_2021_Human", "Reactome_2022", "GO_Biological_Process_2023", 
                   "GO_Molecular_Function_2023", "GO_Cellular_Component_2023", 
                   "OMIM_Disease", "OMIM_Expanded", "DisGeNET", "Human_Phenotype_Ontology")

enriched_result_0_05 = enrich_genes(gene_names_0_05,library_names)
enriched_result_0_001 = enrich_genes(gene_names_0_001,library_names)


# box plot count
boxplot(count_data, 
        main = "Boxplot of Gene Expression Data", 
        ylab = "Expression Levels",
        col = "lightblue", # Optional: Add color to the boxplots
        outline = FALSE, # Removes outliers to make the plot cleaner
        names = colnames(count_data), # Add sample names on the x-axis
        las = 2) # Rotate x-axis labels 90 degrees

# bar plot count
mean_expression <- colMeans(count_data)
max_value <- max(mean_expression) * 1.1 # Increase by 10% to leave some space above the tallest bar

barplot(mean_expression, 
        main = "Barplot of Mean Expression Levels Across Samples", 
        ylab = "Mean Expression Levels",
        col = "lightgreen", 
        names.arg = colnames(count_data), 
        las = 2, 
        ylim = c(0, max_value)) # Set y-axis limits

# heatmap 10 top genes 

top_deg_0_05 = head(order(result_0_05$padj),10)
top_deg_0_001 = head(order(result_0_001$padj),10)

heatmap_data_0_05 = assay(dex)[top_deg_0_05,]
row.names(heatmap_data_0_05) <- rownames(top_10_genes_0_05)

heatmap_data_0_001 = assay(dex)[top_deg_0_001,]
row.names(heatmap_data_0_001) <- rownames(top_10_genes_0_001)


pheatmap(heatmap_data_0_05,cluster_rows = T,cluster_cols = T,annotation_col =  sample_info)
pheatmap(heatmap_data_0_001,cluster_rows = T,cluster_cols = T,annotation_col =  sample_info)

# scatter plot two top genes 
rld <- rlog(dex, blind = FALSE)  # or vst(dex, blind = FALSE)

# Choose two genes from the top differentially expressed genes
gene1_symbol <- "SPARCL1"
gene2_symbol <- "STOM"

# Find Ensembl IDs for the gene symbols
gene1_ensembl <- names(which(ens_to_gene(rownames(dex)) == gene1_symbol))[1]
gene2_ensembl <- names(which(ens_to_gene(rownames(dex)) == gene2_symbol))[1]

# Extract normalized expression values for the two genes using Ensembl IDs
gene1_counts <- assay(rld)[gene1_ensembl, ]
gene2_counts <- assay(rld)[gene2_ensembl, ]

scatter_df <- data.frame(
  Gene1 = gene1_counts,
  Gene2 = gene2_counts
)
ggplot(scatter_df, aes(x = Gene1, y = Gene2)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a linear regression line
  theme_minimal() +
  labs(
    title = paste("Scatter Plot of", gene1_symbol, "vs", gene2_symbol),
    x = paste(gene1_symbol, "Expression (rlog)"),
    y = paste(gene2_symbol, "Expression (rlog)")
  )

cor_test <- cor.test(scatter_df$Gene1, scatter_df$Gene2, method = "pearson")
cat("Pearson Correlation Coefficient:", cor_test$estimate, "\n")
cat("P-value:", cor_test$p.value, "\n")

# correlation analysis of two genes 




