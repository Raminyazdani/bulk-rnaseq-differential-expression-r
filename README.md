# Bulk RNA-seq Differential Expression Analysis

> Automated differential expression analysis pipeline for bulk RNA-seq data using DESeq2 and pathway enrichment

[![R](https://img.shields.io/badge/R-4.x-blue.svg)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.18-green.svg)](https://bioconductor.org/)
[![DESeq2](https://img.shields.io/badge/DESeq2-analysis-orange.svg)](https://bioconductor.org/packages/DESeq2/)

## Overview

A reproducible R-based pipeline for bulk RNA-seq differential expression analysis. This project demonstrates a complete workflow from count data to biological insights, featuring DESeq2 normalization, statistical testing, multi-database pathway enrichment, and comprehensive visualization of gene expression patterns.

**Key capabilities:**
- Robust differential expression analysis with DESeq2
- Pathway enrichment across 9 biological databases
- Comprehensive visualization suite (MA plots, heatmaps, correlation analysis)
- Quality control and expression distribution analysis
- Fully automated gene annotation and enrichment reporting

## Problem & Approach

**Problem**: Analyzing bulk RNA-sequencing data to identify differentially expressed genes between experimental conditions requires normalization, statistical testing, and biological interpretation through pathway analysis.

**Approach**: This pipeline uses DESeq2 for robust differential expression analysis with negative binomial modeling, followed by pathway enrichment using multiple databases (KEGG, Reactome, GO, OMIM, DisGeNET, HPO) via enrichR. The workflow includes quality control, multiple significance thresholds, and correlation analysis.

## Tech Stack

- **R** (4.x) - Statistical computing environment
- **Bioconductor** - Suite of bioinformatics packages
- **DESeq2** - Differential expression analysis
- **airway** - Example RNA-seq dataset
- **enrichR** - Pathway and ontology enrichment
- **ggplot2** - Advanced visualization
- **pheatmap** - Heatmap generation
- **org.Hs.eg.db** - Human gene annotation

## Repository Structure

```
bulk-rnaseq-differential-expression-r/
├── analysis/                    # Analysis scripts and project files
│   ├── analysis.Rproj          # RStudio project file
│   ├── main.R                  # Main analysis pipeline
│   ├── .RData                  # R workspace (gitignored)
│   └── .Rhistory               # R history (gitignored)
├── bulkRNA_project.docx        # Additional documentation
├── .gitignore                  # Git ignore rules for R
├── project_identity.md         # Project metadata
└── README.md                   # This file
```

## Setup & Installation

### Prerequisites
- R version 4.0 or higher
- (Optional) RStudio for interactive analysis

### Install Required Packages

Run the following in R to install all dependencies:

```r
# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c(
    "DESeq2",
    "airway",
    "org.Hs.eg.db"
))

# Install CRAN packages
install.packages(c(
    "GEOquery",
    "Seurat",
    "dplyr",
    "tidyverse",
    "ggplot2",
    "patchwork",
    "enrichR",
    "pheatmap"
))
```

**Note**: The installation script is included in `main.R` but commented out. Uncomment lines 40-41 to auto-install on first run.

## How to Run

### Option 1: RStudio (Recommended)
1. Open `analysis/analysis.Rproj` in RStudio
2. Open `main.R` in the editor
3. Run the entire script: `Ctrl+Alt+R` (Windows/Linux) or `Cmd+Option+R` (Mac)

### Option 2: R Console
```bash
cd bulk-rnaseq-differential-expression-r/analysis
R
```

Then in R:
```r
source("main.R")
```

### Option 3: Command Line
```bash
cd bulk-rnaseq-differential-expression-r/analysis
Rscript main.R
```

**Important**: Always run from within the `analysis/` directory to ensure proper working directory context.

## Data & Inputs

This pipeline uses the **airway** dataset from Bioconductor, which contains RNA-seq data from airway smooth muscle cells treated with dexamethasone.

**Dataset details:**
- 8 samples (4 treated, 4 untreated)
- ~63,000 genes
- Dexamethasone treatment vs control

**For custom data**, you would need:
- Count matrix (genes × samples)
- Sample metadata with experimental conditions
- Gene annotations (or use built-in annotation databases)

## Outputs

The pipeline generates the following outputs (displayed in R session, can be exported):

### Statistical Results
- **Differential expression tables** at α = 0.05 and α = 0.001
- Top 10 differentially expressed genes with annotations
- Summary statistics (log2 fold changes, adjusted p-values)

### Visualizations
1. **MA plots** - Log fold change vs mean expression (3 versions)
2. **Boxplots** - Expression distribution across samples
3. **Barplots** - Mean expression levels
4. **Heatmaps** - Top 10 DEGs with hierarchical clustering (2 versions)
5. **Scatter plots** - Correlation between selected genes

### Enrichment Analysis
- Pathway enrichment results from 9 databases:
  - KEGG_2021_Human
  - Reactome_2022
  - GO_Biological_Process_2023
  - GO_Molecular_Function_2023
  - GO_Cellular_Component_2023
  - OMIM_Disease
  - OMIM_Expanded
  - DisGeNET
  - Human_Phenotype_Ontology

### Annotations
- Gene symbol mapping for top differentially expressed genes
- Correlation statistics for gene pairs

## Reproducibility Notes

- The analysis uses the `airway` dataset which is version-controlled through Bioconductor
- DESeq2 performs internal normalization (no manual scaling required)
- For custom analyses, consider setting a random seed: `set.seed(123)`
- R session info can be captured with `sessionInfo()`
- Bioconductor version: 3.18 recommended

## Methodology Details

### Differential Expression
1. **Pre-filtering**: Remove genes with counts < 3 across all samples
2. **Normalization**: DESeq2 size factor normalization
3. **Testing**: Negative binomial generalized linear model
4. **Multiple testing correction**: Benjamini-Hochberg FDR

### Pathway Enrichment
- Uses enrichR to query 9 biological databases
- Fisher's exact test for over-representation
- Results include adjusted p-values, odds ratios, and gene lists

### Visualization
- rlog transformation for variance-stabilized visualization
- Hierarchical clustering with complete linkage
- Pearson correlation for gene-gene relationships

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Package installation fails | Update BiocManager: `BiocManager::install(version = "3.18")` |
| "Cannot open connection" error | Ensure working directory is `analysis/` |
| Memory issues | Increase memory limit: `options(future.globals.maxSize = 8000 * 1024^2)` |
| Plots not displaying | Check graphics device: `dev.cur()` |
| enrichR connection fails | Check internet connection; enrichR requires online access |
| DESeq2 convergence warnings | Normal for some genes; check specific gene results |

## Performance Considerations

- **Runtime**: ~2-5 minutes on standard hardware (airway dataset)
- **Memory**: ~2-4 GB RAM required
- **Internet**: Required for enrichR pathway queries

## Extension Ideas

- Apply to custom RNA-seq datasets
- Add volcano plots for visualization
- Implement time-series or multi-factor designs
- Export results to CSV/Excel for sharing
- Add PCA visualization for sample clustering
- Integrate with gene set enrichment analysis (GSEA)

## License

This project is provided as-is for educational and research purposes.

## Contact & Contributions

Issues and suggestions are welcome. Please open an issue for any bugs or feature requests.
