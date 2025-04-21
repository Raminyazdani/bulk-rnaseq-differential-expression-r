# Project Identity

## Display Title
Bulk RNA-seq Differential Expression Analysis

## Repo Slug
bulk-rnaseq-differential-expression-r

## Tagline
Automated differential expression analysis pipeline for bulk RNA-seq data using DESeq2 and pathway enrichment

## GitHub Description
A reproducible R-based pipeline for bulk RNA-seq differential expression analysis, featuring DESeq2 normalization, statistical testing, pathway enrichment, and comprehensive visualization of gene expression patterns.

## Primary Stack
- R (4.x)
- Bioconductor
- DESeq2
- enrichR

## Topics/Keywords
- bioinformatics
- rna-seq
- differential-expression
- deseq2
- bioconductor
- genomics
- transcriptomics
- pathway-enrichment
- data-analysis
- r-statistics

## Problem & Approach
**Problem**: Analyzing bulk RNA-sequencing data to identify differentially expressed genes between experimental conditions requires normalization, statistical testing, and biological interpretation through pathway analysis.

**Approach**: This pipeline uses DESeq2 for robust differential expression analysis with negative binomial modeling, followed by pathway enrichment using multiple databases (KEGG, Reactome, GO) via enrichR. Includes quality control visualizations, MA plots, heatmaps, and correlation analysis.

## Inputs & Outputs

**Inputs**:
- Count matrix (genes × samples) or SummarizedExperiment object
- Sample metadata with experimental conditions
- Example uses the airway dataset from Bioconductor

**Outputs**:
- Differential expression results (log2 fold changes, adjusted p-values)
- Normalized count matrices
- MA plots showing differential expression
- Heatmaps of top differentially expressed genes
- Pathway enrichment results from multiple databases
- QC plots (boxplots, expression distributions)
- Correlation analysis of selected genes
