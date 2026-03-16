# Editas CRISPR Analysis Pipeline  

This repository contains Python scripts and workflows I developed during my Co-op at **Editas Medicine** to analyze CRISPR-induced transcript editing outcomes. The codebase is organized into modular folders for clarity.  

## Repository Structure  
- **`pipeline/`** → Core scripts for RNA-seq processing and CRISPR editing quantification  
  - Alignment to cDNA, reference generation, FASTQ-to-BAM conversion, transcript analysis  
- **`visualizations/`** → Data visualization scripts  
  - Volcano plots, DESeq2 plots, CPM analysis, GSEA results, and other matplotlib-based visualizations  
- **`utils/`** → Supporting utilities and automation scripts  
  - Data handling, JIRA integration, average editing calculations, helper functions  

##  Key Highlights  
- **Reusable Python pipeline** for quantifying transcript-level editing outcomes, optimized for speed and scalability  
- **Bulk RNA-seq analysis of 18 samples from AWS S3** using **Snakemake**, covering QC → alignment → downstream analysis  
- **Cloud automation** with **AWS Batch** for large-scale, reproducible execution  
- **Rich visualizations** (heatmaps, volcano plots, GSEA) to compare and interpret guide RNA performance  

---
