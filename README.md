# A Parallel Latent-Space Framework for Interpretable Multi-Omics Integration

This repository contains the code, results, and figures for the paper: 
**"A parallel latent-space framework for interpretable multi-omics integration in microbiome studies"**

## Overview
We propose a novel parallel autoencoder architecture (PLAE) that integrates metagenomic and functional metagenomic data without losing omics-specific features. It uses Maximum Mean Discrepancy (MMD) to align latent spaces and provides interpretable biomarker discovery for Inflammatory Bowel Disease (IBD).

## Requirements
- R (>= 4.3.0)
- R packages: `keras3`, `reticulate`, `pROC`, `curatedMetagenomicData`, `ggplot2`, `cowplot`, `dplyr`, `SummarizedExperiment`
- Python (>= 3.10) with `tensorflow`, `keras`, `numpy`, `pandas`

## Data
We use the iHMP (Integrative Human Microbiome Project) IBD dataset, accessed directly via the `curatedMetagenomicData` Bioconductor package (1,627 samples). 

## Methodology Highlights
- **No Data Leakage:** Feature selection (Wilcoxon rank-sum test, top 200 species & top 100 pathways) and Min-Max scaling are performed strictly *inside* each cross-validation fold using only training data.
- **Strict Patient-Level 5-Fold Cross-Validation:** Prevents the model from memorizing patient-specific microbiome signatures.
- **Interpretability:** First-layer weight magnitudes are used to identify key IBD-associated taxa.

## Results
- **Patient-level 5-Fold CV Mean AUC:** **0.646 ± 0.081**
- **Top Biomarkers identified:** *Collinsella intestinalis*, *Gemmiger formicilis*, *Veillonella spp.*

## How to Run
1. Clone this repository.
2. Open `main_analysis.R` in RStudio.
3. Ensure the Python environment is configured (`Sys.setenv(KERAS_BACKEND = "tensorflow")`).
4. Run the script from top to bottom. All generated figures and CSV tables will be saved in the `results/real_data/` folder.

## Code Availability
All code is available in this repository. For questions, contact shegoftejalal@gmail.com.