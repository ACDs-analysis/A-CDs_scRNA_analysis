# A-CDs_scRNA_analysis

Custom scripts for single-cell RNA sequencing (scRNA-seq) downstream analysis in the A-CDs sepsis study.

## Overview

This repository contains custom R scripts used to analyze single-cell transcriptomic data, focusing on:

- Macrophage–endothelial cell communication (CellChat)
- Endothelial barrier integrity scoring
- Endothelial leakage scoring
- Machine learning-based analysis

## Scripts

- `Cellchat_LR_Bitmap_EC_to_MPhi.R`  
  Ligand–receptor interaction analysis from endothelial cells to macrophages.

- `Cellchat_LR_Bitmap_MPhi_to_EC.R`  
  Ligand–receptor interaction analysis from macrophages to endothelial cells.

- `EC_Barrier_score.R`  
  Calculation of endothelial barrier integrity scores.

- `EC_Leak_score.R`  
  Quantification of endothelial leakage-related signatures.

- `ML_analysis.R`  
  Machine learning workflow for feature selection and predictive modeling.

## Requirements

- R (>= 4.2)
- Seurat
- CellChat
- tidyverse
- ggplot2

## Notes

- Scripts are provided for reproducibility of the analysis presented in the study.
- Input data and intermediate files are not included.

## Author

Panpan Yi
