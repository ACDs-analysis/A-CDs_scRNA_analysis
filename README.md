# A-CDs_scRNA_analysis

Custom R scripts for downstream single-cell RNA sequencing (scRNA-seq) analysis in the A-CDs sepsis study.

## Overview

This repository contains custom scripts used to analyze single-cell transcriptomic data from the A-CDs sepsis project, with a primary focus on macrophage–endothelial interactions, endothelial functional states, macrophage polarization, pathway enrichment, and machine learning-assisted prioritization.

The analyses in this repository include:

- Ligand–receptor communication analysis between macrophages and endothelial cells using CellChat
- Visualization of CellChat-derived gene and transcription factor expression patterns
- Information flow analysis of intercellular signaling pathways
- Endothelial barrier integrity and leakage signature scoring
- Macrophage M1/M2 polarization scoring and visualization
- Pathway enrichment analysis, including TNF signaling and ECM-related programs
- Machine learning-based analysis for feature prioritization and predictive modeling

## Scripts

### CellChat analysis

- `Cellchat_LR_Bitmap_EC_to_MPhi.R`  
  Ligand–receptor interaction analysis from endothelial cells to macrophages.

- `Cellchat_LR_Bitmap_MPhi_to_EC.R`  
  Ligand–receptor interaction analysis from macrophages to endothelial cells.

- `cellchat_ec_to_macrophage_information_flow.R`  
  Information flow analysis of CellChat signaling pathways from endothelial cells to macrophages.

- `cellchat_macrophage_to_ec_information_flow.R`  
  Information flow analysis of CellChat signaling pathways from macrophages to endothelial cells.

- `cellchat_endothelial_gene_expression.R`  
  Visualization of CellChat-related gene expression in endothelial cells.

- `cellchat_macrophage_gene_expression.R`  
  Visualization of CellChat-related gene expression in macrophages.

- `cellchat_macrophage_tf_expression.R`  
  Analysis and visualization of transcription factor expression patterns in macrophages.

### Endothelial functional scoring

- `EC_Barrier_score.R`  
  Calculation of endothelial barrier integrity scores.

- `EC_Leak_score.R`  
  Quantification of endothelial leakage-related signatures.

### Macrophage polarization analysis

- `m1_m2_polarization_heatmap.R`  
  Heatmap visualization of M1/M2 polarization-associated gene expression patterns.

- `m1_m2_score_boxplot_by_treatment.R`  
  Boxplot visualization of M1/M2 scores grouped by treatment condition.

- `m1_m2_score_boxplot_by_treatment_and_cluster.R`  
  Boxplot visualization of M1/M2 scores grouped by both treatment condition and cell cluster.

### Pathway enrichment analysis

- `tnf_pathway_gsea.R`  
  Gene set enrichment analysis of the TNF signaling pathway.

- `ecm_pathway_gsea.R`  
  Gene set enrichment analysis of ECM-related pathways.

### Machine learning analysis

- `ML_analysis.R`  
  Machine learning workflow for feature selection, prioritization, and predictive modeling.

## Requirements

Recommended R environment and commonly used packages include:

- R (>= 4.2)
- Seurat
- CellChat
- monocle / monocle2
- clusterProfiler
- fgsea
- tidyverse
- dplyr
- ggplot2
- pheatmap or ComplexHeatmap
- patchwork

Please note that the exact package requirements may vary across scripts.

## Notes

- These scripts were developed for reproducibility of the downstream analyses performed in the A-CDs sepsis study.
- Raw sequencing data, processed matrices, and intermediate files are not included in this repository.
- Some scripts may require user-specified input paths, metadata tables, marker gene sets, or CellChat objects before execution.
- File names reflect specific downstream analysis tasks and are intended to facilitate script reuse and adaptation.

## Author

Panpan Yi
