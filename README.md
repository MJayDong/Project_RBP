# Project_RBP
This repository contains a fully modularized single-cell transcriptomics analysis framework designed for large-scale multi-sample datasets (e.g., normal → primary → metastasis).

The pipeline supports:
	•	Full Seurat preprocessing
	•	Cell type–specific sub-analysis
	•	RBP-focused regulatory inference
	•	SCENIC / scMLnet / MultiNicheNet network analysis
	•	Cross-cell-type integration
	•	Statistical downstream tests
	•	High-quality figure generation for publication

  00_core — Core Utilities

Contains fundamental helper functions used across the pipeline.
	•	0.BASIC.R
Unified utilities: file loading/saving, message printing, path control.
	•	0.SeuratProject.R
Project initialization and structured directory creation.

⸻

01_preprocessing — Global QC & Integration

Scripts for preprocessing the complete dataset.
	•	0.All_Preprocessing.R
QC filtering, normalization, HVG detection, PCA, Harmony integration.
	•	1.All_Marker_DE.R
Marker discovery & global differential expression across conditions.
	•	2.All_Visualization.R
Global UMAP/tSNE visualization colored by various metadata.

⸻

02_celltype_modules — Per-Cell-Type Analysis

Independent modules that deeply characterize each major cell lineage.
	•	CellType_Basic.R
Subsetting by BigGroup, clustering, within-celltype DE, basic TSNE.
	•	CellType_RBP.R
RBP gene intersection with HVGs and feature sets.
	•	CellType_Relation.R
TF–Target–RBP intersection based on scMLnet+SCENIC outputs.

⸻

03_crosscell_analysis — Cross-Lineage Integration

Integrative analysis of all cell types.
	•	KeyRBP_Identification.R
Multi-feature ranking of key RBPs across the dataset.
	•	Cluster_GVis.R
Cluster-level functional enrichment via ClusterGVis.
	•	GeneSignificance.R
Statistical testing of gene expression across clusters or conditions.

⸻

04_network_inference — Regulatory & Communication Networks

These modules implement major network-biology frameworks.
	•	Ratio_OR_TissuePref.R
Cell-state ratio calculation, enrichment (OR), tissue preference scoring.
	•	SCENIC_Run.R
Full SCENIC workflow to infer TF regulons (GENIE3 → RcisTarget → AUCell).
	•	scMLnet_Run.R
Multi-layer ligand → receptor → TF → target signaling analysis.
	•	MultiNicheNet_Run.R
Multi-celltype ligand activity and ligand–target predictions.

All scripts output .rds and .csv files for downstream interpretation.

⸻

05_figures — Publication-Ready Figure Scripts

High-resolution figure generation for main and supplementary figures.
	•	Fig1A_Tsne_UMAP.R
Clean black-and-white or color tSNE/UMAP plots.
	•	Fig1D_VennPlot.R
Venn plots for gene-set intersections.
	•	SupplyFig2_CustomPlots.R
DotPlots, violin plots, heatmaps for supplementary panels.
