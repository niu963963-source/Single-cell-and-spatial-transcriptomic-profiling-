Single-cell and spatial transcriptomic analysis of tertiary lymphoid structure (TLS) maturation in gastric cancer

This repository contains the R analysis code for the study:


Single-cell and spatial transcriptomic profiling reveals impaired tertiary lymphoid structure maturation associated with T cell exhaustion in gastric cancer
Longteng Niu, Yan Li



The code reproduces the integrated single-cell RNA-seq, 10x Visium spatial transcriptomics, cell–cell communication, T-cell trajectory, three-tier spatial TLS maturity scoring, and TCGA-STAD survival analyses reported in the manuscript.


Data availability

The analysis uses publicly available datasets (not redistributed here):

DatasetTypeSource / AccessionGSE206785Single-cell RNA-seq (gastric cancer)GEO: https://identifiers.org/geo:GSE206785GSE25195010x Visium spatial transcriptomics + H&EGEO: https://identifiers.org/geo:GSE251950TCGA-STADBulk RNA-seq (TPM) + survivalUCSC Xena: https://xenabrowser.net

Download these into the working directory before running. Expected input files include:
, ,
 (SpaceRanger outputs per sample),
, , .GSE206785_scgex.txt.gzGSE206785_metadata.txt.gzGSE251950_RAW/TCGA-STAD.star_tpm.tsvTCGA-STAD.survival.tsvgencode.v36.annotation.gtf.gene.probemap


Repository contents

FileDescriptiongastric_TLS_analysis.RMain analysis pipeline (modules 0–18)TLS_maturity_score.RThree-tier spatial TLS maturity scoring (None/Early/Mature) with adaptive neighbourhood radiusTLS_maturity_sensitivity.R144-combination sensitivity analysis of the maturity scoreREADME.mdThis file


If you keep everything in a single script, name it  and ignore the split above.gastric_TLS_analysis.R




Analysis modules (main pipeline)

ModuleStep0Environment setup, package installation, plotting theme1Data loading (batched read of scRNA-seq matrix)2Quality control (nFeature, mitochondrial %)3Doublet removal (MAD-based score)4Normalization, variable features, PCA5Harmony batch correction6–7UMAP and clustering (resolution 0.2)8Cell-type annotation (SingleR + canonical markers), unified naming9–10Tumor-vs-normal DEG; subtype composition11CellChat intercellular communication12TLS module score, immune-checkpoint score, χ² enrichment / Cramér's V13Spatial integration: label transfer, CXCL13–CXCR5 Ripley's K-cross, CXCL12–CXCR414–15TCGA-STAD survival (TLS / exhaustion / CAF signatures), ROC16NicheNet ligand–target analysis17T-cell sub-clustering and Slingshot pseudotime18Exhausted T-cell visualization+Three-tier spatial TLS maturity scoring and sensitivity analysis


Requirements


R ≥ 4.3
Key packages:  (v5), , , ,  (v2),  (v2), , , , , , , , , , SeuratharmonySingleRcelldexCellChatnichenetrslingshotSingleCellExperimentspatstatsurvivalsurvminerpROCggplot2data.tablepheatmapComplexHeatmap


Install Bioconductor packages via  and CRAN packages via  (handled automatically in Module 0).BiocManager::install()install.packages()


Usage


Place the input data files (see Data availability) in the working directory.
Update the  path at the top of the script to your local directory.setwd()
Run the script module by module in R / RStudio. Intermediate objects are saved under  and figures under .RData/Figures/



Citation

If you use this code, please cite the associated manuscript (citation to be updated upon publication).

License

Released under the MIT License (see ).LICENSE
