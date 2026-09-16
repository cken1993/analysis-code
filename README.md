# Analysis code for TLS-associated chemokine programs in HCC

This repository contains the analysis code for the manuscript:

**Tertiary lymphoid structure-associated chemokine programs define an immune-organized tumor microenvironment with prognostic relevance in hepatocellular carcinoma**

*World Journal of Gastrointestinal Oncology*, Manuscript No. 126440

The code is provided to reproduce the main analysis results presented in the article, including bulk transcriptomic analysis, single-cell RNA sequencing analysis, immune infiltration analysis, and construction/validation of the TLS-related prognostic model.

## Files

| File | Description |
|------|-------------|
| `TLS state class.R` | Defines concordant TLS-positive and TLS-negative samples from histopathological annotations in TCGA-LIHC. |
| `differrntial analysis.R` | Differential expression analysis between TLS-positive and TLS-negative tumors using DESeq2. |
| `geo_validation.ipynb` | Processes GEO datasets (GSE54236, GSE64041, GSE76427), calculates the 12-chemokine TLS signature score, and performs limma differential expression analysis. |
| `immune inflitration.R` | Immune infiltration analysis (ssGSEA/GSVA) and correlation with the TLS signature. |
| `lasso_cox.ipynb` | LASSO-Cox regression to build the seven-gene prognostic model and external validation in ICGC-LIRI-JP. |
| `bootstrap c-index.R` | Bootstrap stability analysis of LASSO-selected genes, Kaplan–Meier curves, C-index, and time-dependent ROC. |
| `scRNA annotation.R` | Seurat processing and cell-type annotation of scRNA-seq data (GSE149614). |
| `scRNA cell communication.R` | CellChat analysis of CCL/CXCL signaling networks and CCL19–CCR7 interactions. |

## Requirements

- R (>= 4.0.4)
- Main R packages: Seurat, SingleR, DESeq2, limma, GSVA, clusterProfiler, glmnet, survival, timeROC, CellChat, IOBR, ComplexHeatmap.

## Data

All datasets analyzed are publicly available. TCGA-LIHC data were obtained from the GDC Data Portal. GEO datasets (GSE54236, GSE64041, GSE149614) were obtained from NCBI GEO. ICGC-LIRI-JP controlled-access data are available through EGA (EGAD00001003547). Histopathological TLS annotations were derived from previously published HCC TLS studies. Please refer to the manuscript for full details.

