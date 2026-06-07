# CRC-MMIS: Colorectal Cancer Malignant-Myeloid Interaction Signature

## Repository Structure

```
├── scripts/          # R 4.6.0 and Python 3.13 analysis scripts
├── data/             # Processed expression matrices and clinical data
├── results/          # Output figures, tables, and analysis reports
├── submission/       # Manuscript, cover letter, supplementary materials
└── README.md
```

## Requirements

- **R 4.6.0** with packages: survival, glmnet, survminer, WGCNA, car, ggplot2, patchwork, pheatmap, readr, dplyr
- **Python 3.13** with: scanpy, scipy, numpy, pandas, matplotlib, h5py, urllib3

## Reproducing the Analysis

1. Clone this repository
2. Run `scripts/01_data_preprocessing.R` → `scripts/02_wgcna.R` → `scripts/03_lasso_signature.R`
3. External validation: `scripts/04_external_validation.R`
4. Generate figures: `scripts/generate_figures_R.R`

## Data Sources

- scRNA-seq: GSE146771, GSE132465, GSE205506 (GEO)
- Bulk RNA-seq training: TCGA-COADREAD (GDC Portal)
- External validation: GSE39582, GSE17536, GSE14333, GSE38832, GSE33113 (GEO)
- RNA-seq validation: CPTAC COAD (LinkedOmics)

## Zenodo DOI

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20582855.svg)](https://doi.org/10.5281/zenodo.20582855)
