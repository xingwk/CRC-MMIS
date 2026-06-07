# CRC-MMIS Submission Data Package

## Contents

### Manuscript
- CRC-MMIS-Manuscript-v3.5.docx - Main manuscript (PeerJ format, Times New Roman 12pt)

### Cover Letter & Highlights
- Cover_Letter.md
- Highlights.md

### figures/
Publication-quality figures (PDF + 300dpi JPG):
- Figure1_scRNA - Single-cell RNA-seq profiling
- Figure2_Validation - Multi-cohort external validation (incl. CPTAC COAD)
- Figure3_CMS - CMS subtype stratification
- Figure4_Adenosine - ENT1 and adenosinergic signaling
- Figure5_TIDE - TIDE immunotherapy prediction
- Figure6_MSI_Survival - MSI-stratified survival
- FigureS10_CellChat - CellChat ligand-receptor interaction analysis (GSE146771)
- FigureS11_GSE205506 - CellChat analysis of PD-1-treated cohort
- FigureS12_CPTAC_Validation - CPTAC COAD RNA-seq external validation

### tables/
- CRC-MMIS_Data_Tables.xlsx - Multi-sheet Excel with Tables 1-3 + supplementary data

### data/
Individual CSV data files:
- risk_scores_seed42.csv / seed123.csv - Risk scores per patient
- survival_analysis.csv - Cox regression results (R survival::coxph)
- km_plot_data.csv - KM curve source data
- tide_stratification.csv - MSI-Risk TIDE response rates
- risk_tide_correlation.csv - Risk-TIDE Spearman correlations
- gse205506_celltype_expression.csv - scRNA-seq gene expression per cell type
- forest_plot_meta.csv - Forest plot meta-analysis data (incl. CPTAC COAD)
- cptac_validation.csv - CPTAC COAD validation data (105 samples, 10 genes, risk scores)
- cellchat_interactions.csv - CellChat ligand-receptor interaction results (86 interactions)
- tcga_clinical.csv - TCGA clinical annotations
- tcga_cms_labels.csv - CMS subtype labels

## Contact
Corresponding author: [to be completed]
