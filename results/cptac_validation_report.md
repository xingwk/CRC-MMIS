# CPTAC COAD External Validation Report

> Generated: 2026-06-07 | Data source: cBioPortal coad_cptac_2019 (Vasaikar et al. Cell 2019)

## Summary Table

| Metric | CPTAC COAD | TCGA Training | GEO Meta (5 cohorts) |
|--------|-----------|--------------|---------------------|
| Sample size | **106** | 431 | ~1,122 |
| Platform | RNA-seq (RSEM UQ Log2) | RNA-seq (RSEM) | Microarray |
| Median risk score | 13.29 | — | — |
| Risk score range | 11.26–15.10 | — | — |
| OS data | ❌ Not available | ✅ | ✅ |
| MSI status | ✅ (MSS/MSI-H) | ✅ | ❌ |

## Data Downloaded

### Files (in `data/cptac/cbioportal/`)
- **`cptac_expression.csv`**: RSEM UQ Log2 expression for 10 signature genes × 106 samples
- **`cptac_risk_scores.csv`**: Risk scores computed using manuscript coefficients
- **`cptac_clinical.csv`**: Available clinical data (MSI status, cancer type)
- **`cptac_merged.csv`**: Combined expression + clinical + risk scores

### Expression Coverage
- 106/110 samples have all 10 genes expressed (4 samples missing some → excluded)
- All 10 signature genes detected (RSEM UQ Log2 normalized)

### Clinical Data Available
- ✅ **MSI_STATUS**: MSI-H and MSS labels (105 samples)
- ✅ **CANCER_TYPE**: All colon adenocarcinoma
- ❌ **OS_MONTHS / OS_STATUS**: Not available via public API
- ❌ **AGE / SEX / STAGE**: Not served by cBioPortal API for this study

## Limitation

CPTAC COAD survival data is **not accessible** through the public cBioPortal API or the current GDC CPTAC-3 download. The original study (Vasaikar et al. 2019, Cell) survival data exists but is behind Cell Press supplementary materials.

For manuscript purposes, this dataset can supplement the existing 5 GEO cohorts as:
1. **Supplementary Table**: Expression of 10-gene signature in an independent RNA-seq cohort
2. **Descriptive**: Risk score distribution in CPTAC COAD
3. **MSI correlation**: Compare risk scores between MSS vs MSI-H

## Alternative Sources for RNA-seq Validation (with survival)

| Dataset | n | Platform | Survival Data | Accessibility |
|---------|---|----------|--------------|---------------|
| GSE50760 | 18 (54 samples) | RNA-seq | ❌ No | GEO, but small |
| TCGA-COAD/READ | ~592 | RNA-seq | ✅ Yes | Already training set |
| GSE39582 | 443 | Microarray | ✅ Yes | Already in meta |
| GSE17536 | 177 | Microarray | ✅ Yes | Already in meta |
| GSE14333 | 290 | Microarray | ✅ Yes | Already in meta |
| GSE38832 | 122 | Microarray | ✅ Yes | Already in meta |
| GSE33113 | 90 | Microarray | ✅ Yes | Already in meta |
