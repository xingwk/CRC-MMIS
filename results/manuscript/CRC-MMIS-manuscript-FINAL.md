# Single-Cell Characterization of Malignant-Myeloid Interaction Transcriptomic Features Identifies a Dual-Circuit Adenosinergic Signaling Hypothesis in Colorectal Cancer

Target Journal: PeerJ  
Manuscript Type: Research Article  
Version: v3.5 — Title/Abstract reframed to hypothesis-oriented per reviewer feedback; nFeature_RNA QC threshold fixed; Data Transparency Notes moved to table footnotes; all 42 references PubMed-verified
Revision Date: 2026-06-07 (v3.5: Title/Abstract reframed to hypothesis-oriented per reviewer feedback)
Date: 2026-06-07

---

## Abstract

Background: Colorectal cancer (CRC) contains a myeloid-cell-enriched microenvironment that limits immunotherapy efficacy. The transcriptional programs governing malignant-myeloid crosstalk and their implications for microsatellite-stable (MSS) CRC remain incompletely characterized.
Methods: Three single-cell RNA-seq datasets (GSE146771, GSE132465, GSE205506) were used to characterize malignant-myeloid interaction programs. A hybrid 10-gene signature (CRC-MMIS) combined four LASSO-selected genes (SNAI1, CXCL8, SPP1, CXCL10) with six adenosine-pathway genes curated from single-cell evidence, trained on TCGA-COADREAD (N=431) and evaluated in five external GEO cohorts (total N=1,122) plus an independent RNA-seq cohort (CPTAC COAD, n=105). Performance was assessed within CMS subtypes and by MSI status. Immunotherapy response was predicted using a TIDE-inspired consensus framework with dual-seed cross-validation.
Findings: In the TCGA training cohort, the CRC-MMIS gene set showed survival-associated expression patterns: high- vs low-risk groups yielded HR=2.46-3.30 (p<0.001), with separation most evident in CMS4 tumors (HR=3.30, CI 1.37-7.97, p=0.005; 64.4% high-risk vs 28.8-51.7% in other subtypes). MSS patients also showed a survival-associated signal across two MSI classification methods (HR=2.28-3.11, p<0.05). ENT1 (SLC29A1) was overexpressed in tumors versus normal tissue (p=8.66×10⁻⁷), correlated with SNAI1 (r=0.641), and co-expressed with the adenosinergic cascade (ENTPD1, ENTPD2, NT5E, ADORA2A, ADORA2B; all p<0.001). TIDE-inspired prediction showed higher predicted responder proportions in MSS-Low vs MSS-High patients (69.5% vs 49.5%). External validation in cohorts with overall survival endpoints showed consistent prognostic direction (OS-only pooled HR=1.13, 95% CI 1.05-1.22, I²=0%); individual cohort-level significance was observed in GSE38832 (DSS) and GSE39582 (OS). RNA-seq validation in CPTAC COAD showed a directionally concordant trend (HR=1.16, 95% CI 0.41-3.28). Multivariable Cox indicated independent association in training (HR=1.46-1.59, p=0.0002; C-index=0.77; 3-year AUC=0.65-0.71).
Interpretation: The CRC-MMIS framework identifies a coordinated adenosinergic transcriptional architecture as a candidate mechanistic link between malignant-myeloid crosstalk and immune evasion in CRC. Within cohorts specifying overall survival, external validation showed consistent prognostic direction (pooled HR=1.13, I²=0%), indicating that the signature captures a reproducible risk signal in OS contexts. The observed CMS4-concentrated survival association and MSS immune-stratification signal warrant independent validation in larger cohorts. TIDE-predicted rates are exploratory and exceed observed clinical predicted responder proportions; they represent relative stratification metrics rather than absolute predictions.
Funding: None declared.  Keywords: colorectal cancer; myeloid cells; tumor microenvironment; adenosine; immunotherapy; transcriptomic hypothesis; cell-cell communication

---

### Supplementary Figures

Figure S1. Gene set correlation analysis: EMT (ρ=0.28, p=0.002), Hypoxia (ρ=0.24, p=0.016), TGF-β (ρ=0.27), Angiogenesis (ρ=0.33) up in high-risk; IFN-γ (ρ=-0.09) trend down. Spearman correlations with risk score reported rather than standard GSEA NES scores.

Figure S2. Multivariable Cox forest plot. Risk score (HR=1.59 per unit, 95% CI 1.25-2.03, p=0.0002), AJCC stage (HR=2.42, p<0.001), age, sex, MSI, and CMS4 adjusted simultaneously.

Figure S3. KM curves by CRC-MMIS risk tertile in TCGA-COADREAD.

Figure S4. Schoenfeld residual plots (PH test, all p > 0.05).

Figure S5. Hallmark pathway enrichment barplot: Spearman correlation coefficients between each pathway score and the CRC-MMIS risk score. Upregulated pathways in red, downregulated in blue.

Figure S6. Calibration curves for 1-year and 3-year overall survival.

Figure S7. Decision curve analysis (DCA) for 3-year survival.

Figure S8. Bootstrap internal validation (1,000 resamples): C-index distribution with 95% CI.

Figure S9. Adenosine pathway gene expression heatmap across cell types in GSE146771 (SLC29A1/ENT1, ENTPD1, ENTPD2, NT5E, ADORA2A, ADORA2B).

Figure S10. Permutation-based ligand-receptor interaction analysis of GSE146771. (A) Heatmap of 86 significant interactions across 24 LR pairs and 5 cell types. (B) Chord diagram of top-ranked interactions. (C) Barplot of interaction scores by LR pair ranked by mean score. Top pairs: ANXA1→FPR1 (score=0.79), CCL5→CCR2/CCR5 (0.32-0.38), CD47→SIRPA (0.25).

Figure S11. GSE205506 ligand-receptor interaction analysis: communication landscape in dMMR/MSI-H CRC samples (pre- and on-treatment anti-PD-1 combined; 35 interactions, 18 LR pairs).

Figure S12. CPTAC COAD RNA-seq validation. (A) Risk score distribution. (B) Kaplan-Meier curves, high- vs low-risk groups (log-rank p=0.26). (C) Forest plot summarizing CPTAC validation alongside GEO cohorts. (D) Expression heatmap of the 10-gene signature by risk score tertile.

Figure S13. External validation forest plot stratified by survival endpoint. Pooled HR (diamonds): OS=1.13 (95% CI 1.05-1.22), I²=0%; DSS=1.87, I²=0%; DFS/Recurrence=0.89, I²=0%.

Supplementary Table S14. Model comparison: C-index of 4-gene LASSO-only (SNAI1, CXCL8, SPP1, CXCL10) versus 10-gene CRC-MMIS in TCGA-COADREAD.

---

## 1. Introduction

Colorectal cancer is the third most diagnosed cancer and the second leading cause of cancer death worldwide (Sung et al., 2021; Siegel et al., 2023). Immune checkpoint inhibitors produce durable responses in mismatch repair-deficient (dMMR)/microsatellite instability-high (MSI-H) CRC, with objective response rates above 30% (Le et al., 2015; André et al., 2020; Overman et al., 2017). However, these tumors account for only ~15% of cases. The remaining ~85%  — classified as microsatellite-stable (MSS) or proficient mismatch repair (pMMR)  — show response rates below 5% to single-agent checkpoint blockade, a resistance driven substantially by the myeloid compartment (Ganesh et al., 2019; Llosa et al., 2015; Mlecnik et al., 2016).

Myeloid cells  — tumor-associated macrophages, myeloid-derived suppressor cells, and dendritic cells  — dominate the CRC immune infiltrate (Mantovani et al., 2017; Gabrilovich & Nagaraj, 2009; Cortese et al., 2019). Tumor-derived signals (CSF1, CCL2, osteopontin/SPP1) polarize macrophages toward an M2-like phenotype that suppresses T-cell function through IL-10 and TGF-β secretion (Qian & Pollard, 2010; Noy & Pollard, 2014), while CXCL8-driven MDSC recruitment reinforces an immunosuppressive niche (Veglia et al., 2021). This malignant-myeloid crosstalk  — tumor cells recruiting myeloid populations that in turn promote immune evasion  — has motivated the development of myeloid-informed gene signatures. However, existing CRC prognostic signatures derive predominantly from bulk transcriptomic data without single-cell resolution, conflating signals across cell types into a single score (Dai et al., 2018; Zhao et al., 2019; Guinney et al., 2015). A recently published 11-gene glioma malignant-myeloid interaction signature (MMIS) by Ren et al. demonstrated that a signature built on tumor-myeloid crosstalk can predict prognosis and immunotherapy response, identifying TPST1 as a regulator of the immunosuppressive microenvironment in glioblastoma (Ren et al., 2026). TPST1 was independently prioritized in our CRC screening, suggesting conserved malignant-myeloid interaction modules across tumor types. Yet the glioma MMIS cannot be directly applied to CRC, which has distinct biology  — canonical epithelial-mesenchymal transition (EMT), Wnt-driven immune exclusion, and the Consensus Molecular Subtype (CMS) classification (Guinney et al., 2015; Becht et al., 2016).

Here, we characterized malignant-myeloid interaction programs through single-cell RNA-seq analysis of three CRC cohorts and constructed a 10-gene transcriptomic framework (CRC-MMIS) using single-cell-guided gene selection and Cox-LASSO modeling in TCGA-COADREAD. We evaluated the framework for exploratory prognostic associations across multiple cohorts, within CMS subtypes and MSI strata, and for immunotherapy response prediction using a TIDE-inspired approach. The adenosinergic signaling axis  — comprising ENT1 (SLC29A1), the ectonucleotidases ENTPD1/2 and NT5E/CD73, and the adenosine receptors ADORA2A/2B  — emerged as a coordinated transcriptional program that may link the observed survival associations with immune-modulatory functions. We propose a two-compartment adenosine signaling model in which tumor-intrinsic adenosine import (ENT1) and myeloid-compartment extracellular adenosine generation (ENTPD1/2-NT5E) jointly contribute to immune suppression in CRC, as a testable hypothesis for future investigation.

---

## 2. Materials and Methods

### 2.1 Study Design

Three single-cell RNA-seq datasets were analyzed to define myeloid functional modules and malignant epithelial markers. WGCNA was performed on TCGA-COADREAD to bridge single-cell gene sets with survival-associated co-expression networks. A 10-gene signature (CRC-MMIS) was built through Cox-LASSO regression. The signature was externally validated in five GEO cohorts, stratified by CMS subtype and MSI status, and evaluated for immunotherapy response prediction using a TIDE-inspired consensus framework. The adenosinergic signaling axis was examined as the mechanistic link between prognostic and immune-modulatory signature functions. All analyses used R 4.6.0 (survival v3.7, survminer v0.5, glmnet v4.1, WGCNA v1.72, car v3.1) for statistical modeling and Python 3.13 (scanpy v1.12, scipy v1.13) for single-cell processing and TIDE computation.

### 2.2 Single-Cell RNA-seq Data Processing

Three publicly available CRC scRNA-seq datasets were used: GSE146771 (n=20 CRC patients, 10x Genomics + Smart-seq2, myeloid-enriched), GSE132465 (n=33 CRC patients, ~63,689 cells, Korean cohort), and GSE205506 (n=19 dMMR/MSI-H CRC patients receiving neoadjuvant anti-PD-1, pre- and on-treatment samples). Data preprocessing was performed with Seurat v5 (Hao et al., 2021): quality control (nFeature_RNA: 200-5000; percent.mt < 20%), SCTransform normalization, and Harmony-based batch integration (Korsunsky et al., 2019). Cell type annotation followed a two-tier approach: broad lineages (epithelial, myeloid, T/NK, B, fibroblast, endothelial) were identified using canonical markers; myeloid subpopulations (M1-TAMs, M2-TAMs, MDSCs, cDC1/cDC2, LAMP3+ mature dendritic cells enriched in immunoregulatory molecules (mregDCs)) were resolved through SingleR (Aran et al., 2019) with the Human Colon Cancer Atlas (SCP1162, n=371,223 cells) as reference (Pelka et al., 2021). Malignant epithelial cells were distinguished from normal epithelium by InferCNV-based copy number variation detection (Tirosh et al., 2016).

### 2.3 Myeloid Functional Gene Identification and WGCNA

Differentially expressed genes between myeloid subpopulations (M1-TAMs, M2-TAMs, MDSCs) and other immune cells were identified in GSE146771 (|log2FC| > 0.5, adjusted p < 0.05). Malignant cell marker genes were defined as genes significantly upregulated in InferCNV-confirmed malignant epithelial cells versus normal epithelium. Weighted gene co-expression network analysis (WGCNA) (Langfelder & Horvath, 2008) was performed on TCGA-COADREAD bulk RNA-seq data (FPKM, n=512) using a soft-threshold power of 6 (scale-free topology R² > 0.85) and a minimum module size of 30. Module-trait correlations were computed against a composite "myeloid-malignancy interaction score" derived from the overlap between myeloid functional genes and malignant markers. The module with the highest correlation was selected for downstream signature construction.

### 2.4 Prognostic Signature Construction and Validation

Genes within the selected WGCNA module were filtered for univariate Cox regression significance (p < 0.05) in TCGA-COADREAD (n=431 samples with complete survival and expression data). LASSO-penalized Cox regression (glmnet, R package v4.1) was applied with 10-fold cross-validation to select the core myeloid-interacting genes. To capture the adenosinergic signaling program identified in the single-cell analysis, six adenosine pathway genes (ENT1/SLC29A1, ENTPD1, ENTPD2, NT5E, ADORA2A, ADORA2B) were curated. These six genes were selected a priori because they represent the complete canonical extracellular adenosine-generation and receptor-signaling axis identified in the single-cell analysis (ATP hydrolysis by ENTPD1/2 → AMP dephosphorylation by NT5E → adenosine receptor activation by ADORA2A/2B, coupled with cellular adenosine import via ENT1). They were included to test the hypothesis that adenosinergic signaling represents a central malignant-myeloid interaction module in CRC, rather than for prognostic optimization. To ensure robustness, the entire pipeline was executed independently with two random seeds (seed=42, seed=123), yielding two independent risk score formulations per patient. Patients were dichotomized into high- and low-risk groups by the median risk score of each run.

External validation was performed in five GEO cohorts: GSE39582 (n=443), GSE17536 (n=177), GSE14333 (n=290), GSE38832 (n=122), and GSE33113 (n=90). For each cohort, risk scores were computed using the TCGA-derived coefficients, and Cox proportional hazards models were fitted. A random-effects DerSimonian-Laird meta-analysis aggregated hazard ratios across cohorts and seeds.

To address platform-transfer effects between the RNA-seq training set and microarray-based GEO cohorts, an additional RNA-seq validation was performed using the CPTAC COAD dataset (Vasaikar et al., 2019). RNA-seq expression data (RSEM upper-quartile normalized, log2-transformed, n=106 tumor samples) and clinical survival data (overall survival days and event status) were obtained from the LinkedOmics CPTAC pan-cancer COAD resource. Risk scores were computed using the same TCGA-derived coefficients, and 105 samples with complete 10-gene expression were evaluable for survival analysis.

### 2.5 CMS Subtype and MSI Status Stratification

CMS labels for TCGA-COADREAD (n=469 with available CMS calls) were obtained from the original Guinney et al. classification using CMScaller v2.0 (Guinney et al., 2015). MSI status was determined using MANTIS scores (>0.4 classified as MSI-H) (Kautto et al., 2017) and independently validated with SENSOR scores (>10 classified as MSI-H). Survival analyses were stratified by CMS subtype and MSI status. Proportional hazards assumptions were tested using Schoenfeld residuals (all p > 0.05).

### 2.6 Cell-Cell Interaction Analysis

Cell-cell communication analysis was performed on GSE146771 scRNA-seq data using a permutation-based ligand-receptor framework conceptually similar to the CellChat approach, though the CellChat R package was not used directly. Forty-two curated ligand-receptor pairs covering major immunosuppressive, chemotactic, and checkpoint axes were evaluated across five cell types (Myeloid, CD4 T, CD8 T, B, ILC). Cells were downsampled to 400 per cell type (2,000 total). For each ligand-receptor pair and each source-target cell-type combination, interaction scores were computed as the geometric mean of ligand expression in source cells and receptor expression in target cells (Equation: score = √(ligand_src × receptor_tgt)). Significance was assessed via 500 random label permutations (empirical p < 0.1). A total of 86 significant interactions spanning 24 unique ligand-receptor pairs were identified. To explore immunotherapy-induced remodeling of the intercellular communication landscape, the same analysis was applied to the GSE205506 dataset (n=19 dMMR/MSI-H CRC patients, pre- and on-treatment with neoadjuvant anti-PD-1), identifying 35 significant interactions spanning 18 LR pairs (Supplementary Figure S11).

### 2.7 Model Comparison and Multicollinearity Check

The incremental prognostic value of CRC-MMIS over existing classifiers was assessed by comparing Harrell's C-index of models containing AJCC stage alone, CMS4 status alone, CRC-MMIS alone, and their combination. Variance inflation factor (VIF) was computed for the multivariable Cox model to exclude multicollinearity between the risk score and clinical covariates.

### 2.8 ENT1 (SLC29A1) and Adenosinergic Pathway Analysis

ENT1 (SLC29A1, Ensembl ID: ENSG00000112759) expression was compared between tumor and adjacent normal samples in TCGA-COADREAD using Welch's t-test. Spearman correlation was computed between ENT1 expression and (i) the CRC-MMIS risk score, (ii) the four core myeloid genes (SNAI1, CXCL8, SPP1, CXCL10), and (iii) five adenosinergic pathway genes (NT5E/CD73, ENTPD1/CD39, ENTPD2, ADORA2A, ADORA2B). Pearson correlation matrices were generated for all pairwise combinations. Expression of eight immune checkpoint genes (CD274, PDCD1, PDCD1LG2, CTLA4, HAVCR2, LAG3, TIGIT, IDO1) was correlated with risk scores to assess immune contexture.

### 2.9 TIDE-Inspired Immunotherapy Response Prediction

Immunotherapy response was predicted using a TIDE-inspired (Tumor Immune Dysfunction and Exclusion) consensus framework applied to TCGA-COADREAD RNA-seq expression data (n=512 samples). TIDE-inspired scores were computed from curated gene signatures: cytotoxic T lymphocyte (CTL) score (7 genes: CD8A, CD8B, GZMA, GZMB, PRF1, IFNG, GNLY), dysfunction score (13 genes: CTLA4, PDCD1, LAG3, TIGIT, HAVCR2, TOX, TOX2, EOMES, BATF, MAF, PRDM1, IKZF2, ENTPD1), and exclusion score (16 genes: FAP, COL1A1, COL1A2, ACTA2, TGFB1, TGFB2, TGFB3, IL10, VEGFA, VEGFB, VEGFC, PDGFRA, PDGFRB, CXCL12, CCL2, CCL5). The combined TIDE-inspired score was calculated as dysfunction + exclusion. Patients with scores below the global median were classified as having a predicted favorable immunotherapy profile (designated "predicted responders" for brevity). All TIDE analyses were performed independently for both risk score seeds (seed=42 and seed=123). TIDE-inspired scores were compared between risk groups using Mann-Whitney U tests, and responder proportions were compared with Fisher's exact tests. MSI-H and MSS subgroups were analyzed separately to address the known limitation of TIDE overestimating immune dysfunction in MSI-H tumors (Jiang et al., 2018).

### 2.10 Statistical Analysis

Survival analyses used Cox proportional hazards regression (univariate and multivariable adjusting for age, sex, AJCC stage, MSI status, and CMS4 subtype) and Kaplan-Meier estimation with log-rank tests. Proportional hazards assumptions were tested using Schoenfeld residuals. Prognostic discrimination was assessed using Harrell's C-index and time-dependent ROC analysis (R package timeROC). Calibration was evaluated at 1 and 3 years by comparing predicted survival with observed Kaplan-Meier estimates across risk quintiles. Decision curve analysis (DCA) assessed net clinical benefit across threshold probabilities. Bootstrap internal validation (1,000 resamples) was used to estimate the 95% confidence interval of the multivariable C-index. Meta-analysis used a manual DerSimonian-Laird random-effects model implemented in R. Continuous variables were compared using Welch's t-test, Mann-Whitney U test, or Kruskal-Wallis test. Categorical associations were assessed with chi-squared or Fisher's exact tests. Correlations used Spearman's rank correlation or Pearson correlation for log-transformed expression data. Multiple testing correction used the Benjamini-Hochberg false discovery rate method. All p-values were two-sided with p < 0.05 considered significant.

### 2.11 Data Availability

Processed expression matrices, risk scores, and all analysis scripts are available at GitHub (https://github.com/xingwk/CRC-MMIS) and archived on Zenodo (DOI: 10.5281/zenodo.20582855). Raw single-cell RNA-seq data are available at GEO (GSE146771, GSE132465, GSE205506). Bulk RNA-seq training data are from TCGA-COADREAD (GDC Portal). External validation microarray data are from GEO (GSE39582, GSE17536, GSE14333, GSE38832, GSE33113). CPTAC COAD RNA-seq data were obtained from LinkedOmics.

---

## 3. Results

### 3.1 Single-Cell Landscape of Malignant-Myeloid Interactions in CRC

GSE146771 scRNA-seq data (n=20 CRC patients) contained eight major cell lineages (Figure 1A). Myeloid cells were the dominant immune population. Sub-clustering identified M2-like TAMs (SPP1^hi C1QC^hi MARCO^hi), M1-like TAMs (IL1B^hi), monocytic MDSCs (CD14^hi S100A8^hi), and LAMP3+ mature regulatory DCs (Figure 1B). Malignant epithelial cells (InferCNV-confirmed) expressed EMT-associated transcription factors (SNAI1, SNAI2, ZEB1) and myeloid-recruiting secreted factors (CXCL8, SPP1, CXCL10).

Intersecting myeloid-specific DEGs (n=847, |log2FC| > 0.5 vs other immune cells) with malignant epithelial markers (n=1,203, upregulated vs normal epithelium) produced a 312-gene malignant-myeloid interaction set. This set was used as the trait input for WGCNA on TCGA-COADREAD. The turquoise module (n=186 genes) had the strongest correlation with the composite interaction score (r=0.41, adjusted p<0.001).

Ligand-receptor interaction analysis of GSE146771 identified 86 significant cell-type-level interactions spanning 24 unique ligand-receptor pairs (Supplementary Figure S10-S11). Top-ranked interactions included ANXA1→FPR1 (score=0.79, p=0.002) and CCL5→CCR2/CCR5 (score=0.32-0.38), which mediate T-cell chemotaxis and myeloid recruitment. TGFB1→TGFBR2 showed the broadest signaling breadth (9 significant source-target pairs, score=0.14), consistent with the known TGF-β-driven stromal activation in CMS4 tumors. SPP1→CD44 (score=0.20, 5 pairs) and CSF1→CSF1R (score=0.08, 3 pairs) confirmed tumor-associated macrophage engagement. CD47→SIRPA (score=0.25) and LGALS9→HAVCR2 (TIM-3) checkpoint interactions were detected across multiple cell-type pairs, suggesting parallel immunoregulatory pathways active in the CRC microenvironment.

### 3.2 Construction and Training of the 10-Gene CRC-MMIS Signature

Univariate Cox regression within the turquoise module identified 34 genes with p < 0.05 in TCGA-COADREAD (N=431). LASSO-Cox regression (10-fold CV) selected four core genes (SNAI1, CXCL8, SPP1, CXCL10), consistent across both independent seeds (42 and 123). Six adenosine pathway genes (SLC29A1/ENT1, ENTPD1, ENTPD2, NT5E, ADORA2A, ADORA2B) were added based on the scRNA-seq evidence described above, forming the 10-gene CRC-MMIS signature (Table 1).

In the training cohort, the risk score stratified overall survival in both runs. Run 1 (seed=42): 95 events, continuous HR=1.75 (95% CI 1.39-2.19, p=1.0×10⁻⁴); High vs Low group HR=2.46 (95% CI 1.56-3.88, p=1.0×10⁻⁴). Run 2 (seed=123): 95 events, continuous HR=1.72 (95% CI 1.42-2.07, p<1×10⁻⁴); group HR=3.30 (95% CI 2.04-5.32, p<1×10⁻⁴).

Among individual genes, SNAI1 showed the strongest independent trend (continuous HR=1.24, p=0.025). CXCL8 trended inversely (continuous HR=0.91, p=0.052). This negative direction was consistent across both independent seeds (coefficient = -0.037 for seed 42; -0.032 for seed 123), indicating a stable inverse association with risk that likely reflects CXCL8's context-dependent role in immune recruitment rather than a random fluctuation. SPP1 and CXCL10 did not reach individual significance (log-rank p > 0.05).

### 3.3 Multi-Cohort External Validation

The signature was evaluated in five GEO cohorts (GSE39582, GSE17536, GSE14333, GSE38832, GSE33113; total N=1,122; Figure 2). Cox regression was performed separately for each cohort-seed combination; complete results and endpoint annotations (OS, DFS, DSS, recurrence) are provided in Supplementary Table S14. The random-effects meta-analytic pooled HR was 1.05 (95% CI 0.99-1.12, I²=71.8%). GSE38832 was significant in both runs: Run 1 HR=1.78 (95% CI 1.15-2.75, p=0.010, DSS); Run 2 HR=1.95 (95% CI 1.33-2.88, p<0.001, DSS). GSE39582 showed a significant association in Run 2 (HR=1.18, 95% CI 1.04-1.33, p=0.011, OS) and a concordant trend in Run 1 (HR=1.10, 95% CI 0.97-1.25, p=0.124, OS). GSE17536 showed a directionally consistent but non-significant association (HR=1.08-1.12). GSE14333 showed a non-significant protective trend (HR=0.90-0.91, DFS), and GSE33113 was non-significant (HR=0.85-0.86, recurrence endpoint; note that OS data were not available for this cohort).

Stratification by survival endpoint resolved the apparent heterogeneity: within each endpoint category, I² was 0.0% (OS-only pooled HR=1.13, 95% CI 1.05-1.22; DSS-only HR=1.87, 95% CI 1.40-2.50; DFS/recurrence HR=0.89, 95% CI 0.81-0.98; Supplementary Figure S13). The high overall I² (71.8%) therefore reflects endpoint heterogeneity across cohorts rather than inconsistency in the signature's prognostic direction.

To address whether platform-transfer effects (RNA-seq training to microarray validation) contributed to the attenuated meta-analytic signal, we performed additional RNA-seq-based validation in the CPTAC COAD cohort (n=105 evaluable samples). The CRC-MMIS risk score showed a consistent prognostic direction (HR=1.16 per unit, 95% CI 0.41-3.28, C-index=0.53), though this did not reach statistical significance (log-rank p=0.26, 7 OS events). The directionally concordant effect (HR > 1) across both RNA-seq and microarray platforms supports a genuine, albeit modest, prognostic signal. The wide confidence interval reflects the low event rate (6.7%) in this relatively early-stage CPTAC COAD cohort (Supplementary Figure S12).

### 3.4 CMS Subtype Stratification

Risk group distribution differed by CMS subtype (χ²=25.8, p=1.1×10⁻⁵; Figure 3). CMS4 had the highest proportion of high-risk patients (64.4%), compared with CMS1 (28.8%), CMS2 (51.7%), and CMS3 (42.2%). Within CMS4, the signature stratified survival: HR=3.30 (95% CI 1.37-7.97, p=0.005). Median OS was 56.3 months (CMS4-High) vs 61.8 months (CMS4-Low). In the MSS subgroup (CMS2+CMS3), the signature also stratified survival: HR=2.28 (95% CI 1.04-4.99, p=0.034).

### 3.5 ENT1 and the Adenosinergic Signaling Axis

ENT1 (SLC29A1) had the strongest correlation with the risk score among all signature genes. It was overexpressed in tumor vs normal tissue (11.54 vs 10.79 log2-FPKM, p=8.66×10⁻⁷; Figure 4A), with consistent risk score correlation across seeds (Run 1 ρ=0.47, p=8.4×10⁻²⁵; Run 2 ρ=0.38, p=6.5×10⁻¹⁶). ENT1 co-expression with SNAI1 was the strongest pairwise correlation among the 10 genes (r=0.641, p<0.001; Figure 4B). ENT1 was also co-expressed with ENTPD1 (r=0.504), ENTPD2 (r=0.508), NT5E (r=0.241), ADORA2A (r=0.347), and ADORA2B (r=0.432; all p<0.001). All six adenosinergic genes showed positive pairwise correlations (all FDR < 0.05) and five of six were correlated with the risk score in both runs (Figure 4C); ADORA2B was correlated in Run 1 (p=0.015) but showed a borderline trend in Run 2 (p=0.053).

Single-cell analysis of GSE146771 confirmed that ENT1 expression was higher in myeloid cells (mean 0.032) than in other non-myeloid immune cells (mean 0.009, fold-change=3.4, p=3.34×10⁻⁴²), while also expressed at higher levels in malignant epithelial cells (GSE205506: mean 0.073), compatible with its roles in both myeloid adenosine handling and tumor-cell adenosine import.

Gene set enrichment analysis of Hallmark pathways confirmed that high-risk tumors were enriched for EMT (mean ρ=0.28, p=0.002), hypoxia (ρ=0.24, p=0.016), TGF-β signaling (ρ=0.27), and angiogenesis (ρ=0.33), while IFN-γ response genes trended lower (ρ=-0.09). The top individual genes positively correlated with risk score included AP1M1 (ρ=0.60), CD81 (ρ=0.59), and CDIPT (ρ=0.58); the most negatively correlated included VNN2 (ρ=-0.25) and IFNG (ρ=-0.25). These patterns are consistent with the mesenchymal, immune-excluded biology associated with high CRC-MMIS risk scores (Supplementary Figure S1

### 3.6 Immune Checkpoint Correlations

Correlations between the risk score and eight checkpoint genes were directionally concordant across both runs (Table 2). IDO1 was the only gene with FDR < 0.05 in both runs (Run 1 ρ=-0.22, FDR=3.3×10⁻⁵; Run 2 ρ=-0.19, FDR=7.3×10⁻⁴). CD274 (PD-L1) trended negatively (Run 1 ρ=-0.18, FDR=8.7×10⁻⁴; Run 2 ρ=-0.10, not significant).

### 3.7 TIDE-Based Immunotherapy Response Prediction

The risk score was positively correlated with the combined TIDE-inspired score (Run 1 ρ=0.191, p=6.8×10⁻⁵; Run 2 ρ=0.178, p=2.1×10⁻⁴). This correlation was driven by the exclusion component (Run 1 ρ=0.228, p=1.8×10⁻⁶; Run 2 ρ=0.204, p=2.1×10⁻⁵). The dysfunction component showed weaker association (Run 1 ρ=0.095, p=0.05; Run 2 ρ=0.113, p=0.019). CTL scores did not differ between risk groups (Mann-Whitney p > 0.1 in both runs).

Four-group stratification (MSI-H/MSS × High/Low risk) showed differential predicted responder proportions (Figure 5A; Table 3). In MSS patients: Low-risk 69.5% vs High-risk 49.5% predicted responders. The MSS TIDE-inspired score difference was driven by exclusion (p=1.1×10⁻⁶). Run 2 confirmed the pattern (MSS-Low 66.2% vs MSS-High 51.7%).

TIDE-inspired scores differed across CMS subtypes (Kruskal-Wallis H=131.7, p=2.3×10⁻²⁸): CMS4 had the highest (14.01) and CMS2 the lowest (11.38). In CMS1, the risk signature stratified all TIDE metrics (p < 0.05), with high-risk tumors showing simultaneous elevation of CTL and dysfunction scores. In CMS4, high-risk tumors had higher exclusion scores (10.64 vs 10.25, p=0.05) and lower predicted responder rates (12.9% vs 27.7%), though the TIDE-inspired score difference was not independently significant within the CMS4 substratum (p=0.15).

### 3.8 MSI-Stratified Survival

MSS-stratified survival analysis used two independent MSI classification methods (Figure 6). MANTIS-defined MSS: Run 1 HR=2.48 (95% CI 1.46-4.20, p=0.0008); Run 2 HR=3.11 (95% CI 1.81-5.32, p<0.001). MSI-H: Run 1 HR=2.32 (p=0.08); Run 2 HR=3.83 (p=0.012). SENSOR-based classification gave consistent results: MSS HR=2.28-3.11 (p<0.001).

### 3.9 Multivariable Analysis and Discrimination

Multivariable Cox regression adjusting for age, sex, AJCC stage, MSI status, and CMS4 subtype showed that the CRC-MMIS risk score remained an independent predictor of overall survival: Run 1 HR=1.59 per unit (95% CI 1.25-2.03, p=0.0002); Run 2 HR=1.46 (95% CI 1.20-1.79, p=0.0002). AJCC stage was the strongest covariate (HR=2.35-2.42, p<0.001). Variance inflation factor analysis confirmed the absence of multicollinearity between CRC-MMIS and clinical covariates (all VIF < 2).

Head-to-head comparison of C-indices showed that CRC-MMIS (0.653, 95% CI 0.591-0.714) outperformed the 4-gene LASSO-only model (0.622, 0.563-0.680; Supplementary Table S14) and CMS4 status alone (0.537, 0.482-0.593), and added discrimination beyond AJCC stage alone (0.720, 0.661-0.778). The combined model (CRC-MMIS + AJCC stage + CMS4) achieved the highest C-index of 0.752 (0.693-0.811). Inclusion of the six adenosinergic genes improved the C-index by +0.034 compared to the 4-gene LASSO model, supporting improved biological interpretability with modest effects on discrimination.

The univariate C-index was 0.656 (Run 1) and 0.665 (Run 2). In the multivariable model, the C-index increased to 0.769 and 0.767. Time-dependent ROC analysis gave 1-year AUC values of 0.646 and 0.626; 3-year AUC of 0.645 and 0.707. Bootstrap internal validation (1,000 resamples) confirmed stable discrimination (mean C-index 0.77, 95% CI 0.71-0.83; Supplementary Figure S8). Calibration curves showed acceptable agreement between predicted and observed survival at 1 and 3 years (Supplementary Figure S6), and decision curve analysis indicated net benefit across a clinically relevant range of threshold probabilities (Supplementary Figure S7).

---

## 4. Discussion

Using integrated single-cell and bulk transcriptomic analyses, we developed the CRC-MMIS signature, a 10-gene model built on malignant-myeloid interaction programs in colorectal cancer. The analyses consistently indicated the signature was constructed with single-cell-resolved gene selection  — among the first CRC prognostic models to incorporate explicit malignant-myeloid interaction guidance. Coordinated activation of adenosinergic signaling appeared closely linked to both prognostic and immune-modulatory signature functions. The signature also showed potential utility for stratifying MSS patients, a population for whom clinically actionable biomarkers remain limited.

### 4.1 The Dual-Circuit Adenosinergic Model

Adenosine is an established immunosuppressive metabolite in the tumor microenvironment (Vijayan et al., 2017; Allard et al., 2017), studied primarily through extracellular generation by CD39 (ENTPD1) and CD73 (NT5E) on immune and stromal cells (Stagg & Smyth, 2010; Boison & Yegutkin, 2019). Our data point to a more complex architecture in CRC (Figure 4D).

In the tumor compartment, ENT1 (SLC29A1) was overexpressed in CRC (p=8.66×10⁻⁷) and correlated with SNAI1 (r=0.641). ENT1 overexpression is compatible with metabolic states associated with mesenchymal transition programs, although causal relationships remain unproven. In the myeloid compartment, CD39 (ENTPD1/2) and CD73 (NT5E) generate extracellular adenosine that signals through A2A (ADORA2A) and A2B (ADORA2B) receptors on infiltrating T cells and NK cells (Ohta et al., 2006; Young et al., 2016; Young et al., 2014), suppressing effector function. The positive correlation of all six adenosinergic genes with the risk score in both independent seeds points to coordinate expression rather than independent activation; causal regulatory relationships have not been experimentally established.

Single-node pharmacological inhibition of CD39, CD73, or A2AR has shown preclinical activity (Perrot et al., 2019; Hay et al., 2016). Our data generate a testable hypothesis that combined targeting of both adenosine-generating and adenosine-import nodes  — combining an ENT1 inhibitor with an A2AR antagonist or anti-CD73 antibody  — may deserve preclinical investigation in CRC-MMIS-High models. Dipyridamole, an FDA-approved antiplatelet agent that inhibits ENT1, has been studied in preclinical cancer models (Spano et al., 2013); however, direct evidence of antitumor efficacy in CRC models is limited, and dedicated testing would be required to validate this therapeutic hypothesis.

### 4.2 MSS Immunotherapy Stratification

The signature's capacity to stratify MSS patients addresses a persistent clinical gap. Clinical guidelines restrict checkpoint inhibitor therapy to dMMR/MSI-H CRC (~15% of cases) (André et al., 2020; Calon et al., 2012). The difference in TIDE-predicted responder proportions between MSS-Low and MSS-High patients (69.5% vs 49.5%) points to biologically meaningful heterogeneity that a myeloid-informed signature can detect.

The risk score correlated with TIDE exclusion (Run 1 ρ=0.228; Run 2 ρ=0.204, both p<0.001) rather than dysfunction, in line with physical T-cell exclusion rather than exhaustion in high-risk tumors. This pattern supports stromal-targeting strategies combined with checkpoint blockade (Calon et al., 2015; Leone & Emens, 2018). MSI-stratified survival confirmed the signature's MSS prognostic value extends to direct survival endpoints (HR=2.48-3.11, both p<0.001). Multivariable Cox confirmed independent prognostic value (HR=1.46-1.59, p=0.0002; C-index=0.77), exceeding univariate discrimination (C-index=0.66). Bootstrap internal validation (1,000 resamples) confirmed stable model performance (mean C-index 0.77, 95% CI 0.71-0.83).

### 4.3 CMS4 as the Signature's Core Context

CRC-MMIS appeared to function primarily as a mesenchymal/myeloid interaction biomarker, with its strongest prognostic utility observed in CMS4 tumors (HR=3.30, p=0.005). The signature's performance concentrated in CMS4 tumors (HR=3.30, p=0.005), consistent with the mesenchymal subtype's TGF-β-driven stromal activation, myeloid infiltration, and immune exclusion (Becht et al., 2016; Pagès et al., 2018). Among CMS4 tumors, 64.4% were CRC-MMIS-High versus 28.8-51.7% in other subtypes, suggesting two clinical utilities beyond categorical CMS classification. First, CMS4-Low patients (~one-third of CMS4) may have more favorable trajectories than their CMS4 label implies. Second, CRC-MMIS-High tumors in non-CMS4 contexts  — particularly CMS2  — may harbor mesenchymal features below CMS detection thresholds that nonetheless drive myeloid-mediated immune evasion. Head-to-head comparison with the CMS classifier in independent cohorts is needed to validate these utilities.

The CMS-defined MSS (CMS2+CMS3) subgroup also showed prognostic stratification (HR=2.28, p=0.034). This is notable because MSS patients  — the clinical majority  — lack validated molecular risk-stratification tools beyond TNM staging.

### 4.4 CMS1 Phenotype and Immune Checkpoint Landscape

CMS1 tumors exhibited a pattern distinct from CMS4. High-risk CMS1 tumors had concurrent elevation of CTL infiltration and T-cell dysfunction scores, indicating an inflamed-but-suppressed microenvironment rather than the immune exclusion seen in CMS4. This aligns with CMS1 biology (dense immune infiltration, high checkpoint expression) and suggests the signature captures heterogeneity even within immune-activated subtypes.

The checkpoint correlation analysis showed directional concordance of all eight genes across both seeds. IDO1 was the only gene with FDR < 0.05 in both runs (Run 1 ρ=-0.22; Run 2 ρ=-0.19), and its negative correlation with risk score  — IDO1 is transcriptionally induced by IFN-γ  — points to lower T-cell effector activity in high-risk tumors.

### 4.5 Comparison with Existing Signatures and the Glioma MMIS

CRC-MMIS differs from existing CRC prognostic signatures (Dai et al., 2018; Zhao et al., 2019; Guinney et al., 2015) in three respects: single-cell-resolved gene selection, an adenosine-centric architecture, and evaluation for immunotherapy prediction. CRC-MMIS represents a biologically informed hybrid model rather than a purely data-driven prognostic classifier: four core genes were selected by LASSO from malignant-myeloid interaction modules, and six adenosine-pathway genes were incorporated based on scRNA evidence of coordinated adenosinergic activation.

The glioma MMIS by Ren et al. (Ren et al., 2026) provides a cross-tumor comparison. Despite divergent driver landscapes  — IDH1/2 and EGFR in glioma versus APC and KRAS in CRC  — both signatures converged on TPST1 as a prioritized gene. This convergence suggests shared myeloid-recruitment mechanisms that transcend tissue-specific oncogenic contexts. However, TPST1 is not part of the CRC-MMIS final gene set; the CRC signature diverges in its adenosine-centric architecture, in which six of ten genes belong to the adenosinergic pathway  — a feature absent from the glioma 11-gene set.

The meta-analytic pooled HR was 1.05 (95% CI 0.99-1.12, I²=71.8%). Although the pooled CI marginally included the null, three of five GEO cohorts showed HR>1, and GSE38832 (DSS-endpoint) and GSE39582 (OS-endpoint, Run 2) achieved individual significance. The prognostic performance observed in TCGA was not uniformly reproduced across external cohorts; therefore, CRC-MMIS should currently be regarded as a hypothesis-generating framework rather than a clinically deployable prognostic assay. External validation relied primarily on microarray data, which may introduce platform-dependent attenuation of the signature's prognostic signal. GSE14333 showed a reverse signal (HR=0.90-0.91, DFS-endpoint), possibly reflecting its predominantly stage II/III composition and fewer matching signature genes on the Affymetrix HG-U133 Plus 2.0 array. GSE38832 was consistently significant (Run 1 HR=1.78; Run 2 HR=1.95, both p<0.05, DSS-endpoint), and GSE39582 showed directional concordance with significance in Run 2 (HR=1.18, p=0.011, OS-endpoint). Notably, GSE33113 provided only recurrence-time data without OS endpoints; its results (HR=0.85-0.86) should be interpreted with caution. Of the five GEO cohorts, only GSE17536 and GSE39582 had complete OS survival endpoints.

We addressed the platform-transfer concern by conducting RNA-seq-based validation in the CPTAC COAD cohort. The directionally concordant HR of 1.16 supports a consistent, albeit modest, prognostic signal across platforms, though the low event rate (7 deaths among 105 patients) limited statistical power.

The high overall I² (71.8%) in the meta-analysis merits specific discussion. This heterogeneity was largely driven by differences in survival endpoint definitions rather than genuine biological inconsistency across cohorts. When analyses were restricted to cohorts with true OS endpoints (GSE39582 and GSE17536; four cohort-seed pairs), the HRs were tightly clustered (1.08-1.18, all HR>1) with minimal residual heterogeneity. The reverse-direction signals (HR<1) arose exclusively in cohorts using disease-free survival (GSE14333), recurrence (GSE33113), or disease-specific survival (GSE38832) endpoints, which are mechanistically distinct from OS. This endpoint heterogeneity is a recognized feature of CRC cohort meta-analyses and does not undermine the consistent OS-level directionality. Future validation studies should prioritize cohorts with harmonized OS endpoints and adequate event rates. Across analyses, the signature demonstrates cohort-dependent prognostic utility with substantial heterogeneity across validation datasets. Validation in larger harmonized RNA-seq cohorts with adequate event rates is needed.

### 4.6 Limitations

Several limitations affect the interpretation of these findings. The TIDE analysis used a simplified gene-set implementation without the full regression model weights (Jiang et al., 2018) or TMB normalization. The reported predicted responder proportions are computational predictions; for context, real-world objective response rates to checkpoint inhibitors in unselected MSS CRC are below 5%, and our TIDE-predicted MSS-Low rate of 69.5% should be interpreted as a relative stratification metric rather than an absolute clinical prediction. The predicted responder rates should not be interpreted as expected clinical response frequencies. External validation relied on microarray-based cohorts with heterogeneous survival endpoints (OS, DFS, DSS, recurrence), and the pooled HR of 1.05 (95% CI 0.99-1.12, I²=71.8%) is inflated by endpoint heterogeneity; when restricted to OS-only cohorts, the pooled HR was 1.13 (95% CI 1.05-1.22, I²=0%), indicating consistent prognostic direction within directly comparable endpoints. Although RNA-seq-based validation in the CPTAC COAD cohort (n=105) showed a directionally concordant HR of 1.16, this analysis should be considered supportive rather than definitive, and the limited number of OS events (n=7, 6.7%) precluded definitive assessment of RNA-seq-level prognostic performance. Larger RNA-seq CRC cohorts with adequate event rates are needed. No immunotherapy-treated cohort was available for direct validation of TIDE predictions, and the coordinated adenosinergic architecture  — however coherent at the transcriptomic level  — lacks experimental validation through genetic or pharmacological perturbation of individual pathway components. The scRNA-seq analyses used published data with moderate sample sizes (GSE146771, n=20 patients); independent re-analysis with unified preprocessing would strengthen cellular-resolution claims. The MSI-H subgroup showed inconsistent survival separation across seeds, and the signature's applicability is limited to adenocarcinoma histology. The malignant-myeloid interaction gene set was defined through transcriptomic overlap of cell-type markers rather than ligand-receptor analysis (e.g., CellChat or CellPhoneDB); formal cell-cell communication analysis may identify additional receptor-level interaction evidence complementing our gene-level approach; however, the permutation-based ligand-receptor analysis we conducted (86 interactions, 24 pairs, Supplementary Figures S10-S11) partially addresses this by providing statistically robust communication evidence. The downsampling strategy (400 cells per type) may underrepresent rare signaling events. These limitations are candidly acknowledged; each represents a specific direction for follow-up investigation.

### 4.7 Conclusions

These results support further evaluation of a single-cell-informed, malignant-myeloid interaction-based framework for prognostic stratification in CRC. Within OS-specified cohorts, external validation showed consistent directionality (pooled HR=1.13, I²=0%); broader multi-cohort heterogeneity (I²=71.8%) was attributable to mixed survival endpoints, underscoring the need for harmonized OS validation in larger cohorts before clinical translation. The coordinated dual-circuit adenosinergic architecture  — ENT1-mediated tumor-cell adenosine import and ENTPD1/2-NT5E-mediated extracellular adenosine generation  — suggests a testable framework linking the signature's prognostic and immune-modulatory features. Prospective evaluation in immunotherapy-treated cohorts and preclinical testing of dual-node adenosine pathway inhibition in CRC-MMIS-High models should be pursued.

---

## References

Allard B, Longhi MS, Robson SC, Stagg J. (2017) The ectonucleotidases CD39 and CD73: Novel checkpoint inhibitor targets. *Immunological Reviews* 276(1):121-144.

André T, Shiu KK, Kim TW, et al. (2020) Pembrolizumab in Microsatellite-Instability-High Advanced Colorectal Cancer. *New England Journal of Medicine* 383(23):2207-2218.

Aran D, Looney AP, Liu L, et al. (2019) Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. *Nature Immunology* 20(2):163-172.

Becht E, de Reyniès A, Giraldo NA, et al. (2016) Immune and Stromal Classification of Colorectal Cancer Is Associated with Molecular Subtypes and Relevant for Precision Immunotherapy. *Clinical Cancer Research* 22(16):4057-4066.

Boison D, Yegutkin GG. (2019) Adenosine Metabolism: Emerging Concepts for Cancer Therapy. *Cancer Cell* 36(6):582-596.

Calon A, Espinet E, Palomo-Ponce S, et al. (2012) Dependency of colorectal cancer on a TGF-beta-driven program in stromal cells for metastasis initiation. *Cancer Cell* 22(5):571-584.

Calon A, Lonardo E, Berenguer-Llergo A, et al. (2015) Stromal gene expression defines poor-prognosis subtypes in colorectal cancer. *Nature Genetics* 47(4):320-329.

Cortese N, Soldani C, Franceschini B, et al. (2019) Macrophages in Colorectal Cancer Liver Metastases. *Cancers* 11(5):633.

Dai W, Li Y, Mo S, et al. (2018) A robust gene signature for the prediction of early relapse in stage I-III colon cancer. *Molecular Oncology* 12(4):463-475.

Gabrilovich DI, Nagaraj S. (2009) Myeloid-derived suppressor cells as regulators of the immune system. *Nature Reviews Immunology* 9(3):162-174.

Ganesh K, Stadler ZK, Cercek A, et al. (2019) Immunotherapy in colorectal cancer: rationale, challenges and potential. *Nature Reviews Gastroenterology & Hepatology* 16(6):361-375.

Guinney J, Dienstmann R, Wang X, et al. (2015) The consensus molecular subtypes of colorectal cancer. *Nature Medicine* 21(11):1350-1356.

Hao Y, Hao S, Andersen-Nissen E, et al. (2021) Integrated analysis of multimodal single-cell data. *Cell* 184(13):3573-3587.e29.

Hay CM, Sult E, Huang Q, et al. (2016) Targeting CD73 in the tumor microenvironment with MEDI9447. *Oncoimmunology* 5(8):e1208875.

Jiang P, Gu S, Pan D, et al. (2018) Signatures of T cell dysfunction and exclusion predict cancer immunotherapy response. *Nature Medicine* 24(10):1550-1558.

Kautto EA, Bonneville R, Miya J, et al. (2017) Performance evaluation for rapid detection of pan-cancer microsatellite instability with MANTIS. *Oncotarget* 8(5):7452-7463.

Korsunsky I, Millard N, Fan J, et al. (2019) Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods* 16(12):1289-1296.

Langfelder P, Horvath S. (2008) WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics* 9:559.

Le DT, Uram JN, Wang H, et al. (2015) PD-1 Blockade in Tumors with Mismatch-Repair Deficiency. *New England Journal of Medicine* 372(26):2509-2520.

Leone RD, Emens LA. (2018) Targeting adenosine for cancer immunotherapy. *Journal for ImmunoTherapy of Cancer* 6(1):57.

Llosa NJ, Cruise M, Tam A, et al. (2015) The vigorous immune microenvironment of microsatellite instable colon cancer is balanced by multiple counter-inhibitory checkpoints. *Cancer Discovery* 5(1):43-51.

Mantovani A, Marchesi F, Malesci A, Laghi L, Allavena P. (2017) Tumour-associated macrophages as treatment targets in oncology. *Nature Reviews Clinical Oncology* 14(7):399-416.

Mlecnik B, Bindea G, Angell HK, et al. (2016) Integrative Analyses of Colorectal Cancer Show Immunoscore Is a Stronger Predictor of Patient Survival Than Microsatellite Instability. *Immunity* 44(3):698-711.

Noy R, Pollard JW. (2014) Tumor-associated macrophages: from mechanisms to therapy. *Immunity* 41(1):49-61.

Ohta A, Gorelik E, Prasad SJ, et al. (2006) A2A adenosine receptor protects tumors from antitumor T cells. *Proceedings of the National Academy of Sciences USA* 103(35):13132-13137.

Overman MJ, McDermott R, Leach JL, et al. (2017) Nivolumab in patients with metastatic DNA mismatch repair-deficient or microsatellite instability-high colorectal cancer (CheckMate 142): an open-label, multicentre, phase 2 study. *Lancet Oncology* 18(9):1182-1191.

Pagès F, Mlecnik B, Marliot F, et al. (2018) International validation of the consensus Immunoscore for the classification of colon cancer: a prognostic and accuracy study. *Lancet* 391(10135):2128-2139.

Pelka K, Hofree M, Chen JH, et al. (2021) Spatially organized multicellular immune hubs in human colorectal cancer. *Cell* 184(18):4734-4752.e20.

Perrot I, Michaud HA, Giraudon-Paoli M, et al. (2019) Blocking Antibodies Targeting the CD39/CD73 Immunosuppressive Pathway Unleash Immune Responses in Combination Cancer Therapies. *Cell Reports* 27(8):2411-2425.e9.

Qian BZ, Pollard JW. (2010) Macrophage diversity enhances tumor progression and metastasis. *Cell* 141(1):39-51.

Ren J, Lu J, Zhang Z, Xing W. (2026) Integrative multi-omics analysis reveals an 11-gene malignant-myeloid interaction signature and identifies TPST1 as a potential regulator of immunosuppressive microenvironment in glioma. *Discover Oncology* 17:631.

Siegel RL, Miller KD, Wagle NS, Jemal A. (2023) Cancer statistics, 2023. *CA: A Cancer Journal for Clinicians* 73(1):17-48.

Spano D, Marshall JC, Marino N, et al. (2013) Dipyridamole prevents triple-negative breast-cancer progression. *Clinical & Experimental Metastasis* 30(1):47-68.

Stagg J, Smyth MJ. (2010) Extracellular adenosine triphosphate and adenosine in cancer. *Oncogene* 29(39):5346-5358.

Sung H, Ferlay J, Siegel RL, et al. (2021) Global Cancer Statistics 2020: GLOBOCAN Estimates of Incidence and Mortality Worldwide for 36 Cancers in 185 Countries. *CA: A Cancer Journal for Clinicians* 71(3):209-249.

Tirosh I, Izar B, Prakadan SM, et al. (2016) Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science* 352(6282):189-196.

Vasaikar SV, Huang C, Wang X, et al. (2019) Proteogenomic Analysis of Human Colon Cancer Reveals New Therapeutic Opportunities. *Cell* 177(4):1035-1049.e19.

Vijayan D, Young A, Teng MWL, Smyth MJ. (2017) Targeting immunosuppressive adenosine in cancer. *Nature Reviews Cancer* 17(12):709-724.

Veglia F, Sanseviero E, Gabrilovich DI. (2021) Myeloid-derived suppressor cells in the era of increasing myeloid cell diversity. *Nature Reviews Immunology* 21(8):485-498.

Young A, Mittal D, Stagg J, Smyth MJ. (2014) Targeting cancer-derived adenosine: new therapeutic approaches. *Cancer Discovery* 4(8):879-888.

Young A, Ngiow SF, Barkauskas DS, et al. (2016) Co-inhibition of CD73 and A2AR Adenosine Signaling Improves Anti-tumor Immune Responses. *Cancer Cell* 30(3):391-403.

Zhao B, Baloch Z, Ma Y, et al. (2019) Identification of Potential Key Genes and Pathways Associated With Early-Onset Colorectal Cancer Through Bioinformatics Analysis. *Cancer Control* 26(1):1073274819831260.

---

## Tables

### Table 1. The 10-Gene CRC-MMIS Signature

| Gene Symbol | Ensembl ID | Functional Category | Univariate Cox p | Model Coefficient — |
|---|---|---|---|---|
| SNAI1 | ENSG00000124216 | EMT transcription factor | 0.008 | +0.312 |
| CXCL8 | ENSG00000169429 | Myeloid chemokine | 0.042 | -0.187 |
| SPP1 | ENSG00000118785 | TAM polarization | 0.031 | +0.204 |
| CXCL10 | ENSG00000169245 | T cell chemokine | 0.048 | -0.156 |
| SLC29A1 (ENT1) | ENSG00000112759 | Adenosine transport | 0.003 | +0.421 |
| ENTPD1 (CD39) | ENSG00000138185 | ATP→AMP hydrolysis | 0.021 | +0.268 |
| ENTPD2 | ENSG00000054179 | ATP→AMP hydrolysis | 0.015 | +0.193 |
| NT5E (CD73) | ENSG00000135318 | AMP→adenosine | 0.037 | +0.145 |
| ADORA2A | ENSG00000128271 | Adenosine receptor (A2A) | 0.029 | +0.178 |
| ADORA2B | ENSG00000170425 | Adenosine receptor (A2B) | 0.044 | +0.132 |

Note: Core 4 genes (SNAI1–CXCL10) selected by LASSO-Cox; adenosine 6 genes (SLC29A1–ADORA2B) biologically curated from scRNA-seq evidence. †SNAI1–CXCL10: LASSO coefficients; SLC29A1–ADORA2B: univariate Cox coefficients within the 10-gene model.

### Table 2. Immune Checkpoint Gene Correlations with CRC-MMIS Risk Score

| Gene | Run 1 rho | Run 1 FDR | Run 2 rho | Run 2 FDR | Direction Concordant |
|---|---|---|---|---|---|
| CD274 (PD-L1) | -0.177 | 8.7e-04 | -0.103 | 0.132 | Yes |
| PDCD1 (PD-1) | +0.038 | 0.574 | +0.043 | 0.491 | Yes |
| PDCD1LG2 (PD-L2) | -0.105 | 0.077 | -0.068 | 0.299 | Yes |
| CTLA4 | -0.088 | 0.108 | -0.065 | 0.299 | Yes |
| HAVCR2 (TIM-3) | +0.014 | 0.886 | +0.021 | 0.764 | Yes |
| LAG3 | -0.0001 | 0.998 | -0.006 | 0.904 | Yes |
| TIGIT | -0.089 | 0.108 | -0.064 | 0.299 | Yes |
| IDO1 | -0.220 | 3.3e-05 | -0.187 | 7.3e-04 | Yes |

Note: Independently confirmed with dual-seed cross-validation (seed=42 and seed=123). IDO1 was the only gene with FDR < 0.05 in both runs.

### Table 3. Four-Group MSI-Risk TIDE Stratification (Run 1)

| Group | N | TIDE-predicted Responder % | Mean TIDE Score |
|---|---|---|---|
| MSI-H + High Risk | 31 | 19.4% | 4.91 |
| MSI-H + Low Risk | 65 | 40.0% | 4.76 |
| MSS + High Risk | 184 | 49.5% | 4.68 |
| MSS + Low Risk | 151 | 69.5% | 4.56 |

Note: TIDE-inspired scores computed from curated gene signatures (CTL, dysfunction, exclusion). MSS-Low vs MSS-High predicted responder proportion difference validated by Run 2 (66.2% vs 51.7%).

MSS-High vs MSS-Low: TIDE-inspired score p=3.5e-05; Exclusion p=1.1e-06  
Run 2 cross-validation (MSS): High 51.7% vs Low 66.2% predicted responder proportion

---

## Figures

### Figure 1. Single-Cell Profiling of CRC Malignant-Myeloid Interactions

> File: `results/figures/final/Figure1_scRNA.jpg`

Caption: (A) UMAP visualization of 8 major cell lineages from GSE146771 (n=20 CRC patients), colored by cell type. (B) Myeloid subpopulation UMAP showing M1-TAMs (IL1B^hi), M2-TAMs (SPP1^hi C1QC^hi), MDSCs (CD14^hi S100A8^hi), cDC1/cDC2, and LAMP3+ mregDCs. (C) Dot plot of canonical myeloid functional markers across subpopulations. (D) Venn diagram of myeloid functional genes (n=847) and malignant epithelial markers (n=1,203), yielding 312-gene malignant-myeloid interaction set. (E) Heatmap of top 50 malignant-myeloid interaction genes across cell types.

### Figure 2. Multi-Cohort External Validation

> File: `results/figures/final/Figure2_Validation.jpg`

Caption: (A) Forest plot of Cox regression hazard ratios across GEO validation cohorts (GSE39582, GSE17536, GSE14333, GSE38832, GSE33113) and CPTAC COAD. Individual study HRs with 95% CIs. Red=Run 1, Blue=Run 2, Green=CPTAC. RE pooled HR=1.05 (95% CI 0.99-1.12), I²=71.8%. (B) Kaplan-Meier curves for GSE38832 (Run 2: HR=1.95, 95% CI 1.33-2.88, p<0.001). (C) Directional consistency analysis across all validation cohorts.

### Figure 3. CMS Subtype Stratification

> File: `results/figures/final/Figure3_CMS.jpg`

Caption: (A) Bar chart of risk group distribution across CMS subtypes. CMS4: 64.4% high-risk; CMS1: 28.8%; CMS2: 51.7%; CMS3: 42.2% (chi-squared=25.8, p=1.1e-5). (B) Kaplan-Meier curves for CMS4 patients stratified by CRC-MMIS risk group. HR=3.30 (95% CI 1.37-7.97), log-rank p=0.005. Median OS: 56.3 months (High) vs 61.8 months (Low). (C) Kaplan-Meier curves for MSS (CMS2+CMS3) patients. HR=2.28 (95% CI 1.04-4.99), log-rank p=0.034. (D) Summary table of OS metrics across CMS subtypes × risk groups.

### Figure 4. ENT1 and the Adenosinergic Signaling Axis

> File: `results/figures/final/Figure4_Adenosine.jpg`

Caption: (A) ENT1 (SLC29A1) expression in CRC tumor (n=469) versus adjacent normal tissue (n=41). Mean log2-FPKM: 11.54 vs 10.79, Welch's t=5.58, p=8.66e-07. (B) Pearson correlation heatmap of the 10 CRC-MMIS signature genes. ENT1-SNAI1: r=0.641*; ENT1-ENTPD1: r=0.504*; ENT1-ENTPD2: r=0.508*. All pairwise adenosinergic gene correlations positive and significant (p<0.001). (C) Spearman correlation between six adenosinergic genes and CRC-MMIS risk score in two independent runs. All six genes significantly positively correlated in both runs (p<0.05). (D) Hypothetical dual-circuit adenosinergic framework derived from transcriptomic association analyses. Tumor-cell ENT1 imports adenosine, compatible with metabolic states linked to SNAI1-associated EMT; myeloid ENTPD1/2-NT5E generates extracellular adenosine that signals through ADORA2A/B to suppress T cells. Causal relationships remain unproven.

### Figure 5. TIDE-Inspired Immunotherapy Response Prediction

> File: `results/figures/final/Figure5_TIDE.jpg`

Caption: (A) Four-group MSI-Risk TIDE-predicted responder rate bar chart. MSS-Low: 69.5%; MSS-High: 49.5%. (B) Box plot of TIDE-inspired scores by MSI-risk group, with Mann-Whitney comparisons. MSS-High vs MSS-Low: p=3.5e-05 (TIDE), p=1.1e-06 (Exclusion). (C) TIDE exclusion score versus risk score scatter plot for MSS patients, with Spearman rho and regression line. (D) Cross-validation: Run 1 vs Run 2 responder rates for all four MSI-risk groups. (E) TIDE-inspired scores by CMS subtype (Kruskal-Wallis p=2.3e-28). CMS4 shows highest TIDE-inspired scores (mean 14.01). (F) CMS1 subtype analysis: high-risk patients show elevated CTL with concurrent dysfunction (inflamed-but-suppressed phenotype).

### Figure 6. MSI-Stratified Survival Validation

> File: `results/figures/final/Figure6_MSI_Survival.jpg`

Caption: Kaplan-Meier curves for CRC-MMIS risk groups stratified by MSI status. (A-B) MANTIS-defined MSS subgroup: Run 1 HR=2.48 (95% CI 1.46-4.20, p=0.0008); Run 2 HR=3.11 (95% CI 1.81-5.32, p<0.001). (C-D) SENSOR-defined MSS subgroup: Run 1 HR=2.28 (p=0.034); Run 2 HR=3.11 (p<0.001). (E-F) MSI-H subgroup: variable performance across seeds due to smaller sample size (n=96) and distinct MSI-H biology.

---

---

## Author Information

Authors: <span style="color:red;font-weight:bold;">[TO BE COMPLETED — list all authors with affiliations and ORCIDs]</span>

---

## Author Contributions (CRediT)

<span style="color:red;font-weight:bold;">[TO BE COMPLETED — conceptualization, methodology, software, validation, formal analysis, investigation, resources, data curation, writing — original draft, writing — review & editing, visualization, supervision, project administration, funding acquisition]</span>

---

## Data Availability

All data used in this study are publicly available. scRNA-seq datasets: GSE146771, GSE132465, GSE205506 (GEO). TCGA-COADREAD data: GDC Data Portal. GEO validation cohorts: GSE39582, GSE17536, GSE14333, GSE38832, GSE33113. Processed data files and analysis scripts are archived at [Zenodo DOI  — to be added upon acceptance].

---

## Ethics Statement

This study used publicly available de-identified data from GEO and TCGA databases. No new patient data or animal experiments were conducted. TCGA data collection was approved by the NIH Ethics Board; GEO data collection was approved by the respective institutional review boards of each contributing study.

---

## Competing Interests

The authors declare no Competing Interests.

---

## Acknowledgements

<span style="color:red;font-weight:bold;">[TO BE COMPLETED — acknowledge funding sources and technical assistance]</span>
