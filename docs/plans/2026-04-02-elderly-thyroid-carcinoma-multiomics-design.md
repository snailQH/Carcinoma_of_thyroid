# Elderly Thyroid Carcinoma Multi-Omics Molecular Subtyping — Research Plan

**Project:** Molecular subtyping of elderly thyroid carcinoma for clinical decision support
**Date:** 2026-04-02
**Status:** Phase 1-7 partially complete, external validation done, elderly vs young done
**Last updated:** 2026-04-03

---

## 0. Analysis Progress Tracker

### Phase 1: Data Acquisition & Preprocessing
- [x] 1. Download TCGA-THCA RNA-seq via TCGAbiolinks/GDC (572 samples, counts + TPM)
- [x] 2. Download TCGA-THCA somatic mutations (MAF, 5,834 variants)
- [x] 3. Download TCGA-THCA clinical data (505 patients)
- [ ] 4. Download TCGA-THCA CNV data (GISTIC2.0)
- [ ] 5. Download TCGA-THCA DNA methylation (450K)
- [ ] 6. Download TCGA-THCA miRNA-seq
- [x] 7. Download GEO scRNA-seq GSE184362 (23 samples, 11 PTC patients)
- [ ] 8. Download GEO scRNA-seq GSE191288
- [ ] 9. Download drug response databases (GDSC2, PRISM, CMap)
- [x] 10. RNA-seq preprocessing: TPM normalization, log2 transform, gene filtering (16,375 genes)
- [x] 11. Gene ID conversion: ENSEMBL -> gene symbols
- [x] 12. Construct elderly cohort (>=60): **120 samples**
- [x] 13. Demographics table: elderly vs young comparison (Table 1)

### Phase 2: Molecular Subtyping (Aim 1)
- [x] 14. Feature selection: MAD top 3000 genes
- [x] 15. Consensus Clustering (k=2-6, 1000 reps, Pearson + ward.D2)
- [ ] 16. NMF complementary clustering
- [ ] 17. iCluster+ / MOVICS multi-omics integration
- [x] 18. Optimal k determination: **k=2** (PAC=0.084, silhouette=0.95)
- [x] 19. Subtype assignment: **C1 (n=41), C2 (n=79)**
- [x] 20. Subtype heatmap (Fig 1)
- [ ] 21. Sensitivity analysis with age >=50 cutoff (n=224)
- [ ] 22. External cohort validation (NTP on GEO datasets)

### Phase 3: Subtype Characterization (Aim 2)
- [x] 23. Driver mutation frequency per subtype (BRAF, RAS, RET, TP53, etc.)
- [x] 24. TMB (tumor mutational burden) per subtype (Fig 2)
- [ ] 25. Co-mutation pattern analysis (Fisher's exact pairwise)
- [ ] 26. OncoPrint visualization (maftools)
- [ ] 27. CNA analysis (GISTIC2.0) -- pending CNV data download
- [x] 28. GSVA pathway scoring: 11 gene sets (EMT, Proliferation, Angiogenesis, Adhesion, Glycolysis, OXPHOS, SASP, Senescence, Immunosenescence, DNA Repair, Apoptosis) (Fig 4)
- [x] 29. Thyroid Differentiation Score (TDS): 10 genes scored (Fig 4)
- [x] 30. GSVA boxplots per pathway per subtype (11 plots)
- [x] 31. Immune marker-based scoring: CD8+ T, Treg, M1/M2 macrophage, NK, B cell, CAF, MDSC
- [x] 32. Immune landscape heatmap (Fig 3)
- [x] 33. Immune cell boxplots per subtype (6 plots)
- [x] 34. Checkpoint expression heatmap: CD274, PDCD1, CTLA4, LAG3, HAVCR2, TIGIT, CD276 (Fig 3)
- [ ] 35. CIBERSORTx deconvolution (LM22)
- [ ] 36. xCell / EPIC / MCPcounter deconvolution
- [ ] 37. ESTIMATE immune/stromal scores
- [ ] 38. Immune phenotype classification (hot/cold/suppressed/excluded)
- [ ] 39. TLS (Tertiary Lymphoid Structure) signature scoring
- [ ] 40. TIDE score computation
- [ ] 41. IPS (immunophenoscore)
- [ ] 42. Metabolic pathway deep analysis (fatty acid, amino acid, iodine, one-carbon)
- [x] 43. Clinical association: AJCC stage by subtype (Fisher p=0.0005)
- [x] 44. LN metastasis by subtype: C1=2.4%, C2=44.3% (Fig 6)
- [x] 45. Kaplan-Meier survival analysis (Fig 6)
- [x] 46. Clinical summary table per subtype
- [ ] 47. Multivariate logistic regression
- [ ] 48. Forest plot visualization

### Phase 4: Single-Cell Analysis (Aim 3)
- [x] 49. GSE184362 data loading: 23 samples (Tumor/Paratumor/LN) in 10x format
- [x] 50. scRNA-seq QC: cell/gene filtering, MT% filtering
- [x] 51. Normalization, HVG selection (3000 genes), PCA
- [x] 52. Harmony batch correction across samples
- [x] 53. UMAP + Leiden clustering (Fig 6)
- [x] 54. Cell type marker dotplot
- [x] 55. Cell type proportion by tissue (Fig 6)
- [x] 56. Tumor vs Paratumor DEG analysis
- [x] 57. Immune checkpoint expression on UMAP (Fig 6)
- [x] 58. Cluster marker gene identification
- [x] 59. Processed h5ad saved (5.3 GB)
- [ ] 60. CellTypist automated annotation (attempted, needs model download)
- [ ] 61. InferCNV / CopyKAT malignant cell identification
- [ ] 62. Tumor cell sub-clustering + state scoring
- [ ] 63. Pseudotime / trajectory analysis (scVelo / Monocle3)
- [ ] 64. CellChat / NicheNet cell-cell communication
- [ ] 65. pySCENIC regulon / TF analysis
- [ ] 66. Map bulk subtypes to single-cell states

### Phase 5: Drug Sensitivity & Targets (Aim 4)
- [x] 67. Drug target expression scoring: 8 drug classes (BRAF-i, MEK-i, lenvatinib, sorafenib, CDK-i, mTOR-i, anti-PD1, RAI)
- [x] 68. Drug sensitivity heatmap (subtypes x drugs) (Fig 7)
- [x] 69. Drug sensitivity boxplots (8 plots)
- [x] 70. TIDE-like dysfunction vs exclusion scatter (Fig 7)
- [x] 71. Actionable target expression table (17 genes)
- [x] 72. Actionable target heatmap (Fig 7)
- [ ] 73. oncoPredict / pRRophetic with GDSC2 training data
- [ ] 74. CMap / L1000 reverse signature query
- [ ] 75. DGIdb / OncoKB cross-reference
- [ ] 76. PPI network of subtype hub genes
- [ ] 77. IPS (immunophenoscore) from TCIA
- [ ] 78. SubMap immunotherapy responder comparison
- [ ] 79. Neoantigen load estimation

### Phase 6: Clinical Model & Validation (Aim 5)
- [x] 80. ANOVA feature selection (100 top discriminating genes)
- [x] 81. Logistic Regression classifier: 97.0% balanced accuracy
- [x] 82. Random Forest classifier: **97.6% balanced accuracy** (best)
- [x] 83. Gradient Boosting classifier: 97.0% balanced accuracy
- [x] 84. Feature importance ranking (top 20 genes) (Fig 8)
- [x] 85. Confusion matrix (5-fold CV) (Fig 8)
- [x] 86. Clinical decision framework table
- [ ] 87. LASSO / elastic-net gene signature (minimal panel)
- [ ] 88. Risk score construction
- [x] 89. External validation on GEO bulk cohorts (GSE29265 n=49, GSE33630 n=105, NTP scoring)
- [ ] 90. Time-dependent ROC / C-index
- [ ] 91. Calibration plot
- [ ] 92. R Shiny / Python scoring tool

### Phase 6A: Aging/Senescence Deep Analysis (NEW — 2026-04-05)
- [x] 92aa. Build 7 senescence gene signatures (SenMayo, SASP Atlas, Cell Cycle Arrest, Telomere, DDR, Immunosenescence, Anti-apoptosis)
- [x] 92ab. GSVA scoring on all 505 tumor samples
- [x] 92ac. Senescence vs age correlation (Spearman, all patients)
- [x] 92ad. Elderly vs young senescence comparison (7 boxplots)
- [x] 92ae. Senescence by AJCC stage (composite score)
- [x] 92af. Senescence by LN metastasis (N0 vs N1)
- [x] 92ag. **Multivariate logistic regression: Senescence OR=4.374 (p=8.6e-5) strongest predictor of LN met**
- [x] 92ah. ROC analysis: Senescence AUC=0.636, Combined AUC=0.672
- [x] 92ai. Senescence landscape heatmap (all samples, 7 signatures)
- [x] 92aj. Senescence by elderly subtype (C1 vs C2, 7 signatures + composite)
- [x] 92ak. Key senescence gene expression heatmap (20 genes)
- [x] 92al. Survival by senescence group (all patients + elderly)
- [x] 92am. Senescence scatter plots vs age (7 signatures)
- [x] 92an. Forest plot for multivariate predictors

- [x] 92ao. C1/C2 specificity analysis: BRAF/RAS effect amplification in elderly (RAS-like elderly = 0% LN met)
- [x] 92ao2. Transcriptomic age predictor (14 up + 12 down aging genes)
- [x] 92ap. Age acceleration: C2 = +0.280, C1 = -0.461, p = 4.5e-12
- [x] 92aq. Age acceleration vs senescence correlation (rho=0.28, p=1.8e-10)
- [x] 92ar. Age acceleration by LN status and AJCC stage
- [x] 92as. Aging score scatter plots (6 figures)

### Phase 6B: Elderly vs Young Adult Comparison (NEW)
- [x] 92a. Elderly vs young demographics comparison (120 vs 385 patients)
- [x] 92b. DEG analysis: elderly vs young (volcano plot, t-test per gene)
- [x] 92c. GSVA pathway comparison: elderly vs young (11 pathways, boxplots)
- [x] 92d. Immune cell scoring comparison: elderly vs young (6 cell types)
- [x] 92e. Mutation frequency comparison: elderly vs young (10 drivers, BRAF Fisher p=0.83)
- [x] 92f. Clinical comparison: LN metastasis (30% vs 41.6%), Stage III/IV (51.7% vs 22.3%)
- [x] 92g. Survival comparison: elderly vs young (KM curve)
- [x] 92h. TDS comparison: elderly vs young

### Phase 7: Outputs & Manuscript
- [x] 93. Summary PPTX generated (17 slides with speaker notes, figures, tables embedded)
- [x] 94. All main figures (Fig 1-8) generated as PDF
- [x] 95. All supplementary boxplots generated (30+ plots)
- [x] 96. All summary tables generated (11 CSV files)
- [ ] 97. Publication-quality figure polishing (multi-panel composites)
- [ ] 98. Supplementary figures compilation
- [ ] 99. Methods section draft
- [ ] 100. Results section draft
- [ ] 101. Manuscript assembly

### Additional Novel Analyses (Section 11 — Future)
- [ ] 102. Epigenetic clock analysis (Horvath, needs methylation data)
- [ ] 103. Alternative splicing (rMATS/SUPPA2, needs BAM files)
- [ ] 104. Computational spatial inference (SpaCET/CytoSPACE)
- [ ] 105. Pan-cancer elderly comparison
- [ ] 106. Circulating biomarker candidates (secretome filtering)
- [ ] 107. RAI resistance prediction (TDS + NIS deep analysis)
- [ ] 108. Immune Age Score construction
- [ ] 109. Cancer stemness scoring (DepMap/CRISPR)
- [ ] 110. Deep learning drug sensitivity (PASO)
- [ ] 111. Mendelian randomization for causal targets
- [ ] 112. Immunosenescence plasma biomarker cross-reference

- [x] 97a. Manuscript draft v1 with embedded figures and elderly vs young section
- [x] 97b. External validation section added to manuscript

**Summary: 84 of 140 analyses completed (60%)**

---

## 1. Scientific Question

> Does elderly thyroid carcinoma (age >=60) harbor distinct molecular subtypes that differ in aggressiveness, immune microenvironment, metabolic rewiring, and drug sensitivity — and can these subtypes guide personalized clinical decision-making (surgery vs. active surveillance vs. targeted therapy)?

---

## 1.1 Literature Landscape & Novelty Justification (as of April 2026)

**Key finding: No published study combines elderly-specific molecular subtyping with clinical decision modeling in thyroid cancer.** This is our unique niche.

### Recent relevant work to position against:

| Study | Year | Key contribution | How we differ |
|-------|------|-----------------|---------------|
| Proteogenomic subtypes of advanced DTC (Cell Reports Medicine) | 2026 | 3 subtypes (canonical, stromal, immunogenic) from 113 advanced DTC via proteogenomics | We focus on elderly-specific PTC, use TCGA (larger n), and add clinical decision framework |
| Multi-omics consensus clustering (Genes & Immunity) | 2025 | 2 prognostic subtypes (CS1/CS2) from 539 pts using DNA methylation + mRNA + lncRNA + miRNA + mutations | We add elderly-specific angle, senescence biology, drug sensitivity prediction |
| Cell senescence-related signature in PTC (Aging) | 2024 | 4 senescence markers (SNAI1, CDKN2A, HDAC4, NDRG1) for prognosis | We integrate senescence into multi-dimensional subtyping, not just prognostic genes |
| E2F1 as senescence-linked prognostic gene (Frontiers in Genetics) | 2025 | Single gene focus | We build comprehensive senescence + immune age scoring |
| Aging drives aggressiveness via 1p/1q alterations (Aging and Disease) | 2023 | Age 55+ independent of BRAF drives aggressiveness through chromosomal changes, depleted CD8+ T cells | Most directly relevant; we extend with subtyping, drug sensitivity, clinical framework |
| Immunosenescence plasma protein biomarkers (Frontiers in Oncology) | 2025 | 6 plasma proteins for thyroid cancer risk via GWAS + proteomics + transcriptomics | Complementary; we could cross-reference our subtype markers with their biomarkers |
| Protein-based XGBoost classifier (EMBO Mol Med) | 2025 | 24-protein classifier (AUROC 0.899) for follicular adenoma vs carcinoma | Different question (diagnosis vs subtyping); our ML classifier is for subtype assignment |
| scRNA-seq TIME atlas — 29 subpopulations in PTC (Front Immunol) | 2025 | High-resolution cellular atlas | We leverage similar data but link to elderly subtypes + drug sensitivity |

### Emerging methodologies to consider:

- **Spatial transcriptomics** for TIME analysis is trending (Cancers, 2025) — our computational spatial inference (SpaCET/CytoSPACE) aligns with this
- **Emerging immune checkpoints**: LAG-3 and TIM-3 alongside PD-1 for exhaustion profiling — include in our checkpoint analysis
- **Deep learning drug sensitivity** (PASO model: transformer + multi-scale CNNs) — could supplement oncoPredict as a modern alternative
- **Mendelian randomization** for drug target prioritization (Endocrine, 2024) — potential add-on for causal inference
- **Cancer stemness scoring** from CRISPR data — novel angle for subtype characterization

### Our competitive advantages:

1. **Elderly-specific focus** — underexplored niche with clinical relevance (surgical risk in elderly)
2. **Multi-dimensional subtyping** — mutation + immune + senescence + metabolism + drug sensitivity (most studies cover 1-2 dimensions)
3. **Clinical decision orientation** — not just subtypes, but actionable management guidance
4. **Drug sensitivity integration** — connects subtypes to therapeutic opportunities
5. **Senescence + immune age** — novel angle bridging aging biology and cancer immunology

---

## 2. Study Design Overview

```
TCGA-THCA (bulk multi-omics, n~500)
    |
    v
Filter elderly cohort (>=60, sensitivity: >=50)
    |
    v
[Aim 1] Unsupervised molecular subtyping (RNA-seq + optional multi-omics integration)
    |
    +--> [Aim 2] Subtype characterization
    |       - Mutation landscape (BRAF/RAS/TERT/TMB)
    |       - Tumor biology (GSVA: EMT, proliferation, adhesion, angiogenesis)
    |       - Immune microenvironment (CIBERSORTx, xCell, EPIC)
    |       - Aging/senescence signatures
    |       - Metabolic reprogramming
    |       - Clinical association (stage, LN metastasis, ETE, recurrence)
    |
    +--> [Aim 3] Single-cell mechanistic validation (GEO scRNA-seq)
    |       - Cell-type deconvolution & tumor cell states
    |       - Cell-cell communication (CellChat / NicheNet)
    |       - Trajectory / pseudotime analysis
    |
    +--> [Aim 4] Drug sensitivity & therapeutic target prediction
    |       - CMap / L1000 connectivity mapping
    |       - GDSC / PRISM drug response prediction
    |       - Immune checkpoint expression profiling
    |       - Actionable target identification
    |
    +--> [Aim 5] Clinical decision model & validation
            - LASSO / elastic-net gene signature
            - Risk score construction
            - External validation (GEO bulk cohorts)
            - ML classifier for subtype assignment
            - Decision framework: subtype -> clinical recommendation
```

---

## 3. Data Sources

### 3.1 Primary: TCGA-THCA (Bulk Multi-Omics)

| Data type | Source | Access |
|-----------|--------|--------|
| RNA-seq (HTSeq counts + TPM) | GDC | Public |
| Somatic mutations (MAF) | GDC | Public |
| Copy number variation (SNP6) | GDC | Public |
| DNA methylation (450K) | GDC | Public |
| miRNA-seq | GDC | Public |
| Clinical & follow-up | GDC / TCGA-CDR | Public |

**Download tools:** `TCGAbiolinks` (R) or `gdc-client` (CLI)

### 3.2 Single-Cell RNA-seq (Mechanism Validation)

| Dataset | Description | Samples |
|---------|-------------|---------|
| GSE184362 | PTC: tumor + LN metastasis + normal | Multi-tissue |
| GSE191288 | Thyroid cancer immune microenvironment | Tumor + normal |
| (search for more) | ATC/PDTC if available | Progression |

**Tools:** `scanpy` (Python) or `Seurat` (R)

### 3.3 External Validation (Bulk)

- GEO bulk RNA-seq/microarray datasets (e.g., GSE29265, GSE33630, GSE60542)
- UCSC Xena preprocessed matrices (quick access)

### 3.4 Drug Response Databases

| Database | Content | Use |
|----------|---------|-----|
| GDSC (Sanger) | Drug IC50 across cell lines | pRRophetic / oncoPredict |
| PRISM (Broad) | Drug sensitivity, 4,518 compounds | Broader drug screening |
| CMap / L1000 | Gene expression connectivity | Reverse signature matching |
| DrugBank | Drug-target mapping | Actionable target annotation |
| ClinicalTrials.gov | Ongoing trials | Translational relevance |

### 3.5 Supplementary

- **cBioPortal**: Oncoprint, co-mutation visualization
- **Human Protein Atlas**: Protein-level validation of key biomarkers
- **TISCH**: Tumor immune single-cell hub (pre-annotated)
- **STRING**: Protein-protein interaction networks

---

## 4. Detailed Aim Design

### Aim 1: Molecular Subtyping of Elderly Thyroid Carcinoma

**Objective:** Identify reproducible molecular subtypes within elderly (>=60) PTC patients.

**Methods:**

1. **Cohort construction**
   - Filter TCGA-THCA for age >= 60 at diagnosis
   - Record sample size; if < 80, expand to >= 50 with sensitivity analysis
   - Compare elderly vs. young demographic/clinical characteristics (Table 1)

2. **Gene expression preprocessing**
   - RNA-seq counts -> TPM -> log2(TPM+1)
   - Remove low-expression genes (median TPM < 1)
   - Batch correction if merging with external data (ComBat-seq)

3. **Feature selection**
   - Median absolute deviation (MAD) top 3000-5000 genes
   - Alternatively: highly variable genes via `scanpy`-style dispersion

4. **Clustering approaches** (run all, compare)
   - **Consensus Clustering** (ConsensusClusterPlus, R) — primary
   - **NMF** (non-negative matrix factorization) — complementary
   - **iCluster+** (if integrating mutation + methylation + expression)
   - **MOVICS** (multi-omics integration framework, recommended)
   - Determine optimal k via CDF, PAC, silhouette, gap statistic

5. **Subtype validation**
   - Bootstrap resampling (1000x) stability
   - Silhouette width per sample
   - External cohort reproducibility (NTP — nearest template prediction)

**Expected output:** 2-4 stable molecular subtypes with distinct expression profiles.

---

### Aim 2: Multi-Dimensional Subtype Characterization

**Objective:** Define the biological identity of each subtype across mutation, immune, metabolic, and clinical dimensions.

#### 2A. Mutation Landscape

- **Driver mutations**: BRAF V600E, RAS family (NRAS/HRAS/KRAS), TERT promoter, RET fusions, PAX8-PPARG
- **Tumor mutational burden (TMB)**: mutations/Mb per subtype
- **Co-mutation patterns**: Fisher's exact test for pairwise enrichment
- **Visualization**: OncoPrint (ComplexHeatmap / maftools), lollipop plots

#### 2B. Copy Number Alterations

- GISTIC2.0 for recurrent CNA identification per subtype
- Arm-level vs. focal events
- Key loci: 1q gain, 22q loss (known in thyroid cancer)

#### 2C. Tumor Biology (Pathway Activity)

**Method:** GSVA / ssGSEA with MSigDB hallmark + curated gene sets

| Pathway category | Key gene sets |
|-----------------|---------------|
| Proliferation | MKI67 signature, E2F targets, G2M checkpoint |
| EMT | Hallmark EMT, Taube EMT signature |
| Cell adhesion | Integrin signaling, focal adhesion, ECM-receptor |
| Angiogenesis | VEGF pathway, hallmark angiogenesis |
| DNA damage repair | Hallmark DNA repair, HR/NHEJ signatures |
| Apoptosis | Hallmark apoptosis, BCL2 family |
| Thyroid differentiation | TDS (thyroid differentiation score): TG, TPO, NIS, TSHR, PAX8, FOXE1, NKX2-1 |

**NEW — Thyroid Differentiation Score (TDS):**
A well-established metric in thyroid cancer (Cancer Genome Atlas Research Network, 2014). Calculate per sample and compare across subtypes — this directly links to de-differentiation and aggressiveness.

#### 2D. Aging & Senescence Signatures (Novel angle)

**Rationale:** Elderly-specific study — aging biology is directly relevant and rarely explored in thyroid cancer subtyping.

| Signature | Source | Genes/method |
|-----------|--------|--------------|
| Cellular senescence | SenMayo gene set (Saul et al. 2022) | 125 genes |
| SASP (senescence-associated secretory phenotype) | Coppe et al. | IL6, IL8, CXCL1, MMP3, etc. |
| Telomere maintenance | TERT expression + ALT markers | TERT, ATRX, DAXX |
| Oxidative stress | Hallmark ROS pathway | |
| Epigenetic aging | Horvath / Hannum clock (if methylation used) | CpG-based |
| Immunosenescence | T cell exhaustion + senescent T markers | KLRG1, CD57, TIGIT |

**Analysis:**
- Score each signature per sample (ssGSEA)
- Correlate with subtypes
- Test: do elderly subtypes differ in senescence burden?
- Hypothesis: one subtype may be "senescence-driven" with paradoxically better prognosis (senescent cells = less proliferative)

#### 2E. Immune Microenvironment (Deep Profiling)

**Deconvolution methods** (use >= 2 for robustness):
- CIBERSORTx (LM22 signature)
- xCell (64 cell types)
- EPIC
- MCPcounter
- ESTIMATE (immune/stromal scores)

**Immune phenotype classification:**
- Immune-hot (high CD8+, high IFN-gamma)
- Immune-cold (low infiltration overall)
- Immune-suppressed (high Treg, M2 macrophage, MDSC)
- Immune-excluded (high stromal, low intratumoral immune)

**Checkpoint & immunotherapy markers:**
- PD-L1 (CD274), PD-1 (PDCD1), CTLA-4, LAG3, TIM-3 (HAVCR2), TIGIT, VISTA, B7-H3
- TIDE score (tumor immune dysfunction and exclusion)
- IPS (immunophenoscore)
- MSI status (if applicable)

**NEW — Tertiary Lymphoid Structure (TLS) signature:**
TLS are increasingly recognized as predictive of immunotherapy response. Score TLS using established gene signatures (e.g., Cabrita et al. 2020, Helmink et al. 2020).

#### 2F. Metabolic Reprogramming (Novel angle)

**Rationale:** Thyroid is a metabolically active organ; metabolic shifts in elderly tumors are understudied.

| Pathway | Gene set source |
|---------|----------------|
| Glycolysis | Hallmark glycolysis |
| OXPHOS | Hallmark oxidative phosphorylation |
| Fatty acid metabolism | Hallmark fatty acid metabolism |
| Amino acid metabolism | KEGG arginine/proline, tryptophan |
| Iodine metabolism | NIS (SLC5A5), TPO, TG, DIO1/2 |
| One-carbon metabolism | MTHFR, SHMT1/2, TYMS |

**Analysis:**
- ssGSEA scoring per pathway
- Compare across subtypes
- Correlate with thyroid differentiation score
- Hypothesis: de-differentiated subtypes lose iodine metabolism -> RAI resistance

#### 2G. Clinical Association

**Primary endpoints** (surrogate for survival):
- AJCC stage (I-IV)
- Lymph node metastasis (N0 vs N1)
- Extrathyroidal extension (ETE)
- Tumor size
- Multifocality
- Recurrence / disease-free interval (if available)
- BRAF V600E status

**Secondary** (if data sufficient):
- Overall survival (Kaplan-Meier + log-rank, acknowledging limited events)
- Disease-specific survival

**Statistical tests:**
- Chi-square / Fisher's exact for categorical
- Kruskal-Wallis / ANOVA for continuous
- Multivariate logistic regression for independent predictors
- Forest plot visualization

---

### Aim 3: Single-Cell Mechanistic Validation

**Objective:** Resolve subtype-driving mechanisms at single-cell resolution.

#### 3A. Data Processing

- Standard scRNA-seq pipeline: QC -> normalization -> HVG -> PCA -> UMAP -> clustering
- Tools: `scanpy` or `Seurat`
- Doublet removal: `scrublet` or `DoubletFinder`
- Cell type annotation: `celltypist` or manual markers + reference mapping

#### 3B. Cell Type Composition

- Compare proportions across tumor vs. normal vs. metastasis
- Key populations:
  - Thyrocytes (tumor cells): differentiated vs. de-differentiated states
  - T cells: CD8+ effector, exhausted (TOX+, PDCD1+), Treg (FOXP3+)
  - Macrophages: M1 vs. M2 polarization
  - CAFs (cancer-associated fibroblasts): inflammatory vs. myofibroblastic
  - Endothelial cells: tip cells, lymphatic endothelium
  - DCs: cDC1 (cross-presentation), cDC2, pDC

#### 3C. Tumor Cell Heterogeneity

- **InferCNV** / **CopyKAT**: distinguish malignant from non-malignant cells
- Sub-cluster tumor cells -> identify transcriptional states
- Score for: EMT, differentiation, proliferation, stemness
- **Pseudotime / trajectory** (Monocle3 or scVelo): model de-differentiation trajectory
- Map bulk subtypes onto single-cell states using gene signatures (NTP or correlation)

#### 3D. Cell-Cell Communication

- **CellChat**: ligand-receptor interaction network
- **NicheNet**: predict which ligands drive target gene expression in receiver cells
- Focus on:
  - Tumor <-> immune suppressive signals
  - CAF <-> tumor (ECM remodeling, growth factors)
  - Endothelial <-> tumor (adhesion, intravasation)

#### 3E. Transcription Factor Regulon Analysis (NEW)

- **SCENIC** (pySCENIC): infer active transcription factor regulons
- Identify master regulators per subtype
- Known thyroid TFs: PAX8, NKX2-1 (TTF1), FOXE1
- Oncogenic TFs: MYC, TWIST1, SNAI1/2, ZEB1 (EMT drivers)

---

### Aim 4: Drug Sensitivity & Therapeutic Target Prediction (NEW)

**Objective:** Predict drug responses per subtype and identify actionable therapeutic targets.

#### 4A. Drug Sensitivity Prediction from Bulk Expression

**Method 1 — oncoPredict / pRRophetic:**
- Train ridge regression on GDSC2 / PRISM cell line data
- Predict IC50 per TCGA sample
- Compare predicted drug sensitivity across subtypes
- Key drugs to examine:
  - **BRAF inhibitors**: vemurafenib, dabrafenib
  - **MEK inhibitors**: trametinib, selumetinib
  - **Multi-kinase inhibitors**: sorafenib, lenvatinib (FDA-approved for thyroid)
  - **Immune checkpoint inhibitors**: anti-PD1 response prediction (TIDE)
  - **CDK inhibitors**: palbociclib (if proliferative subtype)
  - **mTOR inhibitors**: everolimus
  - **Ferroptosis inducers**: erastin, RSL3 (if metabolic subtype)

**Method 2 — CMap / L1000 Connectivity Mapping:**
- Extract subtype-specific DEG signatures (up + down)
- Query CMap/L1000 for compounds that reverse the disease signature
- Rank candidate drugs by connectivity score
- Cross-reference with DrugBank for clinical feasibility

#### 4B. Actionable Target Identification

- Cross-reference subtype driver genes with:
  - **DGIdb** (Drug-Gene Interaction database)
  - **OncoKB** (actionable alterations)
  - **CIViC** (Clinical Interpretation of Variants in Cancer)
- Identify subtype-specific druggable targets
- Network analysis: PPI network of subtype hub genes -> druggable nodes

#### 4C. Immune Checkpoint & Immunotherapy Prediction

- TIDE score per subtype (predict anti-PD1/anti-CTLA4 response)
- IPS (immunophenoscore) from TCIA
- SubMap analysis: compare subtypes to immunotherapy responder/non-responder profiles
- Neoantigen load estimation (from mutation data)

#### 4D. Visualization

- Drug sensitivity heatmap (subtypes x drugs)
- Violin plots of predicted IC50 per subtype
- Network diagrams of druggable targets
- Sankey diagram: subtype -> drug class -> specific agents

---

### Aim 5: Clinical Decision Model & Validation

**Objective:** Build a practical risk stratification tool and decision framework.

#### 5A. Gene Signature Development

- LASSO / elastic-net Cox regression (or logistic for surrogate endpoints)
- Input: subtype-defining DEGs
- Output: minimal gene panel (target: 10-30 genes)
- Risk score = weighted sum of gene expression

#### 5B. Risk Score Validation

- Internal: LOOCV / 10-fold CV in TCGA
- External: apply signature to GEO validation cohorts
- Calibration plot + C-index / AUC
- Time-dependent ROC (if survival used)

#### 5C. Machine Learning Classifier

- Train RF / SVM / XGBoost to assign new patients to subtypes
- Feature importance ranking
- Evaluate: accuracy, kappa, balanced accuracy
- Deploy as a simple scoring tool (R Shiny or Python script)

#### 5D. Clinical Decision Framework (Core Innovation)

| Subtype | Molecular features | Clinical profile | Suggested management |
|---------|--------------------|-----------------|---------------------|
| A (Indolent) | High TDS, low TMB, immune-hot | Low stage, no LN met | Active surveillance, follow-up |
| B (Immune-suppressed) | High Treg/M2, SASP+, senescence | Intermediate risk | Consider immunotherapy, close monitoring |
| C (Aggressive) | High EMT, adhesion, de-differentiated | High stage, LN met+, ETE | Surgical intervention + adjuvant therapy |
| D (Drug-sensitive)* | Specific pathway dependency | Variable | Targeted therapy (BRAF/MEK/lenvatinib) |

*Number of subtypes to be determined by data; this is illustrative.

**Important framing:** "These subtypes suggest potential stratification for personalized clinical decision-making in elderly thyroid carcinoma patients" — not direct clinical recommendations.

---

## 5. Figure Layout (Publication-Ready)

| Figure | Content | Key message |
|--------|---------|-------------|
| Fig 1 | Study design flowchart + cohort demographics + clustering heatmap | Subtypes exist in elderly PTC |
| Fig 2 | Mutation landscape: OncoPrint + TMB + co-mutation | Subtypes have distinct genetic drivers |
| Fig 3 | Immune microenvironment: deconvolution + checkpoint + TLS + TIDE | Subtypes differ in immune context |
| Fig 4 | Tumor biology: GSVA heatmap + TDS + EMT + senescence | Biological programs define subtypes |
| Fig 5 | Metabolic reprogramming + iodine metabolism | Metabolic vulnerability per subtype |
| Fig 6 | Single-cell: UMAP + cell proportions + CellChat + SCENIC | Cellular mechanisms driving subtypes |
| Fig 7 | Drug sensitivity: predicted IC50 + CMap hits + druggable targets | Therapeutic opportunities per subtype |
| Fig 8 | Clinical decision framework: risk score + validation + decision tree | Translational value |

**Supplementary:** Table S1 (clinical characteristics), Fig S1-S3 (sensitivity analyses: age cutoffs, clustering stability, method comparisons)

---

## 6. Computational Pipeline & Tools

### 6.1 Environment

- **Server:** `bioinfo@192.168.100.127`, path: `/data2/projects/Carcinoma_of_thyroid/`
- **Languages:** R (primary for TCGA/immune), Python (scRNA-seq, ML)
- **Environment management:** conda / mamba

### 6.2 Key R Packages

| Task | Package |
|------|---------|
| TCGA data download | TCGAbiolinks |
| Consensus clustering | ConsensusClusterPlus |
| Multi-omics integration | MOVICS, iCluster+ |
| Mutation analysis | maftools |
| Pathway scoring | GSVA, clusterProfiler |
| Immune deconvolution | immunedeconv (wraps CIBERSORT, xCell, EPIC, MCP) |
| Survival analysis | survival, survminer |
| Visualization | ComplexHeatmap, ggplot2, ggpubr |
| Drug prediction | oncoPredict |

### 6.3 Key Python Packages

| Task | Package |
|------|---------|
| scRNA-seq pipeline | scanpy |
| Cell annotation | celltypist |
| CNV inference | infercnvpy / CopyKAT |
| Cell communication | cellchat (R) or commot (Python) |
| Regulon analysis | pyscenic |
| Trajectory | scvelo, monocle3 (R) |
| ML classifier | scikit-learn, xgboost |
| CMap querying | cmapPy |

### 6.4 Directory Structure (Proposed)

```
/data2/projects/Carcinoma_of_thyroid/
├── data/
│   ├── tcga_thca/          # TCGA downloads (GDC)
│   │   ├── rnaseq/
│   │   ├── mutation/
│   │   ├── cnv/
│   │   ├── methylation/
│   │   └── clinical/
│   ├── geo_scrna/           # scRNA-seq datasets
│   │   ├── GSE184362/
│   │   └── GSE191288/
│   ├── geo_bulk_validation/ # External bulk validation
│   └── drug_databases/      # GDSC, PRISM, CMap
├── scripts/
│   ├── 01_data_download/
│   ├── 02_preprocessing/
│   ├── 03_subtyping/
│   ├── 04_characterization/
│   ├── 05_scrna_analysis/
│   ├── 06_drug_prediction/
│   ├── 07_clinical_model/
│   └── utils/
├── results/
│   ├── figures/
│   ├── tables/
│   └── reports/
├── envs/                    # Conda environment files
└── docs/
```

---

## 7. Analysis Execution Order

> See **Section 0 (Analysis Progress Tracker)** at the top of this document for detailed checkbox status of all 112 sub-analyses.

---

## 8. Risk Mitigation

| Risk | Impact | Mitigation |
|------|--------|------------|
| Small elderly cohort (n<80) | Underpowered clustering | Expand to >=50; use age as continuous covariate in sensitivity |
| Unstable subtypes | Unreproducible results | Bootstrap 1000x; validate in external cohort; report PAC/silhouette |
| TCGA survival data too weak | Cannot show prognostic value | Use surrogate endpoints (stage, LN met, ETE); frame carefully |
| scRNA-seq age info missing | Cannot map elderly specifically | Use as general mechanism validation; discuss limitation |
| Drug prediction = computational only | Overinterpretation risk | Frame as "hypothesis-generating"; cross-validate with cell line data |
| Clinical recommendation overreach | Reviewer pushback | Use "potential stratification" language; no direct clinical advice |

---

## 9. Expected Outputs

1. **Molecular subtypes** of elderly thyroid carcinoma with biological annotations
2. **Comprehensive characterization** across mutation, immune, metabolic, aging dimensions
3. **Single-cell validation** of subtype-driving mechanisms
4. **Drug sensitivity profiles** and candidate therapeutic targets per subtype
5. **Gene signature + risk score** for subtype classification
6. **Clinical decision framework** linking subtypes to management strategies
7. **Publication-ready figures** (8 main + supplementary)
8. **Summary PPTX** with all key results

---

## 10. Target Journals

| Tier | Journal | IF range | Fit |
|------|---------|----------|-----|
| High | Clinical Cancer Research | ~12 | If external validation strong |
| High | iMeta | ~20 | Multi-omics integration focus |
| Mid-High | Briefings in Bioinformatics | ~9 | Methods + biology |
| Mid | Cancers | ~5 | Solid multi-omics cancer paper |
| Mid | Frontiers in Oncology | ~4 | Good backup |
| Mid | Frontiers in Immunology | ~7 | If immune story is primary |

---

## 11. Additional Novel Ideas (Beyond Original Plan)

### 11.1 Epigenetic Clock Analysis
If DNA methylation data is used: calculate biological age (Horvath clock) per sample. Test whether "accelerated aging" (biological age >> chronological age) correlates with aggressive subtypes. This is highly novel for thyroid cancer.

### 11.2 Alternative Splicing Events
Use TCGA RNA-seq BAM files to identify subtype-specific alternative splicing (rMATS or SUPPA2). Splicing dysregulation is emerging in cancer but barely explored in thyroid.

### 11.3 Tumor Microenvironment Spatial Context (Computational)
Use SpaCET or CytoSPACE to computationally infer spatial organization from bulk/scRNA data. Predict whether immune cells are infiltrating or excluded.

### 11.4 Pan-Cancer Elderly Comparison
Compare elderly thyroid subtypes with elderly subsets from other TCGA cancers (e.g., BRCA, LUAD). Ask: is the "indolent elderly subtype" a pan-cancer phenomenon or thyroid-specific?

### 11.5 Circulating Biomarker Candidates
From subtype-defining DEGs, filter for genes encoding secreted proteins (using Human Protein Atlas "secretome" annotation). These could serve as non-invasive blood biomarker candidates — highly relevant for elderly patients where minimizing procedures matters.

### 11.6 RAI (Radioactive Iodine) Resistance Prediction
Thyroid differentiation score + NIS expression → predict which subtypes would respond to RAI therapy vs. those that are RAI-refractory. Directly clinically actionable.

### 11.7 Immune Age Score
Construct a composite "immune age" score combining immunosenescence markers (exhausted T cells, reduced naive T cells, increased MDSC) and correlate with subtypes. Novel concept bridging aging immunology and cancer.

### 11.8 Cancer Stemness Scoring
Use CRISPR-based stemness indices (from DepMap/CRISPR screens) to score cancer stemness per sample. A 2025 study used stemness models across 13 transcriptomic datasets for risk stratification — this could add a differentiation hierarchy dimension to our subtypes.

### 11.9 Deep Learning Drug Sensitivity (PASO)
Supplement oncoPredict with the PASO model (transformer encoder + multi-scale CNNs), which predicts cell-line drug sensitivity from omics + molecular structures. More modern than ridge regression approaches.

### 11.10 Mendelian Randomization for Causal Target Validation
Use MR (two-sample, using GWAS summary stats + eQTL) to test whether subtype-defining genes have causal relationships with thyroid cancer risk. Adds causal inference rigor beyond correlative analyses.

### 11.11 Cross-Reference with Immunosenescence Plasma Biomarkers
A 2025 study (Frontiers in Oncology) identified 6 immunosenescence-associated plasma protein biomarkers for thyroid cancer risk. Cross-reference our subtype DEGs with these biomarkers to validate and extend their findings in the elderly-specific context.

---

## 12. Key References

1. Proteogenomic subtypes of advanced DTC — *Cell Reports Medicine*, 2026 (DOI: 10.1016/j.xcrm.2026.00078)
2. Multi-omics consensus clustering in thyroid cancer — *Genes & Immunity*, 2025 (DOI: 10.1038/s41435-025-00322-w)
3. Cell senescence-related signature (CSRS) in PTC — *Aging*, 2024 (DOI: 10.18632/aging.205520)
4. E2F1 senescence gene in PTC — *Frontiers in Genetics*, 2025 (DOI: 10.3389/fgene.2025.1605385)
5. Aging drives aggressiveness independent of BRAF — *Aging and Disease*, 2023 (PMC10187705)
6. Immunosenescence plasma biomarkers — *Frontiers in Oncology*, 2025 (DOI: 10.3389/fonc.2024.1525767)
7. Protein-based XGBoost classifier — *EMBO Molecular Medicine*, 2025 (DOI: 10.1038/s44321-025-00242-2)
8. scRNA-seq thyroid TIME atlas (29 subpopulations) — *Frontiers in Immunology*, 2025 (DOI: 10.3389/fimmu.2025.1738583)
9. Spatial transcriptomics in thyroid TIME — *Cancers*, 2025 (DOI: 10.3390/cancers17050794)
10. MTC multi-omics + integrative ML — *Nature Communications*, 2026 (DOI: 10.1038/s41467-025-67533-7)
11. Drug target prioritization via Mendelian randomization — *Endocrine*, 2024 (DOI: 10.1007/s12020-024-03933-x)
12. Cancer stemness model — *Discover Oncology*, 2025 (DOI: 10.1007/s12672-025-02813-8)
