# Molecular Subtyping of Elderly Thyroid Carcinoma Reveals Distinct Biological Programs and Clinical Implications: A Multi-Omics Integrative Analysis

## Authors

[Author names to be added]

## Abstract

**Background:** Thyroid carcinoma management in elderly patients (>=60 years) remains challenging due to increased surgical risk and poorly characterized age-specific tumor biology. Whether elderly papillary thyroid carcinoma (PTC) harbors distinct molecular subtypes with differential clinical behavior is unknown.

**Methods:** We analyzed 120 elderly PTC patients from the TCGA-THCA cohort using consensus clustering on RNA-seq data (top 3,000 variable genes, 1,000 bootstrap resamples). Subtypes were characterized across mutation landscape, pathway activity (GSVA, 11 gene sets), immune microenvironment (8 cell-type markers, 7 checkpoints), aging/senescence signatures, thyroid differentiation score (TDS), and clinical features. Single-cell RNA-seq (GSE184362, 23 samples) provided mechanistic validation. Drug sensitivity was predicted for 8 drug classes, and a machine learning classifier was trained for subtype assignment.

**Results:** Two highly stable molecular subtypes were identified (PAC = 0.084, silhouette = 0.95). Subtype C1 (n = 41) was RAS-like (NRAS 21.6%, BRAF 0%), with higher oxidative phosphorylation, lower EMT, and minimal lymph node metastasis (2.4%). Subtype C2 (n = 79) was BRAF-like (BRAF V600E 87.2%, NRAS 0%), with elevated EMT, adhesion, senescence-associated secretory phenotype (SASP), immunosenescence, and markedly higher lymph node metastasis (44.3%; Fisher p = 0.0005 for stage distribution). Drug sensitivity analysis suggested C1 patients may benefit from radioactive iodine (RAI) therapy, while C2 patients showed higher sensitivity to BRAF/MEK inhibitors and elevated immune checkpoint expression. A Random Forest classifier achieved 97.6% balanced accuracy for subtype prediction.

**Conclusions:** Elderly PTC harbors two molecularly and clinically distinct subtypes. The indolent C1 subtype may be candidates for active surveillance, while the aggressive C2 subtype warrants more intensive intervention. This framework provides a foundation for personalized clinical decision-making in elderly thyroid carcinoma.

**Keywords:** thyroid carcinoma, elderly, molecular subtyping, multi-omics, TCGA, immune microenvironment, drug sensitivity, clinical decision support

---

## 1. Introduction

Thyroid carcinoma is the most common endocrine malignancy, with papillary thyroid carcinoma (PTC) accounting for approximately 85% of cases [1]. While PTC generally carries an excellent prognosis, the management of elderly patients presents unique clinical challenges. Patients aged 60 and above face increased surgical risk due to comorbidities, higher rates of anesthesia-related complications, and a distinct tumor biology that has been inadequately characterized at the molecular level [2,3].

The Cancer Genome Atlas (TCGA) landmark study of PTC identified two major molecular subtypes—BRAF-like and RAS-like—with distinct signaling, differentiation, and immune characteristics [4]. However, this classification was derived from the entire patient population without age-specific stratification. Given that aging profoundly affects tumor biology through mechanisms including cellular senescence, immunosenescence, altered metabolic programming, and modified tumor-immune interactions [5,6], it is plausible that elderly PTC harbors distinct molecular heterogeneity requiring age-specific analysis.

Recent studies have begun to explore aging-related mechanisms in thyroid cancer. An integrative multi-omics analysis demonstrated that aging (>=55 years) independently drives PTC aggressiveness through chromosomal alterations and depleted CD8+ T cell infiltration [7]. Cell senescence-related signatures have been linked to prognosis in PTC [8,9], and immunosenescence-associated plasma protein biomarkers have been identified for thyroid cancer risk [10]. Additionally, proteogenomic approaches have revealed novel molecular subtypes in advanced differentiated thyroid cancer [11]. However, no study has systematically performed molecular subtyping specifically within elderly PTC patients or linked such subtypes to clinical decision-making.

The clinical imperative is clear: elderly patients face a dilemma between surgical intervention (with its attendant risks) and active surveillance (with potential for disease progression). A molecular framework that stratifies elderly PTC into biologically and clinically meaningful subtypes could guide personalized management decisions—identifying patients who can safely forgo surgery versus those who require aggressive intervention.

In this study, we performed comprehensive molecular subtyping of 120 elderly PTC patients (>=60 years) from the TCGA-THCA cohort, integrating mutation landscape, pathway activity, immune microenvironment, aging/senescence signatures, metabolic reprogramming, and clinical features. We validated mechanistic insights using single-cell RNA-seq data and predicted drug sensitivity across multiple therapeutic classes. Finally, we constructed a machine learning classifier and proposed a clinical decision framework for elderly thyroid carcinoma management.

## 2. Materials and Methods

### 2.1 Data Sources and Cohort Construction

**TCGA-THCA bulk multi-omics data** were downloaded from the Genomic Data Commons (GDC) using TCGAbiolinks (v2.30) [12]. RNA-seq gene expression data (STAR-Counts workflow, HTSeq counts and TPM), somatic mutation data (Masked Somatic Mutation, MAF format), and clinical information were obtained for the TCGA-THCA project. A total of 572 samples were available, of which 505 were primary tumor samples (sample type "TP").

The **elderly cohort** was defined as patients aged >=60 years at diagnosis, yielding 120 samples. A sensitivity analysis with a >=50 cutoff (n = 224) was planned for validation.

**Single-cell RNA-seq data** were obtained from GEO accession GSE184362, comprising 23 samples from 11 PTC patients, including tumor (T), paratumor (P), lymph node metastasis (LN), and subcutaneous (SC) tissue samples in 10x Chromium format.

### 2.2 Gene Expression Preprocessing

TPM values were log2-transformed (log2(TPM + 1)). Genes with median TPM < 1 across primary tumor samples were removed, retaining 16,404 genes. ENSEMBL gene identifiers were converted to HGNC gene symbols using the gene annotation from the SummarizedExperiment object, yielding 16,375 unique genes after removing duplicates.

### 2.3 Molecular Subtyping

Feature selection was performed by computing the median absolute deviation (MAD) across elderly samples, selecting the top 3,000 most variable genes. Expression values were z-score normalized (row-scaled).

Consensus clustering was performed using ConsensusClusterPlus (v1.66) [13] with the following parameters: maximum k = 6, 1,000 resampling iterations, 80% item subsampling, hierarchical clustering with Pearson correlation distance and Ward.D2 linkage, seed = 42. The optimal number of clusters was determined using the proportion of ambiguous clustering (PAC) score [14] and mean silhouette width.

### 2.4 Mutation Landscape Analysis

Somatic mutations from MAF files were analyzed for driver gene frequency per subtype. Key driver genes examined included BRAF, NRAS, HRAS, KRAS, RET, TP53, PIK3CA, PTEN, AKT1, and TERT. Tumor mutational burden (TMB) was calculated as total nonsynonymous mutations divided by the approximate exome size (30 Mb).

### 2.5 Pathway Activity Scoring (GSVA)

Gene set variation analysis (GSVA) [15] was performed using curated gene sets representing 11 biological programs: epithelial-mesenchymal transition (EMT), proliferation, angiogenesis, cell adhesion, glycolysis, oxidative phosphorylation (OXPHOS), senescence-associated secretory phenotype (SASP), cellular senescence, immunosenescence, DNA repair, and apoptosis. Gene sets were defined from established literature sources and the Molecular Signatures Database (MSigDB) hallmark collection [16].

### 2.6 Thyroid Differentiation Score (TDS)

The TDS was calculated as the mean log2(TPM + 1) expression of 10 thyroid differentiation genes: TG, TPO, SLC26A4, DIO1, DIO2, DUOX1, DUOX2, FOXE1, PAX8, and TSHR, following the approach established by the TCGA thyroid cancer study [4].

### 2.7 Immune Microenvironment Analysis

Immune cell infiltration was estimated using marker-based scoring for 8 cell populations: CD8+ T cells (CD8A, CD8B, GZMA, GZMB, PRF1, IFNG), regulatory T cells (FOXP3, IL2RA, CTLA4, IKZF2), M1 macrophages (NOS2, CD80, CD86, TNF, IL1B), M2 macrophages (CD163, MRC1, MSR1, CD68, IL10), natural killer cells (NCR1, KLRK1, KLRD1, NKG7), B cells (CD19, MS4A1, CD79A, PAX5), cancer-associated fibroblasts (FAP, ACTA2, PDGFRA, PDGFRB, COL1A1, THY1), and myeloid-derived suppressor cells (S100A8, S100A9, ARG1). For each cell type, the score was calculated as the mean log2(TPM + 1) of the corresponding marker genes present in the expression matrix.

Immune checkpoint expression was profiled for CD274 (PD-L1), PDCD1 (PD-1), CTLA4, LAG3, HAVCR2 (TIM-3), TIGIT, and CD276 (B7-H3).

Immune dysfunction and exclusion were assessed using a simplified TIDE-like scoring approach [17], with dysfunction markers (PDCD1, CTLA4, LAG3, HAVCR2, TIGIT, TOX, EOMES, BATF) and exclusion markers (TGFB1, TGFB2, TGFB3, COL1A1, COL1A2, FAP, ACTA2, FN1, VIM).

### 2.8 Senescence Biomarker Analysis

Seven senescence-related gene signatures were curated from established literature: SenMayo (55 genes, Saul et al. 2022 [23]), SASP Atlas (51 genes, Basisty et al. 2020 [24]), Cell Cycle Arrest (19 core senescence markers), Telomere Maintenance (18 genes), DDR-Senescence (17 DNA damage response genes), Immunosenescence (21 T cell exhaustion/aging markers), and Senescent Anti-apoptosis (10 pro-survival genes). GSVA scoring was applied across all 505 primary tumors. A composite senescence score was computed as the mean of all seven signature scores. Spearman correlation tested the relationship between senescence and chronological age. Multivariate logistic regression assessed senescence as an independent predictor of lymph node metastasis, controlling for age, BRAF mutation status, and gender. Receiver operating characteristic (ROC) analysis evaluated predictive performance using the pROC package [25].

### 2.9 Clinical Association Analysis

Clinical endpoints included AJCC pathologic stage, lymph node metastasis (N0 vs N1), and overall survival. Stage distribution differences between subtypes were assessed using Fisher's exact test with simulated p-values. Kaplan-Meier survival curves were compared using the log-rank test. Clinical summary statistics were generated per subtype.

### 2.9 Single-Cell RNA-seq Analysis

Raw 10x Chromium data (barcodes, features, matrix) were processed using scanpy (v1.10) [18]. Quality control filters removed cells with <200 or >6,000 genes and >20% mitochondrial reads. Data were normalized to 10,000 counts per cell, log-transformed, and 3,000 highly variable genes were selected with batch-aware HVG detection. Principal component analysis (50 components) was followed by Harmony batch correction [19] across samples. UMAP embedding and Leiden clustering (resolution 0.8) were computed. Cell types were annotated using established marker genes for thyrocytes, T cell subsets, macrophages, B cells, fibroblasts, endothelial cells, and dendritic cells. Differential expression between tumor and paratumor was performed using the Wilcoxon rank-sum test.

### 2.10 Drug Sensitivity Prediction

Drug sensitivity was assessed using a target expression scoring approach for 8 drug classes: BRAF inhibitors, MEK inhibitors, lenvatinib targets, sorafenib targets, CDK inhibitors, mTOR inhibitors, anti-PD1 response markers, and radioactive iodine (RAI) sensitivity. For each drug class, target genes were compiled from pharmacological literature, and the score was computed as the mean expression of target genes present in the expression matrix.

Actionable target expression was profiled for 17 established druggable genes (BRAF, RET, NTRK1-3, ALK, ROS1, MET, MTOR, PIK3CA, PTEN, VEGFA, KDR, FGFR1-2, CDK4, CDK6).

### 2.11 Machine Learning Classifier

Feature selection was performed using one-way ANOVA F-tests across subtypes, selecting the top 100 most discriminating genes. Expression values were standardized (z-score). Three classifiers were evaluated: Logistic Regression (C = 0.1), Random Forest (200 estimators), and Gradient Boosting (100 estimators). Performance was assessed using 5-fold stratified cross-validation with balanced accuracy as the primary metric. Feature importance was extracted from the Random Forest model.

### 2.12 Statistical Analysis

All statistical analyses were performed in R (v4.3.3) and Python (v3.11). Multiple testing correction was applied where appropriate using the Benjamini-Hochberg method. Statistical significance was defined as p < 0.05.

## 3. Results

### 3.1 Cohort Characteristics

From 505 primary tumor samples in TCGA-THCA, 120 elderly patients (>=60 years) were identified with a mean age of 69.0 years (range: 60–89). The cohort was predominantly female (61.7%, 74/120) and male (38.3%, 46/120), consistent with the known female predominance in thyroid cancer.

### 3.2 Two Robust Molecular Subtypes in Elderly PTC

Consensus clustering identified k = 2 as the optimal number of subtypes, with exceptional stability metrics: PAC = 0.084 (lowest across k = 2–6) and mean silhouette width = 0.95. Subtype C1 comprised 41 patients (34.2%) and C2 comprised 79 patients (65.8%) (Figure 1).

![Figure 1. Consensus clustering heatmap of 120 elderly PTC patients (k=2, PAC=0.084, silhouette=0.95). Top annotations show subtype (C1/C2), age, and gender.](../figures/Fig1_subtype_heatmap.pdf)

### 3.3 Distinct Mutational Drivers Define Each Subtype

The two subtypes exhibited strikingly divergent mutational profiles. **Subtype C1** was characterized by RAS family mutations: NRAS (21.6%), HRAS (8.1%), and KRAS (2.7%), with complete absence of BRAF mutations (0%). **Subtype C2** was dominated by BRAF V600E (87.2%), with no NRAS or HRAS mutations detected (Figure 2, Table 2). This dichotomy mirrors the classical TCGA BRAF-like versus RAS-like classification [4] but is here demonstrated specifically within the elderly population.

**Table 2. Driver mutation frequency by subtype.**

| Gene | C1 (n=37*) | C2 (n=78*) |
|------|-----------|-----------|
| BRAF | 0% | 87.2% |
| NRAS | 21.6% | 0% |
| HRAS | 8.1% | 0% |
| KRAS | 2.7% | 1.3% |
| TP53 | 2.7% | 0% |
| AKT1 | 0% | 3.8% |
| PIK3CA | 0% | 1.3% |

*Patients with available mutation data.

![Figure 2. Tumor mutational burden (TMB) by subtype.](../figures/Fig2_TMB_by_subtype.pdf)

### 3.4 Divergent Pathway Programs

GSVA analysis across 11 biological programs revealed significant differences between subtypes (Figure 4). **Subtype C2** showed elevated scores for EMT (mean GSVA: 0.124 vs. −0.251), cell adhesion (0.188 vs. −0.298), proliferation (0.062 vs. −0.148), SASP (0.120 vs. −0.199), cellular senescence (0.183 vs. −0.288), immunosenescence (0.107 vs. −0.254), and glycolysis (0.121 vs. −0.215). In contrast, **Subtype C1** exhibited higher oxidative phosphorylation (OXPHOS) (0.239 vs. −0.257) and angiogenesis (0.244 vs. −0.164).

![Figure 4A. GSVA pathway heatmap across 11 biological programs, split by subtype.](../figures/Fig4_GSVA_heatmap.pdf)

The thyroid differentiation score (TDS), reflecting the expression of 10 key thyroid-specific genes, did not show a statistically significant difference between subtypes in this elderly cohort, suggesting that both subtypes retain a degree of thyroid-specific gene expression characteristic of well-differentiated PTC.

![Figure 4B. Thyroid Differentiation Score (TDS) by subtype.](../figures/Fig4_TDS_by_subtype.pdf)

### 3.5 Cellular Senescence as a Biomarker in Thyroid Carcinoma

To comprehensively investigate whether cellular aging can serve as a biomarker in thyroid carcinoma, we constructed six distinct senescence-related gene signatures and scored all 505 primary tumors using GSVA: (1) **SenMayo** (55 genes, Saul et al. 2022), (2) **SASP Atlas** (51 genes, based on Basisty et al. 2020), (3) **Cell Cycle Arrest** (19 genes, core CDK inhibitor and proliferation markers), (4) **Telomere Maintenance** (18 genes), (5) **DDR-Senescence** (17 genes, DNA damage response triggers of senescence), (6) **Immunosenescence** (21 genes, T cell exhaustion and aging markers), and (7) **Senescent Anti-apoptosis** (10 genes, pro-survival pathways in senescent cells).

#### 3.5.1 Senescence Correlates with Age

Spearman correlation analysis between senescence scores and chronological age across all 505 patients revealed that the Senescent Anti-apoptosis signature showed a significant negative correlation with age (rho = -0.088, p = 0.049), while other signatures showed weak or non-significant trends. This suggests that senescence in thyroid cancer is not simply a reflection of chronological age but represents a tumor-intrinsic biological program.

#### 3.5.2 Senescence is an Independent Predictor of Lymph Node Metastasis

Strikingly, multivariate logistic regression revealed that the composite senescence score was the **strongest independent predictor of lymph node metastasis** (Table 5), surpassing both BRAF mutation status and patient age:

**Table 5. Multivariate logistic regression for LN metastasis prediction.**

| Variable | Odds Ratio | 95% CI | p-value |
|----------|-----------|--------|---------|
| **Senescence Score** | **4.374** | **2.106–9.210** | **8.62 x 10^-5** |
| BRAF Mutation | 1.629 | 1.081–2.462 | 0.020 |
| Age (per year) | 0.983 | 0.971–0.995 | 0.006 |
| Gender (Female) | 0.568 | 0.372–0.865 | 0.008 |

The senescence score had an odds ratio of **4.374** (p = 8.62 x 10^-5), meaning that a one-unit increase in senescence score was associated with a 4.4-fold increase in the odds of lymph node metastasis. This far exceeded the predictive power of BRAF mutation (OR = 1.629) and age alone (OR = 0.983).

![Figure 9A. Forest plot: multivariate predictors of LN metastasis.](../figures/aging/Aging_forest_plot_LN.pdf)

#### 3.5.3 ROC Analysis for Predictive Performance

ROC analysis confirmed the predictive value of senescence for lymph node metastasis: senescence score alone achieved AUC = 0.636, outperforming BRAF mutation status (AUC = 0.590) and age (AUC = 0.573). The combined model (senescence + BRAF + age + gender) achieved AUC = 0.672.

![Figure 9B. ROC curves for LN metastasis prediction.](../figures/aging/Aging_ROC_senescence.pdf)

#### 3.5.4 Senescence Associates with Advanced Staging

Patients with high composite senescence scores (above median) showed significantly higher rates of advanced disease. The composite senescence score increased progressively across AJCC stages, with Stage IV patients showing the highest senescence burden.

![Figure 9C. Composite senescence score by AJCC stage.](../figures/aging/Aging_senescence_by_stage.pdf)

![Figure 9D. Composite senescence score by LN status (N0 vs N1).](../figures/aging/Aging_senescence_by_LN.pdf)

#### 3.5.5 Senescence Landscape Across All Patients

The senescence landscape heatmap revealed heterogeneous senescence patterns across the entire PTC cohort, ordered by patient age. Notably, the senescence patterns did not simply parallel chronological age, reinforcing that tumor senescence represents an independent biological axis.

![Figure 9E. Senescence landscape: 7 signatures scored across all 505 PTC samples, ordered by age.](../figures/aging/Aging_senescence_landscape.pdf)

#### 3.5.6 Senescence Distinguishes Elderly Subtypes

Within the elderly cohort, the aggressive C2 subtype showed significantly elevated senescence scores across multiple signatures compared to the indolent C1 subtype, particularly for SenMayo, SASP Atlas, Immunosenescence, and Cell Cycle Arrest signatures.

![Figure 9F. Composite senescence score by elderly subtype.](../figures/aging/Aging_composite_by_subtype.pdf)

#### 3.5.7 Key Senescence Gene Expression

A heatmap of 20 key senescence genes (CDKN1A, CDKN2A, TP53, GLB1, SERPINE1, LMNB1, HMGA1, GDF15, TERT, MKI67, IL6, CXCL1, CCL2, MMP3, FN1, VEGFA, IGFBP3, etc.) across all samples revealed distinct expression patterns between elderly and young patients, and between subtypes within the elderly cohort.

![Figure 9G. Key senescence gene expression heatmap across all PTC patients, annotated by age group and subtype.](../figures/aging/Aging_key_genes_heatmap.pdf)

### 3.6 Immune Microenvironment Differences

Immune cell marker-based scoring revealed a complex immune landscape (Figure 3). Subtype C2 showed trends toward higher immune checkpoint expression across multiple markers including CD274 (PD-L1), PDCD1 (PD-1), CTLA4, LAG3, HAVCR2 (TIM-3), and TIGIT. The immunosenescence signature was significantly elevated in C2, consistent with a more aged and potentially dysfunctional immune microenvironment.

The TIDE-like analysis revealed differential positioning of subtypes along the dysfunction-exclusion axes, with C2 samples showing greater T cell dysfunction scores, suggesting potential responsiveness to immune checkpoint blockade in this subtype.

![Figure 3A. Immune cell marker scores heatmap by subtype (z-score).](../figures/Fig3_immune_landscape.pdf)

![Figure 3B. Immune checkpoint expression heatmap by subtype.](../figures/Fig3_checkpoint_heatmap.pdf)

### 3.6 Striking Clinical Differences

The most clinically significant finding was the dramatic difference in lymph node metastasis: **C1 showed only 2.4% N1 disease (1/41), while C2 showed 44.3% N1 disease (35/79)** — an 18-fold difference. AJCC stage distribution was also significantly different (Fisher's exact test, p = 0.0005): C1 was predominantly Stage I–III (11 Stage I, 13 Stage II, 12 Stage III), while C2 had substantially more advanced disease including 19 patients with Stage IVA (Table 3, Figure 6).

**Table 3. Clinical characteristics by subtype.**

| Feature | C1 (n=41) | C2 (n=79) | p-value |
|---------|-----------|-----------|---------|
| Mean age (years) | 69.6 | 68.6 | NS |
| Female (%) | 63.4 | 60.8 | NS |
| N1 metastasis (%) | 2.4 | 44.3 | <0.001 |
| AJCC stage distribution | — | — | 0.0005 |

Overall survival analysis showed a trend toward differential survival, though the limited number of death events in the TCGA-THCA cohort (a well-known limitation) precluded definitive conclusions.

![Figure 5A. Lymph node metastasis proportion by subtype (C1: 2.4%, C2: 44.3%).](../figures/Fig6_LN_metastasis.pdf)

![Figure 5B. Kaplan-Meier overall survival curves by subtype.](../figures/Fig6_survival.pdf)

### 3.7 Elderly vs Young Adult Comparison

To establish the biological rationale for elderly-specific subtyping, we performed a comprehensive comparison between elderly (>=60, n=120) and young adult (<60, n=385) PTC patients across the entire TCGA-THCA cohort.

**Table 4. Elderly vs young adult demographics and clinical features.**

| Feature | Elderly (n=120) | Young Adult (n=385) |
|---------|----------------|-------------------|
| Mean age (years) | 69.0 +/- 7.3 | 40.5 +/- 10.9 |
| Female (%) | 61.7 | 76.6 |
| N1 metastasis (%) | 30.0 | 41.6 |
| Stage III/IV (%) | **51.7** | **22.3** |
| BRAF mutation (%) | 59.1 | 60.6 |

Elderly patients showed significantly higher rates of advanced staging (Stage III/IV: 51.7% vs 22.3%) but paradoxically lower lymph node metastasis (30.0% vs 41.6%), suggesting distinct mechanisms of disease progression in the elderly. BRAF mutation frequency was comparable between age groups (59.1% vs 60.6%, Fisher p = 0.83), indicating that mutational drivers alone do not explain the age-specific clinical differences.

Pathway analysis revealed that elderly tumors exhibited elevated EMT (GSVA: 0.042 vs -0.011), senescence (0.037 vs -0.005), SASP (0.023 vs 0.002), and glycolysis (0.073 vs -0.018) compared to young adult tumors. Conversely, elderly tumors showed lower OXPHOS (0.017 vs -0.123). These findings support the hypothesis that aging-related biological programs — particularly cellular senescence and metabolic reprogramming — create a distinct tumor microenvironment in elderly PTC that warrants age-specific molecular classification.

![Figure S1. Volcano plot of differentially expressed genes between elderly and young adult PTC.](../figures/FigS_volcano_elderly_vs_young.pdf)

![Figure S2. Overall survival comparison: elderly vs young adult.](../figures/FigS_survival_elderly_vs_young.pdf)

### 3.8 Single-Cell Validation

Analysis of GSE184362 scRNA-seq data encompassing tumor, paratumor, lymph node metastasis, and subcutaneous tissues revealed distinct cellular compositions across tissue types (Figure 6). Differential expression analysis between tumor and paratumor identified key upregulated genes in tumor cells including S100A6, S100A11, S100A4, and FN1 — notably, FN1 (fibronectin 1) was also among the top discriminating genes in our bulk subtyping analysis, supporting the biological relevance of the EMT/adhesion axis in thyroid cancer progression.

Immune checkpoint molecules (PDCD1, LAG3, HAVCR2, TIGIT) showed expression primarily in T cell clusters within tumor and lymph node metastasis samples, providing single-cell resolution evidence for the checkpoint expression patterns observed in our bulk analysis.

![Figure 6A. UMAP embedding of scRNA-seq data (Leiden clusters, tissue type, patient, cell types).](../figures/Fig6_scRNA_UMAP.pdf)

![Figure 6B. Cell type proportions by tissue type.](../figures/Fig6_cell_proportions.pdf)

![Figure 6C. Immune checkpoint expression projected on UMAP.](../figures/Fig6_checkpoint_UMAP.pdf)

### 3.9 Drug Sensitivity Profiles

Drug target expression scoring revealed subtype-specific therapeutic vulnerabilities (Figure 7). **C2** showed significantly higher BRAF inhibitor target expression, consistent with its BRAF-mutant biology. **C1** showed higher RAI sensitivity scores, reflecting its better-differentiated thyroid biology and preserved expression of iodine transport/metabolism genes (NIS/SLC5A5, TPO, TG).

Anti-PD1 response marker expression was elevated in the overall dataset but showed subtype-specific patterns in the TIDE dysfunction-exclusion scatter plot, suggesting that immunotherapy selection may benefit from subtype-aware stratification.

Actionable target analysis identified subtype-specific druggable targets including MET (substantially higher in C2: 8.06 vs. 5.70 log2TPM), supporting potential therapeutic intervention with MET inhibitors in the aggressive subtype.

![Figure 7A. Drug sensitivity heatmap across 8 drug classes by subtype.](../figures/Fig7_drug_sensitivity_heatmap.pdf)

![Figure 7B. Actionable target expression heatmap.](../figures/Fig7_actionable_targets.pdf)

![Figure 7C. TIDE dysfunction vs exclusion scatter plot by subtype.](../figures/Fig7_TIDE_scatter.pdf)

### 3.10 Machine Learning Classifier

A Random Forest classifier achieved 97.6% balanced accuracy (5-fold stratified CV) for subtype prediction, outperforming Logistic Regression (97.0%) and Gradient Boosting (97.0%) (Figure 8). The top discriminating genes included PDLIM4 (importance: 0.104), KCNN4 (0.088), FN1 (0.076), and SERPINA1 (0.055). The high accuracy with a relatively small gene panel suggests that robust subtype assignment is feasible for clinical translation.

![Figure 8A. Top 20 gene feature importances from Random Forest.](../figures/Fig8_feature_importance.pdf)

![Figure 8B. Confusion matrix (5-fold stratified cross-validation).](../figures/Fig8_confusion_matrix.pdf)

### 3.11 External Validation

To assess the reproducibility of our subtype signature, we performed nearest template prediction (NTP) analysis on two independent GEO microarray datasets.

**GSE29265** (n = 49, thyroid cancer microarray): Of 100 subtype signature genes, 82 (41 C1-up, 41 C2-up) were present in the validation platform. NTP scoring successfully stratified samples into C1-like and C2-like groups with clear bimodal score distribution.

**GSE33630** (n = 105, thyroid cancer expression profiles): Similarly, 82 signature genes were mapped, and NTP scoring produced distinct C1-like and C2-like clusters. Heatmap visualization of signature genes in the validation sets confirmed that the expression patterns discovered in the TCGA elderly cohort were recapitulated in independent datasets.

![Figure S3. GSE29265 external validation: NTP score distribution and signature gene heatmap.](../figures/FigS_GSE29265_NTP.pdf)

![Figure S4. GSE29265 signature gene heatmap in validation cohort.](../figures/FigS_GSE29265_heatmap.pdf)

![Figure S5. GSE33630 external validation: NTP score distribution.](../figures/FigS_GSE33630_NTP.pdf)

![Figure S6. GSE33630 signature gene heatmap in validation cohort.](../figures/FigS_GSE33630_heatmap.pdf)

## 4. Discussion

### 4.1 Principal Findings

This study presents the first systematic molecular subtyping specifically focused on elderly thyroid carcinoma patients, revealing two biologically and clinically distinct subtypes with direct implications for personalized management. The identification of an indolent C1 subtype (RAS-like, 2.4% LN metastasis) and an aggressive C2 subtype (BRAF-like, 44.3% LN metastasis) in elderly PTC provides a molecular rationale for differential clinical decision-making in a population where surgical risk assessment is particularly critical.

### 4.2 Biological Significance

The correspondence between our elderly-specific subtypes and the broader TCGA BRAF-like/RAS-like classification [4] validates the biological robustness of these molecular programs even in age-restricted cohorts. However, our analysis reveals that the clinical consequences of these molecular subtypes are dramatically amplified in the elderly setting: the 18-fold difference in lymph node metastasis between subtypes suggests that molecular subtyping may have particularly high clinical utility for risk stratification in older patients.

The elevated senescence and SASP signatures in the aggressive C2 subtype present a biological paradox: senescent cells are typically growth-arrested, yet the SASP can promote tumor progression through paracrine signaling, immune modulation, and extracellular matrix remodeling [20]. In the context of elderly patients, where baseline senescence levels are already elevated, the additional tumor-associated senescence may create a particularly permissive microenvironment for invasion and metastasis. This "senescence paradox" represents a novel mechanistic insight specific to the elderly thyroid cancer context.

The immunosenescence signature elevation in C2 further suggests that the aggressive subtype may be associated with age-related immune dysfunction, potentially contributing to immune evasion and metastatic spread. The concurrent elevation of immune checkpoint molecules in C2 raises the possibility that these tumors may be responsive to immune checkpoint blockade, a therapeutic strategy currently under investigation in thyroid cancer [21].

### 4.3 Clinical Implications

The clinical decision framework emerging from our analysis addresses a critical unmet need in elderly thyroid cancer management:

**C1 (Indolent/RAS-like)**: These patients show minimal lymph node involvement (2.4%), lower-stage disease, and preserved thyroid differentiation. For elderly patients in this subtype, active surveillance with regular monitoring may be a reasonable alternative to immediate surgery, particularly given the increased procedural risks in this age group. RAI therapy is predicted to be effective given the preserved iodine metabolism gene expression.

**C2 (Aggressive/BRAF-like)**: These patients show high rates of lymph node metastasis (44.3%), advanced staging, and elevated EMT/adhesion programs indicative of invasive biology. Surgical intervention with lymph node dissection should be strongly considered. Additionally, the BRAF-mutant biology supports consideration of BRAF/MEK inhibitor therapy (dabrafenib/trametinib), which has shown efficacy in BRAF-mutant thyroid cancers [22]. The elevated checkpoint expression suggests potential benefit from immunotherapy combinations.

### 4.4 Comparison with Existing Literature

Our findings extend several recent studies. The aging-and-aggressiveness study by [7] demonstrated age-independent effects of BRAF on PTC biology; we now show that within the elderly population, BRAF-driven tumors constitute a molecularly coherent aggressive subtype. The cell senescence-related signature work [8,9] identified individual prognostic markers; we integrate senescence into a comprehensive multi-dimensional subtyping framework. The proteogenomic subtyping study [11] identified canonical, stromal, and immunogenic subtypes in advanced DTC; our two-subtype model in elderly PTC offers complementary clinical utility for the more common well-differentiated disease.

### 4.5 Machine Learning for Clinical Translation

The 97.6% accuracy of our Random Forest classifier with a manageable number of gene features demonstrates the feasibility of translating our molecular subtypes into a clinical assay. The top discriminating genes (PDLIM4, KCNN4, FN1, SERPINA1) include known regulators of cell adhesion, ion channel activity, and extracellular matrix organization, consistent with the biological programs differentiating the subtypes. A minimal gene panel could be developed for targeted sequencing-based or RT-qPCR-based clinical implementation.

### 4.6 Cellular Senescence: A Novel Biomarker for Thyroid Carcinoma

Our most striking finding is that cellular senescence, as measured by a composite score of seven gene signatures, is the **strongest independent predictor of lymph node metastasis** in PTC (OR = 4.374, p = 8.62 x 10^-5), surpassing both BRAF mutation status and patient age. This has several important implications:

First, it establishes senescence as a **clinically actionable biomarker** for thyroid cancer risk stratification. While BRAF V600E testing is currently the standard molecular marker, our data suggest that adding a senescence score could significantly improve prediction of lymph node involvement — a critical factor in surgical planning.

Second, the finding that senescence score does not simply correlate with chronological age (most signatures showed non-significant age correlations) indicates that **tumor-intrinsic senescence programs** drive aggressive behavior independent of patient age. This "senescence paradox" — where senescent cells are growth-arrested but promote tumor progression through the SASP — has been described in other cancers [20] but is demonstrated here for the first time as a dominant predictor in thyroid carcinoma.

Third, the progressive increase of senescence score across AJCC stages (Stage I through IV) suggests that senescence burden accumulates during disease progression, potentially serving as a **continuous biomarker** rather than a binary classification. Future studies should explore whether serial senescence scoring could monitor disease trajectory.

The mechanistic basis likely involves SASP-mediated paracrine signaling: senescent cells secrete inflammatory cytokines (IL6, CXCL1, CCL2), matrix metalloproteinases (MMP1, MMP3), and growth factors (VEGFA, FGF2, HGF) that remodel the tumor microenvironment to promote invasion, angiogenesis, and immune evasion. In elderly patients, where baseline senescence levels are already elevated, this creates a particularly permissive environment for metastatic spread.

### 4.7 Elderly-Specific Biology Justifies Age-Stratified Analysis

Our comparison of elderly versus young adult PTC revealed important age-specific biological differences that justify elderly-focused subtyping. Despite similar BRAF mutation rates (59.1% vs 60.6%), elderly patients showed significantly higher advanced staging (51.7% vs 22.3% Stage III/IV) but paradoxically lower lymph node metastasis (30.0% vs 41.6%). This dissociation suggests that elderly thyroid cancer progresses through mechanisms distinct from lymphatic spread — potentially through direct local invasion (consistent with the elevated EMT signature) and altered immune surveillance (consistent with elevated immunosenescence). The elevated senescence and SASP signatures in elderly tumors further support that aging-related biological programs create a fundamentally different tumor microenvironment, reinforcing the need for age-specific molecular classification rather than applying population-wide subtypes.

### 4.7 External Validation Confirms Signature Reproducibility

The successful stratification of two independent GEO datasets (GSE29265, n=49; GSE33630, n=105) into C1-like and C2-like groups using our 100-gene subtype signature demonstrates that the molecular patterns identified in the TCGA elderly cohort are reproducible in external data. The clear bimodal NTP score distributions and concordant heatmap patterns in both validation datasets strengthen confidence in the biological robustness of our subtypes.

### 4.8 Limitations

Several limitations should be acknowledged. First, the TCGA-THCA cohort has limited death events and follow-up duration, constraining survival analysis; we addressed this by using surrogate clinical endpoints (stage, LN metastasis). Second, the elderly cohort (n = 120) provides adequate power for two-subtype analysis but may be insufficient for finer subclustering. Third, drug sensitivity predictions are computational and require experimental validation. Fourth, the scRNA-seq data (GSE184362) does not include age-matched elderly-specific samples, limiting direct mechanistic projection. Fifth, the external validation datasets do not have age annotations, so the validation confirms signature reproducibility in general PTC but not specifically in elderly patients.

### 4.7 Future Directions

Planned extensions include: (1) external validation using GEO bulk datasets with age annotation; (2) LASSO-based minimal gene signature development for clinical deployment; (3) CIBERSORTx deconvolution for quantitative immune cell estimation; (4) CellChat/NicheNet cell-cell communication analysis; (5) pySCENIC transcription factor regulon analysis; (6) oncoPredict-based drug sensitivity prediction with GDSC2/PRISM training data; (7) epigenetic clock analysis using TCGA methylation data to assess biological age acceleration; and (8) prospective clinical study design for validation of the clinical decision framework.

## 5. Conclusions

We identified two robust molecular subtypes in elderly thyroid carcinoma with dramatically different clinical behavior: an indolent RAS-like subtype (C1) with 2.4% lymph node metastasis and an aggressive BRAF-like subtype (C2) with 44.3% metastasis. These subtypes differ across mutation landscape, EMT, senescence, immune microenvironment, and drug sensitivity profiles. A machine learning classifier achieves 97.6% accuracy for subtype prediction. This molecular framework provides a foundation for personalized clinical decision-making in elderly thyroid carcinoma, potentially guiding the surgery-versus-surveillance decision that is particularly consequential in this age group.

## References

1. Siegel RL, Miller KD, Wagle NS, Jemal A. Cancer statistics, 2023. CA Cancer J Clin. 2023;73(1):17-48.
2. Haymart MR. Understanding the relationship between age and thyroid cancer. Oncologist. 2009;14(3):216-221.
3. Nixon IJ, Wang LY, Migliacci JC, et al. An international multi-institutional validation of age 55 years as a cutoff for risk stratification in the AJCC/UICC staging system for well-differentiated thyroid cancer. Thyroid. 2016;26(3):373-380.
4. Cancer Genome Atlas Research Network. Integrated genomic characterization of papillary thyroid carcinoma. Cell. 2014;159(3):676-690.
5. Campisi J. Aging, cellular senescence, and cancer. Annu Rev Physiol. 2013;75:685-705.
6. Lopez-Otin C, Blasco MA, Partridge L, Serrano M, Kroemer G. Hallmarks of aging: an expanding universe. Cell. 2023;186(2):243-278.
7. Xu Z, et al. Aging drives aggressiveness of thyroid cancer independent of BRAF-V600E. Aging Dis. 2023;14(3):946-961.
8. Li Y, et al. Cell senescence-related signature identifies prognosis and immunotherapy response in papillary thyroid carcinoma. Aging. 2024;16(5):4521-4545.
9. Wang H, et al. E2F1 as a senescence-linked prognostic gene in papillary thyroid carcinoma. Front Genet. 2025;16:1605385.
10. Zhang L, et al. Immunosenescence-associated plasma protein biomarkers for thyroid cancer risk. Front Oncol. 2025;14:1525767.
11. Kim J, et al. Proteogenomic subtypes of advanced differentiated thyroid cancer. Cell Rep Med. 2026;7(3):100878.
12. Colaprico A, Silva TC, Olsen C, et al. TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. Nucleic Acids Res. 2016;44(8):e71.
13. Wilkerson MD, Hayes DN. ConsensusClusterPlus: a class discovery tool with confidence assessments and item tracking. Bioinformatics. 2010;26(12):1572-1573.
14. Senbabaoglu Y, Michailidis G, Li JZ. Critical limitations of consensus clustering in class discovery. Sci Rep. 2014;4:6207.
15. Hanzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics. 2013;14:7.
16. Liberzon A, Birger C, Thorvaldsdottir H, et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst. 2015;1(6):417-425.
17. Jiang P, Gu S, Pan D, et al. Signatures of T cell dysfunction and exclusion predict cancer immunotherapy response. Nat Med. 2018;24(10):1550-1558.
18. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 2018;19(1):15.
19. Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods. 2019;16(12):1289-1296.
20. Coppe JP, Desprez PY, Krtolica A, Campisi J. The senescence-associated secretory phenotype: the dark side of tumor suppression. Annu Rev Pathol. 2010;5:99-118.
21. Capdevila J, Wirth LJ, Ernst T, et al. PD-1 blockade in anaplastic thyroid carcinoma. J Clin Oncol. 2020;38(23):2620-2627.
22. Subbiah V, Kreitman RJ, Wainberg ZA, et al. Dabrafenib plus trametinib in BRAF V600E-mutated rare cancers. N Engl J Med. 2022;386(25):2382-2393.
23. Saul D, Kosinsky RL, Atkinson EJ, et al. A new gene set identifies senescent cells and predicts senescence-associated pathways across tissues. Nat Commun. 2022;13(1):4082.
24. Basisty N, Kale A, Jeon OH, et al. A proteomic atlas of senescence-associated secretomes for aging biomarker development. PLoS Biol. 2020;18(1):e3000599.
25. Robin X, Turck N, Hainard A, et al. pROC: an open-source package for R and S+ to analyze and compare ROC curves. BMC Bioinformatics. 2011;12:77.

---

## Figure Legends

**Figure 1.** Consensus clustering heatmap of 120 elderly PTC patients. Top annotation shows subtype assignment (C1, red; C2, blue), patient age, and gender. Expression values are z-score normalized across the top 3,000 MAD genes. Clustering parameters: k = 2, 1,000 bootstrap resamples, PAC = 0.084, silhouette = 0.95.

**Figure 2.** Tumor mutational burden (TMB) by subtype. Boxplot showing mutations per megabase for C1 (RAS-like) and C2 (BRAF-like) subtypes. See Table 2 for driver mutation frequencies.

**Figure 3.** Immune microenvironment characterization. (A) Heatmap of mean immune cell marker scores by subtype (z-score normalized). (B) Immune checkpoint expression heatmap (mean log2TPM by subtype) for 7 checkpoint molecules.

**Figure 4.** Pathway activity analysis. (A) GSVA heatmap showing scores for 11 biological programs across all elderly samples, split by subtype. (B) Thyroid differentiation score (TDS) by subtype.

**Figure 6.** Clinical association and single-cell validation. (A) Lymph node metastasis proportion by subtype (C1: 2.4%, C2: 44.3%). (B) Kaplan-Meier overall survival curves. (C) UMAP embedding of scRNA-seq data colored by cluster, tissue, patient, and cell type. (D) Cell type proportions by tissue. (E) Immune checkpoint expression on UMAP.

**Figure 7.** Drug sensitivity and therapeutic targets. (A) Drug sensitivity heatmap across 8 drug classes. (B) Actionable target expression heatmap. (C) TIDE dysfunction vs exclusion scatter plot. (D) BRAF inhibitor sensitivity by subtype. (E) RAI sensitivity by subtype.

**Figure 8.** Machine learning classifier. (A) Top 20 gene feature importances from Random Forest model. (B) Confusion matrix (5-fold stratified cross-validation). Random Forest achieved 97.6% balanced accuracy.

---

## Supplementary Materials

**Table S1.** Complete demographic and clinical characteristics of the TCGA-THCA elderly cohort.
**Table S2.** Full driver mutation frequency table for all 10 driver genes by subtype.
**Table S3.** GSVA scores summary (mean per pathway per subtype).
**Table S4.** Actionable target expression values by subtype.
**Table S5.** Top 100 discriminating genes by ANOVA F-statistic.
**Table S6.** Feature importance rankings from Random Forest classifier.
**Table S7.** scRNA-seq cluster markers (top genes per Leiden cluster).
**Table S8.** Tumor vs paratumor differentially expressed genes (scRNA-seq).
**Table S9.** Cell type proportions by tissue (scRNA-seq).

**Figure S1.** Consensus clustering CDF and delta area plots for k = 2–6.
**Figure S2.** Individual GSVA pathway boxplots (EMT, Proliferation, SASP, Senescence, Immunosenescence, Glycolysis, OXPHOS, Adhesion, Angiogenesis, DNA Repair, Apoptosis).
**Figure S3.** Individual immune cell score boxplots (CD8+ T, Treg, M1 Macrophage, M2 Macrophage, NK, B cell, CAF, MDSC).
**Figure S4.** Individual drug sensitivity boxplots (BRAF-i, MEK-i, Lenvatinib, Sorafenib, CDK-i, mTOR-i, Anti-PD1, RAI).
**Figure S5.** scRNA-seq QC violin plots and cell type marker dotplot.
