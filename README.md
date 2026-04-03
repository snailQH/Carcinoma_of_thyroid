# Elderly Thyroid Carcinoma Multi-Omics Molecular Subtyping

Molecular subtyping of elderly thyroid carcinoma using multi-omics data for clinical decision support.

## Background

Thyroid carcinoma is the most common endocrine malignancy. While most papillary thyroid cancers (PTC) have an excellent prognosis, elderly patients face unique challenges: higher surgical risk, increased comorbidities, and a distinct tumor biology that is poorly characterized at the molecular level. This project identifies molecular subtypes within elderly thyroid cancer patients to guide personalized clinical decision-making.

## Objectives

1. **Molecular Subtyping** -- Identify reproducible subtypes in elderly (>=60) PTC using unsupervised clustering on TCGA-THCA RNA-seq data
2. **Subtype Characterization** -- Profile subtypes across mutation landscape, immune microenvironment, aging/senescence signatures, metabolic reprogramming, and clinical features
3. **Single-Cell Validation** -- Validate subtype-driving mechanisms using published scRNA-seq data (GEO)
4. **Drug Sensitivity Prediction** -- Predict drug responses per subtype using GDSC/PRISM/CMap and identify actionable therapeutic targets
5. **Clinical Decision Model** -- Build a gene signature-based classifier and clinical decision framework (surgery vs. surveillance vs. targeted therapy)

## Data Sources

| Source | Data type | Samples |
|--------|-----------|---------|
| TCGA-THCA (GDC) | RNA-seq, mutation, CNV, methylation, clinical | ~500 PTC |
| GEO (GSE184362, etc.) | scRNA-seq (tumor + normal + metastasis) | Multiple |
| GDSC2 / PRISM | Drug sensitivity (cell lines) | Training data |
| CMap / L1000 | Connectivity mapping | Compound signatures |

## Analysis Progress (as of 2026-04-03)

**Overall: 63 / 120 analyses completed (53%)**

| Phase | Status | Key results |
|-------|--------|-------------|
| Phase 1: Data Acquisition | **Done** | TCGA-THCA: 572 samples (RNA-seq + mutation + clinical); GEO GSE184362: 23 scRNA-seq samples |
| Phase 2: Molecular Subtyping | **Done** | **2 subtypes** identified (k=2, PAC=0.084, silhouette=0.95); C1 n=41, C2 n=79 |
| Phase 3: Characterization | **Partially done** | GSVA (11 pathways), immune scoring (8 cell types), TDS, checkpoints, TMB, clinical association done; CIBERSORTx, CNA, TLS pending |
| Phase 4: Single-Cell | **Partially done** | QC, clustering, UMAP, DEGs, checkpoint UMAP done; CellChat, SCENIC, trajectory pending |
| Phase 5: Drug Sensitivity | **Partially done** | Target expression scoring (8 drugs), TIDE, actionable targets done; oncoPredict, CMap pending |
| Phase 6: Clinical Model | **Partially done** | RF classifier 97.6% accuracy, external validation done (GSE29265, GSE33630); LASSO signature pending |
| Phase 6B: Elderly vs Young | **Done** | DEGs, GSVA, immune, mutations, clinical, survival comparison (8 sub-analyses) |
| Phase 7: Outputs | **In progress** | 55+ figures, 15+ tables, 17-slide PPTX with speaker notes, manuscript draft v1 with embedded figures |

### Key Findings

| Subtype | N | LN Metastasis | Stage | Clinical Implication |
|---------|---|---------------|-------|---------------------|
| **C1 (Indolent)** | 41 | 2.4% | Mostly I-III | Active surveillance candidate |
| **C2 (Aggressive)** | 79 | 44.3% | More III-IV | Surgical intervention candidate |

- Stage distribution significantly different between subtypes (**Fisher p = 0.0005**)
- ML classifier achieves **97.6% balanced accuracy** (Random Forest, 5-fold CV)
- **Elderly vs Young**: elderly show higher Stage III/IV (51.7% vs 22.3%), higher EMT/senescence, but lower LN metastasis (30% vs 41.6%)
- **External validation**: subtype signature successfully stratifies GSE29265 (n=49) and GSE33630 (n=105)

For detailed checkbox tracking of all 120 sub-analyses, see [`docs/plans/2026-04-02-elderly-thyroid-carcinoma-multiomics-design.md`](docs/plans/2026-04-02-elderly-thyroid-carcinoma-multiomics-design.md).

---

## Project Structure

```
Carcinoma_of_thyroid/
├── docs/
│   ├── Data_Structure.md                    # Server and data inventory
│   ├── carcinoma_of_thyroid_study_ideas.md  # Initial study design notes
│   └── plans/                               # Detailed research plans with progress tracker
├── scripts/
│   ├── 01_data_download/                    # TCGA + GEO download scripts
│   ├── 02_preprocessing/                    # QC, normalization, cohort construction, subtyping
│   ├── 04_characterization/                 # Mutation, immune, GSVA, senescence, clinical
│   ├── 05_scrna_analysis/                   # Single-cell pipeline (scanpy)
│   ├── 06_drug_prediction/                  # Drug sensitivity, actionable targets
│   ├── 07_clinical_model/                   # ML classifier, PPTX generation
│   └── run_all.sh                           # Master pipeline script
├── research/
│   ├── figures/                             # All generated figures (45 PDFs)
│   ├── tables/                              # All summary tables (11 CSVs)
│   ├── pptx/                                # Summary presentations
│   └── manuscripts/                         # Manuscript drafts (future)
├── envs/                                    # Conda environment YAML files
├── CLAUDE.md                                # AI assistant context
└── README.md
```

**Server data directory** (`bioinfo@192.168.100.127:/data2/projects/Carcinoma_of_thyroid/`):
- `data/tcga_thca/` — TCGA-THCA downloads (RNA-seq, mutation, clinical) ~2.9 GB
- `data/geo_scrna/GSE184362/` — scRNA-seq raw + processed h5ad ~7.3 GB
- `results/` — mirrors `research/` folder in this repo

## Quick Start

```bash
# SSH into the analysis server
ssh bioinfo@192.168.100.127

# Run the full pipeline
cd /data2/projects/Carcinoma_of_thyroid
bash scripts/run_all.sh

# Or run individual steps:
conda run -n thyroid_r Rscript scripts/01_data_download/01_download_tcga.R
conda run -n thyroid_r Rscript scripts/02_preprocessing/02_preprocess_and_subtype.R
conda run -n thyroid_r Rscript scripts/04_characterization/03_characterization.R
conda run -n thyroid_r Rscript scripts/06_drug_prediction/04_drug_prediction.R
conda run -n thyroid_py python scripts/05_scrna_analysis/05_scrna_analysis.py
conda run -n thyroid_py python scripts/07_clinical_model/06_clinical_model_and_summary.py
```

## Methods Summary

- **Subtyping**: Consensus Clustering (k=2-6, 1000 bootstrap replicates) on top 3000 MAD genes
- **Pathway scoring**: GSVA/ssGSEA with Hallmark + curated gene sets
- **Immune deconvolution**: Marker-based scoring (CD8+ T, Treg, M2 macrophage, CAF, etc.)
- **Senescence**: SenMayo gene set, SASP, telomere maintenance, immunosenescence markers
- **Drug prediction**: Target expression scoring + oncoPredict (GDSC2/PRISM ridge regression)
- **Classifier**: Random Forest / Gradient Boosting with 5-fold stratified CV
- **Single-cell**: scanpy pipeline with CellTypist annotation

## Key Outputs

| Output | Path | Description |
|--------|------|-------------|
| Subtype heatmap | `results/figures/Fig1_subtype_heatmap.pdf` | Consensus clustering visualization |
| OncoPrint | `results/figures/Fig2_oncoplot.pdf` | Mutation landscape per subtype |
| Immune landscape | `results/figures/Fig3_immune_landscape.pdf` | Immune cell composition |
| GSVA heatmap | `results/figures/Fig4_GSVA_heatmap.pdf` | Pathway activity scores |
| scRNA UMAP | `results/figures/Fig6_scRNA_UMAP.pdf` | Single-cell clusters |
| Drug sensitivity | `results/figures/Fig7_drug_sensitivity_heatmap.pdf` | Predicted drug response |
| Classifier | `results/figures/Fig8_feature_importance.pdf` | ML model performance |
| Summary PPTX | `results/reports/Elderly_Thyroid_Carcinoma_Summary.pptx` | Full presentation |
| Demographics | `results/tables/Table1_demographics.csv` | Cohort characteristics |

## Dependencies

**R environment** (`thyroid_r`): R 4.3, TCGAbiolinks, ConsensusClusterPlus, GSVA, maftools, ComplexHeatmap, survival, survminer, clusterProfiler, ggplot2, ggpubr

**Python environment** (`thyroid_py`): Python 3.11, scanpy, scikit-learn, xgboost, python-pptx, celltypist, decoupler, matplotlib, seaborn

## License

For academic research use only.
