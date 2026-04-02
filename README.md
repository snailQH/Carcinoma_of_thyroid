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

## Project Structure

```
Carcinoma_of_thyroid/
├── docs/
│   ├── Data_Structure.md              # Server and data inventory
│   ├── carcinoma_of_thyroid_study_ideas.md  # Initial study design notes
│   └── plans/                         # Detailed research plans
├── CLAUDE.md                          # AI assistant context
└── README.md
```

**Server analysis directory** (`bioinfo@192.168.100.127:/data2/projects/Carcinoma_of_thyroid/`):

```
├── data/
│   ├── tcga_thca/{rnaseq,mutation,cnv,methylation,clinical}
│   ├── geo_scrna/{GSE184362,GSE191288}
│   ├── geo_bulk_validation/
│   └── drug_databases/
├── scripts/
│   ├── 01_data_download/        # TCGA + GEO download scripts
│   ├── 02_preprocessing/        # QC, normalization, cohort construction, subtyping
│   ├── 04_characterization/     # Mutation, immune, GSVA, senescence, clinical
│   ├── 05_scrna_analysis/       # Single-cell pipeline (scanpy)
│   ├── 06_drug_prediction/      # Drug sensitivity, actionable targets
│   ├── 07_clinical_model/       # ML classifier, PPTX generation
│   └── run_all.sh               # Master pipeline script
├── results/
│   ├── figures/                 # Publication-ready plots (PDF)
│   ├── tables/                  # Summary tables (CSV)
│   └── reports/                 # PPTX summary, pipeline log
└── envs/                        # Conda environment YAML files
```

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
