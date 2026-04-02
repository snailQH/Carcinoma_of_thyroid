# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a bioinformatics research project on **elderly thyroid carcinoma molecular subtyping**. The goal is to use multi-omics data (TCGA-THCA) and published single-cell RNA-seq data to identify molecular subtypes in elderly thyroid cancer patients (≥60 or ≥50), characterize their immune microenvironment, tumor biology, and mutation landscape, and build a clinical decision-oriented classification framework.

## Research Design

The project follows a 5-aim structure:
1. **Aim 1 — Molecular subtyping**: Consensus clustering / NMF / iCluster+ / MOVICS on elderly TCGA-THCA RNA-seq → k=2–4 subtypes
2. **Aim 2 — Subtype characterization**: Mutation landscape, CNA, pathway activity (GSVA), immune deconvolution, aging/senescence signatures, metabolic reprogramming, thyroid differentiation score (TDS), clinical association
3. **Aim 3 — Single-cell mechanism**: scRNA-seq for cell-type composition, tumor cell states, CellChat/NicheNet communication, pySCENIC regulon analysis, pseudotime trajectory
4. **Aim 4 — Drug sensitivity & targets**: oncoPredict/pRRophetic (GDSC/PRISM), CMap/L1000 connectivity mapping, TIDE/IPS immunotherapy prediction, DGIdb/OncoKB actionable targets
5. **Aim 5 — Clinical decision model**: LASSO gene signature, risk score, ML classifier, external validation, decision framework (subtype → management)

## Key Data Sources

- **TCGA-THCA**: ~500 PTC samples, RNA-seq/mutation/CNV/methylation/clinical via GDC
- **GEO scRNA-seq**: e.g. GSE184362 (tumor + LN metastasis + normal)
- **Validation**: GEO bulk datasets, Human Protein Atlas (protein-level)
- **Drug response**: GDSC2 (Sanger), PRISM (Broad), CMap/L1000 (connectivity mapping)
- **Supplementary**: cBioPortal, UCSC Xena, Human Protein Atlas, TISCH, STRING, DGIdb, OncoKB

## Detailed Plan

See `docs/plans/2026-04-02-elderly-thyroid-carcinoma-multiomics-design.md` for the full research plan including methods, timeline, figure layout, tools, and novel ideas.

## Remote Server

- **Host**: `bioinfo@192.168.100.127`
- **Data path**: `/data2/projects/Carcinoma_of_thyroid/`
- **Access**: passwordless SSH, sudo-enabled
- All data downloading and computational analysis should run on this server

## Important Caveats (from study design)

- TCGA-THCA has **weak survival data** (low death events, short follow-up) — use surrogate clinical endpoints (stage, metastasis, recurrence) instead of OS/PFS as primary outcomes
- Elderly cohort sample size may be limited — plan sensitivity analyses across age cutoffs (50+ vs 60+)
- Subtyping stability must be validated with bootstrap + external cohort
- Clinical recommendations must be framed as "potential stratification for clinical decision-making", not direct treatment advice
