#!/bin/bash
# =============================================================================
# run_all.sh — Master pipeline for Elderly Thyroid Carcinoma analysis
# =============================================================================
set -e

export PATH=/opt/conda/bin:$PATH
BASE="/data2/projects/Carcinoma_of_thyroid"
SCRIPTS="$BASE/scripts"
LOG="$BASE/results/reports/pipeline.log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG"; }

log "=== Starting Elderly Thyroid Carcinoma Analysis Pipeline ==="

# Phase 1: Data download
log "Phase 1: Downloading TCGA-THCA data..."
conda run -n thyroid_r Rscript "$SCRIPTS/01_data_download/01_download_tcga.R" 2>&1 | tee -a "$LOG"

# Phase 2: Preprocessing + Subtyping
log "Phase 2: Preprocessing and molecular subtyping..."
conda run -n thyroid_r Rscript "$SCRIPTS/02_preprocessing/02_preprocess_and_subtype.R" 2>&1 | tee -a "$LOG"

# Phase 3: Characterization
log "Phase 3: Subtype characterization..."
conda run -n thyroid_r Rscript "$SCRIPTS/04_characterization/03_characterization.R" 2>&1 | tee -a "$LOG"

# Phase 4: Drug prediction
log "Phase 4: Drug sensitivity prediction..."
conda run -n thyroid_r Rscript "$SCRIPTS/06_drug_prediction/04_drug_prediction.R" 2>&1 | tee -a "$LOG"

# Phase 5: Single-cell analysis (if data available)
log "Phase 5: Single-cell RNA-seq analysis..."
conda run -n thyroid_py python "$SCRIPTS/05_scrna_analysis/05_scrna_analysis.py" 2>&1 | tee -a "$LOG" || log "WARNING: scRNA-seq analysis skipped (data may not be downloaded yet)"

# Phase 6: Clinical model + PPTX
log "Phase 6: Clinical model and summary generation..."
conda run -n thyroid_py python "$SCRIPTS/07_clinical_model/06_clinical_model_and_summary.py" 2>&1 | tee -a "$LOG"

log "=== Pipeline complete ==="
log "Results: $BASE/results/"
log "Figures: $BASE/results/figures/"
log "Tables:  $BASE/results/tables/"
log "PPTX:    $BASE/results/reports/Elderly_Thyroid_Carcinoma_Summary.pptx"
