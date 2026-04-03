#!/usr/bin/env python3
"""
05_scrna_analysis.py — Aim 3: Single-cell RNA-seq analysis (GSE184362)
Multi-sample PTC: Tumor, Paratumor, LN metastasis
"""

import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

base_dir = "/data2/projects/Carcinoma_of_thyroid"
data_dir = os.path.join(base_dir, "data/geo_scrna/GSE184362")
results_dir = os.path.join(base_dir, "results")
fig_dir = os.path.join(results_dir, "figures")

sc.settings.figdir = fig_dir
sc.settings.set_figure_params(dpi=150, facecolor='white', frameon=False)

# ============================================================================
# 1. Load all samples from 10x format
# ============================================================================
print("=== Loading multi-sample 10x data ===")

# Find all unique samples (each has barcodes/features/matrix)
import glob
barcode_files = sorted(glob.glob(os.path.join(data_dir, "*_barcodes.tsv.gz")))

adatas = []
for bf in barcode_files:
    basename = os.path.basename(bf).replace("_barcodes.tsv.gz", "")
    # Extract sample name (after GSMxxxxxxx_)
    parts = basename.split("_", 1)
    gsm_id = parts[0]
    sample_name = parts[1] if len(parts) > 1 else basename

    matrix_file = bf.replace("_barcodes.tsv.gz", "_matrix.mtx.gz")
    features_file = bf.replace("_barcodes.tsv.gz", "_features.tsv.gz")

    if not (os.path.exists(matrix_file) and os.path.exists(features_file)):
        print(f"  Skipping {sample_name}: missing files")
        continue

    # Create temp dir structure for read_10x_mtx
    tmp_dir = os.path.join(data_dir, f"tmp_{sample_name}")
    os.makedirs(tmp_dir, exist_ok=True)
    os.system(f"ln -sf '{matrix_file}' '{tmp_dir}/matrix.mtx.gz'")
    os.system(f"ln -sf '{features_file}' '{tmp_dir}/features.tsv.gz'")
    os.system(f"ln -sf '{bf}' '{tmp_dir}/barcodes.tsv.gz'")

    try:
        adata = sc.read_10x_mtx(tmp_dir, var_names='gene_symbols', cache=False)
        adata.obs['sample'] = sample_name
        adata.obs['gsm_id'] = gsm_id

        # Parse tissue type
        if '_T' in sample_name and 'LN' not in sample_name:
            tissue = 'Tumor'
        elif '_P' in sample_name:
            tissue = 'Paratumor'
        elif 'LN' in sample_name:
            tissue = 'LN_metastasis'
        elif '_SC' in sample_name:
            tissue = 'Subcutaneous'
        else:
            tissue = 'Unknown'
        adata.obs['tissue'] = tissue

        # Patient ID
        patient = sample_name.split('_')[0]
        adata.obs['patient'] = patient

        adatas.append(adata)
        print(f"  Loaded {sample_name}: {adata.n_obs} cells, {adata.n_vars} genes [{tissue}]")
    except Exception as e:
        print(f"  Error loading {sample_name}: {e}")

    # Clean up symlinks
    os.system(f"rm -rf '{tmp_dir}'")

# Concatenate
print(f"\nConcatenating {len(adatas)} samples...")
adata = sc.concat(adatas, join='outer', fill_value=0)
adata.obs_names_make_unique()
adata.var_names_make_unique()
print(f"Total: {adata.n_obs} cells, {adata.n_vars} genes")
print(f"Tissue distribution:\n{adata.obs['tissue'].value_counts()}")

# ============================================================================
# 2. QC and preprocessing
# ============================================================================
print("\n=== QC and preprocessing ===")

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# QC plots
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
sc.pl.violin(adata, ['n_genes_by_counts'], groupby='tissue', ax=axes[0], show=False, rotation=45)
sc.pl.violin(adata, ['total_counts'], groupby='tissue', ax=axes[1], show=False, rotation=45)
sc.pl.violin(adata, ['pct_counts_mt'], groupby='tissue', ax=axes[2], show=False, rotation=45)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "scrna_QC_violin.pdf"))
plt.close()

# Filter
n_before = adata.n_obs
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=6000)
adata = adata[adata.obs.pct_counts_mt < 20, :].copy()
sc.pp.filter_genes(adata, min_cells=10)
print(f"Cells: {n_before} -> {adata.n_obs} after QC")
print(f"Genes: {adata.n_vars}")

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()

# HVG, scale, PCA
sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key='sample')
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata, n_comps=50)

# Batch correction with Harmony
print("Running Harmony batch correction...")
try:
    sc.external.pp.harmony_integrate(adata, key='sample')
    sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_pcs=30)
except:
    print("Harmony not available, using PCA directly")
    sc.pp.neighbors(adata, n_pcs=30)

# UMAP and clustering
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.8)
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_0.5')

print(f"Clusters (res=0.8): {adata.obs['leiden'].nunique()}")

# ============================================================================
# 3. Cell type annotation
# ============================================================================
print("\n=== Cell type annotation ===")

marker_genes = {
    'Thyrocyte': ['TG', 'TPO', 'PAX8', 'TSHR', 'SLC5A5'],
    'T_CD8': ['CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1', 'NKG7'],
    'T_CD4': ['CD4', 'IL7R', 'CCR7', 'LEF1', 'TCF7'],
    'Treg': ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2'],
    'T_exhausted': ['PDCD1', 'LAG3', 'HAVCR2', 'TOX', 'TIGIT', 'ENTPD1'],
    'B_cell': ['CD19', 'MS4A1', 'CD79A', 'PAX5'],
    'Plasma': ['JCHAIN', 'MZB1', 'SDC1', 'XBP1', 'IGHG1'],
    'Macrophage_M1': ['CD68', 'NOS2', 'CD80', 'CD86', 'TNF'],
    'Macrophage_M2': ['CD68', 'CD163', 'MRC1', 'MSR1', 'IL10'],
    'DC': ['CLEC9A', 'CD1C', 'ITGAX', 'BATF3', 'FLT3'],
    'NK': ['KLRD1', 'NKG7', 'NCAM1', 'GNLY', 'KLRK1'],
    'Endothelial': ['PECAM1', 'CDH5', 'VWF', 'ERG', 'FLT1'],
    'Fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM', 'FAP', 'ACTA2'],
    'Mast': ['KIT', 'TPSAB1', 'TPSB2', 'CPA3']
}

# Score each cell type
for ct, genes in marker_genes.items():
    present = [g for g in genes if g in adata.raw.var_names]
    if present:
        sc.tl.score_genes(adata, gene_list=present, score_name=f'{ct}_score', use_raw=True)

# Try celltypist
try:
    import celltypist
    from celltypist import models
    models.download_models(force_update=False, model='Immune_All_Low.pkl')
    predictions = celltypist.annotate(adata.raw, model='Immune_All_Low.pkl', majority_voting=True)
    adata.obs['celltypist'] = predictions.predicted_labels['majority_voting']
    print("CellTypist annotation complete.")
    cell_type_col = 'celltypist'
except Exception as e:
    print(f"CellTypist unavailable: {e}")
    cell_type_col = 'leiden'

# ============================================================================
# 4. Visualizations
# ============================================================================
print("\n=== Generating plots ===")

# Main UMAP panels
fig, axes = plt.subplots(2, 2, figsize=(16, 14))
sc.pl.umap(adata, color='leiden', ax=axes[0, 0], show=False, title='Leiden Clusters', legend_loc='on data', legend_fontsize=8)
sc.pl.umap(adata, color='tissue', ax=axes[0, 1], show=False, title='Tissue Type')
sc.pl.umap(adata, color='patient', ax=axes[1, 0], show=False, title='Patient')
if cell_type_col in adata.obs.columns:
    sc.pl.umap(adata, color=cell_type_col, ax=axes[1, 1], show=False, title='Cell Types', legend_fontsize=6)
else:
    sc.pl.umap(adata, color='n_genes_by_counts', ax=axes[1, 1], show=False, title='Genes/Cell')
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "Fig6_scRNA_UMAP.pdf"), dpi=150, bbox_inches='tight')
plt.close()

# Dotplot of markers
all_markers = []
for ct, genes in marker_genes.items():
    present = [g for g in genes if g in adata.raw.var_names]
    all_markers.extend(present[:3])  # Top 3 per type

if all_markers:
    sc.pl.dotplot(adata, var_names=list(dict.fromkeys(all_markers)),  # deduplicate preserving order
                  groupby='leiden', standard_scale='var',
                  save='_cell_markers.pdf')

# ============================================================================
# 5. Cell type proportions by tissue
# ============================================================================
print("\n=== Cell type proportions ===")

proportions = pd.crosstab(adata.obs['tissue'], adata.obs[cell_type_col], normalize='index')
proportions.to_csv(os.path.join(results_dir, "tables", "cell_type_proportions.csv"))

fig, ax = plt.subplots(figsize=(12, 6))
proportions.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')
ax.set_ylabel('Proportion')
ax.set_title('Cell Type Proportions by Tissue')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=7)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "Fig6_cell_proportions.pdf"), dpi=150, bbox_inches='tight')
plt.close()

# ============================================================================
# 6. Tumor vs Normal DEGs
# ============================================================================
print("\n=== Differential expression: Tumor vs Paratumor ===")

adata_tn = adata.raw.to_adata()[adata.obs['tissue'].isin(['Tumor', 'Paratumor'])].copy()
if adata_tn.n_obs > 100:
    sc.tl.rank_genes_groups(adata_tn, groupby='tissue', reference='Paratumor', method='wilcoxon', n_genes=200)
    sc.pl.rank_genes_groups(adata_tn, n_genes=20, save='_Tumor_vs_Paratumor.pdf')

    deg_df = sc.get.rank_genes_groups_df(adata_tn, group='Tumor')
    deg_df.to_csv(os.path.join(results_dir, "tables", "scrna_DEG_Tumor_vs_Paratumor.csv"), index=False)
    print(f"Top DEGs saved: {len(deg_df)} genes")

# ============================================================================
# 7. Immune checkpoint expression in T cells
# ============================================================================
print("\n=== Immune checkpoint analysis ===")

checkpoints = ['PDCD1', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT', 'TOX', 'CD274', 'CD276']
present_cp = [g for g in checkpoints if g in adata.raw.var_names]

if present_cp:
    fig, axes = plt.subplots(2, len(present_cp) // 2 + 1, figsize=(4 * (len(present_cp) // 2 + 1), 8))
    axes = axes.flatten()
    for i, gene in enumerate(present_cp):
        if i < len(axes):
            sc.pl.umap(adata, color=gene, ax=axes[i], show=False, title=gene,
                       use_raw=True, vmin=0, cmap='Reds')
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "Fig6_checkpoint_UMAP.pdf"), dpi=150, bbox_inches='tight')
    plt.close()

# ============================================================================
# 8. DEGs per cluster
# ============================================================================
print("\n=== Cluster marker genes ===")

sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', n_genes=100)
sc.pl.rank_genes_groups(adata, n_genes=15, save='_cluster_markers.pdf')

deg_results = pd.DataFrame()
for group in adata.obs['leiden'].cat.categories:
    try:
        degs = sc.get.rank_genes_groups_df(adata, group=group)
        degs['cluster'] = group
        deg_results = pd.concat([deg_results, degs])
    except:
        pass
deg_results.to_csv(os.path.join(results_dir, "tables", "scrna_cluster_markers.csv"), index=False)

# ============================================================================
# 9. Save processed data
# ============================================================================
print("\n=== Saving processed data ===")
adata.write(os.path.join(data_dir, "processed_adata.h5ad"))
print(f"Saved: {adata.n_obs} cells, {adata.n_vars} genes")

# Summary stats
summary = {
    'Total cells': adata.n_obs,
    'Total genes': adata.raw.n_vars,
    'HVGs used': adata.n_vars,
    'Clusters': adata.obs['leiden'].nunique(),
    'Patients': adata.obs['patient'].nunique(),
    'Tissues': dict(adata.obs['tissue'].value_counts()),
}
print("\n=== Summary ===")
for k, v in summary.items():
    print(f"  {k}: {v}")

print("\n=== scRNA-seq analysis complete ===")
