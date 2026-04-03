#!/usr/bin/env Rscript
# =============================================================================
# 02_preprocess_and_subtype.R — Preprocess TCGA-THCA + elderly cohort + subtyping
# =============================================================================

library(SummarizedExperiment)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(survival)
library(survminer)

base_dir <- "/data2/projects/Carcinoma_of_thyroid"
data_dir <- file.path(base_dir, "data/tcga_thca")
results_dir <- file.path(base_dir, "results")

# ============================================================================
# PART 1: Load and preprocess
# ============================================================================
cat("=== Loading data ===\n")
se <- readRDS(file.path(data_dir, "rnaseq", "TCGA_THCA_rnaseq_SE.rds"))
tpm <- readRDS(file.path(data_dir, "rnaseq", "TCGA_THCA_tpm.rds"))
clinical <- readRDS(file.path(data_dir, "clinical", "TCGA_THCA_clinical.rds"))
sample_info <- readRDS(file.path(data_dir, "clinical", "TCGA_THCA_sample_info.rds"))

# Keep only primary tumor samples (sample type 01)
tumor_idx <- which(sample_info$shortLetterCode == "TP")
tpm_tumor <- tpm[, tumor_idx]
sample_info_tumor <- sample_info[tumor_idx, ]
cat("Primary tumor samples:", ncol(tpm_tumor), "\n")

# Log2 transform
log2tpm <- log2(tpm_tumor + 1)

# Filter low expression genes (median TPM >= 1)
gene_median <- apply(tpm_tumor, 1, median)
keep_genes <- gene_median >= 1
log2tpm_filtered <- log2tpm[keep_genes, ]
cat("Genes after filtering:", nrow(log2tpm_filtered), "\n")

# ============================================================================
# PART 2: Construct elderly cohort
# ============================================================================
cat("\n=== Constructing elderly cohort ===\n")

# Extract age from sample_info (age_at_diagnosis is in days)
sample_info_tumor$age_years <- as.numeric(sample_info_tumor$age_at_index)

# Elderly definition: >=60
elderly_idx <- which(sample_info_tumor$age_years >= 60)
elderly_samples <- sample_info_tumor[elderly_idx, ]
elderly_expr <- log2tpm_filtered[, elderly_idx]
cat("Elderly (>=60) samples:", ncol(elderly_expr), "\n")

# Also prepare >=50 for sensitivity analysis
elderly50_idx <- which(sample_info_tumor$age_years >= 50)
cat("Elderly (>=50) samples:", length(elderly50_idx), "\n")

# ============================================================================
# PART 3: Demographics table (Table 1)
# ============================================================================
cat("\n=== Generating demographics table ===\n")

# Merge clinical info
sample_info_tumor$patient_id <- substr(sample_info_tumor$patient, 1, 12)
clinical$patient_id <- clinical$submitter_id

merged <- merge(
  as.data.frame(sample_info_tumor),
  clinical,
  by.x = "patient_id",
  by.y = "patient_id",
  all.x = TRUE
)

merged$age_group <- ifelse(merged$age_years >= 60, "Elderly (>=60)", "Young (<60)")

# Table 1: demographics
demographics <- merged %>%
  group_by(age_group) %>%
  summarise(
    n = n(),
    age_mean = round(mean(age_years, na.rm = TRUE), 1),
    age_sd = round(sd(age_years, na.rm = TRUE), 1),
    female_n = sum(gender.x == "female", na.rm = TRUE),
    male_n = sum(gender.x == "male", na.rm = TRUE),
    stage_I = sum(grepl("stage i$|stage i ", ajcc_pathologic_stage, ignore.case = TRUE), na.rm = TRUE),
    stage_II = sum(grepl("stage ii", ajcc_pathologic_stage, ignore.case = TRUE), na.rm = TRUE),
    stage_III = sum(grepl("stage iii", ajcc_pathologic_stage, ignore.case = TRUE), na.rm = TRUE),
    stage_IV = sum(grepl("stage iv", ajcc_pathologic_stage, ignore.case = TRUE), na.rm = TRUE),
    .groups = "drop"
  )

write.csv(demographics, file.path(results_dir, "tables", "Table1_demographics.csv"), row.names = FALSE)
cat("Demographics table saved.\n")
print(demographics)

# Save processed data
saveRDS(elderly_expr, file.path(data_dir, "rnaseq", "elderly_log2tpm_filtered.rds"))
saveRDS(elderly_samples, file.path(data_dir, "clinical", "elderly_sample_info.rds"))
saveRDS(log2tpm_filtered, file.path(data_dir, "rnaseq", "all_tumor_log2tpm_filtered.rds"))
saveRDS(sample_info_tumor, file.path(data_dir, "clinical", "all_tumor_sample_info.rds"))
saveRDS(merged, file.path(data_dir, "clinical", "merged_clinical.rds"))

# ============================================================================
# PART 4: Molecular subtyping — Consensus Clustering
# ============================================================================
cat("\n=== Running Consensus Clustering ===\n")

# Select top variable genes by MAD
mad_values <- apply(elderly_expr, 1, mad)
top_genes <- names(sort(mad_values, decreasing = TRUE))[1:min(3000, sum(mad_values > 0))]
expr_for_clustering <- elderly_expr[top_genes, ]
cat("Genes for clustering:", length(top_genes), "\n")

# Scale
expr_scaled <- t(scale(t(expr_for_clustering)))

# Consensus Clustering
cc_dir <- file.path(results_dir, "figures", "consensus_clustering")
dir.create(cc_dir, recursive = TRUE, showWarnings = FALSE)

cc_results <- ConsensusClusterPlus(
  d = as.matrix(expr_scaled),
  maxK = 6,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "hc",
  distance = "pearson",
  innerLinkage = "ward.D2",
  finalLinkage = "ward.D2",
  seed = 42,
  plot = "pdf",
  title = cc_dir
)

# Calculate cluster metrics for each k
icl <- calcICL(cc_results, title = cc_dir, plot = "pdf")

saveRDS(cc_results, file.path(results_dir, "consensus_clustering_results.rds"))
saveRDS(icl, file.path(results_dir, "consensus_clustering_ICL.rds"))

# ============================================================================
# PART 5: Determine optimal k and assign subtypes
# ============================================================================
cat("\n=== Evaluating optimal k ===\n")

# Calculate PAC (Proportion of Ambiguous Clustering) for k=2:6
pac_scores <- sapply(2:6, function(k) {
  cons_mat <- cc_results[[k]]$consensusMatrix
  cdf <- ecdf(cons_mat[lower.tri(cons_mat)])
  pac <- cdf(0.9) - cdf(0.1)
  return(pac)
})
names(pac_scores) <- paste0("k=", 2:6)
cat("PAC scores:\n")
print(round(pac_scores, 4))

# Silhouette width
library(cluster)
sil_scores <- sapply(2:6, function(k) {
  clusters <- cc_results[[k]]$consensusClass
  cons_mat <- cc_results[[k]]$consensusMatrix
  dist_mat <- as.dist(1 - cons_mat)
  sil <- silhouette(clusters, dist_mat)
  mean(sil[, "sil_width"])
})
names(sil_scores) <- paste0("k=", 2:6)
cat("Mean silhouette scores:\n")
print(round(sil_scores, 4))

# Select optimal k (lowest PAC)
optimal_k <- which.min(pac_scores) + 1
cat("\nOptimal k by PAC:", optimal_k, "\n")

# Assign subtypes
subtype_assignments <- cc_results[[optimal_k]]$consensusClass
names(subtype_assignments) <- colnames(elderly_expr)

elderly_samples$subtype <- paste0("C", subtype_assignments)
saveRDS(elderly_samples, file.path(data_dir, "clinical", "elderly_sample_info_with_subtypes.rds"))

# ============================================================================
# PART 6: Subtype heatmap (Fig 1)
# ============================================================================
cat("\n=== Generating subtype heatmap ===\n")

# Top annotation
subtype_colors <- c("C1" = "#E41A1C", "C2" = "#377EB8", "C3" = "#4DAF4A",
                     "C4" = "#984EA3", "C5" = "#FF7F00", "C6" = "#A65628")
subtype_colors <- subtype_colors[1:optimal_k]

ha_top <- HeatmapAnnotation(
  Subtype = elderly_samples$subtype,
  Age = elderly_samples$age_years,
  Gender = elderly_samples$gender,
  col = list(
    Subtype = subtype_colors,
    Age = colorRamp2(c(60, 75, 90), c("lightyellow", "orange", "red")),
    Gender = c("female" = "#FF69B4", "male" = "#4169E1")
  ),
  annotation_name_side = "left"
)

# Order samples by subtype
sample_order <- order(subtype_assignments)

pdf(file.path(results_dir, "figures", "Fig1_subtype_heatmap.pdf"), width = 14, height = 10)
ht <- Heatmap(
  expr_scaled[1:500, sample_order],
  name = "Z-score",
  top_annotation = ha_top[sample_order],
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_split = factor(elderly_samples$subtype[sample_order]),
  use_raster = TRUE,
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
)
draw(ht)
dev.off()
cat("Fig1 heatmap saved.\n")

# ============================================================================
# PART 7: Subtype summary stats
# ============================================================================
cat("\n=== Subtype summary ===\n")
subtype_summary <- elderly_samples %>%
  group_by(subtype) %>%
  summarise(
    n = n(),
    age_mean = round(mean(age_years, na.rm = TRUE), 1),
    female_pct = round(100 * sum(gender == "female", na.rm = TRUE) / n(), 1),
    .groups = "drop"
  )
print(subtype_summary)
write.csv(subtype_summary, file.path(results_dir, "tables", "subtype_summary.csv"), row.names = FALSE)

cat("\n=== Preprocessing and subtyping complete ===\n")
cat("Optimal k:", optimal_k, "\n")
cat("Subtype sizes:", table(elderly_samples$subtype), "\n")
