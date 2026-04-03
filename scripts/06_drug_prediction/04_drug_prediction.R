#!/usr/bin/env Rscript
# =============================================================================
# 04_drug_prediction.R — Aim 4: Drug sensitivity & therapeutic target prediction
# =============================================================================

library(ggplot2)
library(ggpubr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

base_dir <- "/data2/projects/Carcinoma_of_thyroid"
data_dir <- file.path(base_dir, "data/tcga_thca")
results_dir <- file.path(base_dir, "results")

# Load data
cat("=== Loading data ===\n")
elderly_expr <- readRDS(file.path(data_dir, "rnaseq", "elderly_log2tpm_filtered.rds"))
elderly_info <- readRDS(file.path(data_dir, "clinical", "elderly_sample_info_with_subtypes.rds"))

# ============================================================================
# 4A: Drug sensitivity prediction using oncoPredict
# ============================================================================
cat("\n=== 4A: Drug sensitivity prediction ===\n")

# Try oncoPredict; if not available, use pRRophetic or manual ridge regression
tryCatch({
  if (requireNamespace("oncoPredict", quietly = TRUE)) {
    library(oncoPredict)
    cat("Using oncoPredict package.\n")

    # GDSC2 training data (need to download)
    # oncoPredict provides built-in training data
    calcPhenotype(
      trainingExprData = GDSC2_Expr,
      trainingPtype = GDSC2_Res,
      testExprData = as.matrix(elderly_expr),
      batchCorrect = "eb",
      powerTransformPhenotype = TRUE,
      removeLowVaryingGenes = 0.2,
      minNumSamples = 10,
      selection = -1,
      printOutput = TRUE,
      removeLowVaringGenesFrom = "homogenizeData",
      report_pc = FALSE,
      cc = FALSE,
      rsq = FALSE
    )

    # Read results
    drug_pred <- read.csv(file.path(getwd(), "calcPhenotype_Output/DrugPredictions.csv"),
                          row.names = 1)
    cat("Drug predictions generated for", ncol(drug_pred), "drugs.\n")
  }
}, error = function(e) {
  cat("oncoPredict not available, using manual approach.\n")
})

# ============================================================================
# Manual drug sensitivity scoring (fallback / complementary)
# ============================================================================
cat("\n=== Manual drug sensitivity scoring ===\n")

# Drug target gene sets (expression-based proxy for sensitivity)
drug_target_sets <- list(
  BRAF_inhibitor_sensitivity = c("BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3"),
  MEK_inhibitor_sensitivity = c("MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "DUSP6"),
  Lenvatinib_targets = c("FLT1", "KDR", "FLT4", "FGFR1", "FGFR2", "FGFR3", "FGFR4",
                          "PDGFRA", "RET", "KIT"),
  Sorafenib_targets = c("BRAF", "RAF1", "FLT3", "KIT", "PDGFRB", "KDR", "FLT4", "RET"),
  CDK_inhibitor_sensitivity = c("CDK4", "CDK6", "CCND1", "RB1", "CDKN2A"),
  mTOR_inhibitor_sensitivity = c("MTOR", "RPTOR", "RICTOR", "AKT1", "PIK3CA", "PTEN"),
  Anti_PD1_response_markers = c("CD274", "PDCD1LG2", "IFNG", "CXCL9", "CXCL10",
                                 "GZMA", "GZMB", "PRF1", "CD8A", "STAT1"),
  RAI_sensitivity = c("SLC5A5", "TPO", "TG", "TSHR", "DIO1", "DIO2", "DUOX2", "PAX8")
)

gene_names <- rownames(elderly_expr)
drug_scores <- data.frame(row.names = colnames(elderly_expr))

for (drug in names(drug_target_sets)) {
  targets_present <- intersect(drug_target_sets[[drug]], gene_names)
  if (length(targets_present) >= 2) {
    drug_scores[[drug]] <- colMeans(elderly_expr[targets_present, , drop = FALSE])
  }
}
drug_scores$Subtype <- elderly_info$subtype

saveRDS(drug_scores, file.path(results_dir, "drug_sensitivity_scores.rds"))

# ============================================================================
# Drug sensitivity heatmap (Fig 7)
# ============================================================================
cat("\n=== Generating drug sensitivity heatmap ===\n")

drug_mat <- as.matrix(drug_scores[, -ncol(drug_scores)])
drug_means <- aggregate(drug_mat, by = list(Subtype = drug_scores$Subtype), FUN = mean)
rownames(drug_means) <- drug_means$Subtype
drug_means$Subtype <- NULL

pdf(file.path(results_dir, "figures", "Fig7_drug_sensitivity_heatmap.pdf"), width = 10, height = 6)
Heatmap(
  as.matrix(scale(t(drug_means))),
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 12),
  row_labels = gsub("_", " ", rownames(t(drug_means)))
)
dev.off()

# Boxplots per drug by subtype
for (drug in colnames(drug_mat)) {
  p <- ggboxplot(drug_scores, x = "Subtype", y = drug,
                  fill = "Subtype", palette = "jco") +
    stat_compare_means() +
    labs(title = gsub("_", " ", drug), y = "Score (mean log2 TPM)") +
    theme_minimal()
  ggsave(file.path(results_dir, "figures", paste0("Drug_", drug, "_boxplot.pdf")),
         p, width = 7, height = 5)
}

# ============================================================================
# 4C: TIDE immunotherapy prediction
# ============================================================================
cat("\n=== 4C: Immunotherapy prediction markers ===\n")

# TIDE-like scoring (simplified)
# TIDE = T cell dysfunction + exclusion
tide_dysfunction_genes <- c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT",
                             "TOX", "EOMES", "BATF")
tide_exclusion_genes <- c("TGFB1", "TGFB2", "TGFB3", "COL1A1", "COL1A2",
                           "FAP", "ACTA2", "FN1", "VIM")

dys_present <- intersect(tide_dysfunction_genes, gene_names)
exc_present <- intersect(tide_exclusion_genes, gene_names)

if (length(dys_present) >= 2 && length(exc_present) >= 2) {
  elderly_info$TIDE_dysfunction <- colMeans(elderly_expr[dys_present, ])
  elderly_info$TIDE_exclusion <- colMeans(elderly_expr[exc_present, ])

  p_tide <- ggplot(elderly_info, aes(x = TIDE_dysfunction, y = TIDE_exclusion, color = subtype)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "T Cell Dysfunction vs Exclusion by Subtype",
         x = "Dysfunction Score", y = "Exclusion Score") +
    theme_minimal() +
    stat_ellipse()
  ggsave(file.path(results_dir, "figures", "Fig7_TIDE_scatter.pdf"), p_tide, width = 8, height = 6)
}

# ============================================================================
# 4B: Actionable target summary table
# ============================================================================
cat("\n=== 4B: Actionable targets summary ===\n")

actionable_genes <- c("BRAF", "RET", "NTRK1", "NTRK2", "NTRK3", "ALK", "ROS1",
                       "MET", "MTOR", "PIK3CA", "PTEN", "VEGFA", "KDR",
                       "FGFR1", "FGFR2", "CDK4", "CDK6")

targets_present <- intersect(actionable_genes, gene_names)
target_expr <- as.data.frame(t(elderly_expr[targets_present, ]))
target_expr$Subtype <- elderly_info$subtype

# Mean expression per subtype
target_means <- aggregate(. ~ Subtype, target_expr, mean)
write.csv(target_means, file.path(results_dir, "tables", "actionable_target_expression.csv"), row.names = FALSE)

# Heatmap
rownames(target_means) <- target_means$Subtype
target_means$Subtype <- NULL

pdf(file.path(results_dir, "figures", "Fig7_actionable_targets.pdf"), width = 10, height = 6)
Heatmap(
  as.matrix(scale(t(target_means))),
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 10)
)
dev.off()

cat("\n=== Drug sensitivity and target prediction complete ===\n")
