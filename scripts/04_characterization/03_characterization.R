#!/usr/bin/env Rscript
# =============================================================================
# 03_characterization.R — Aim 2: Subtype characterization (fixed)
# =============================================================================

library(GSVA)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
library(dplyr)

base_dir <- "/data2/projects/Carcinoma_of_thyroid"
data_dir <- file.path(base_dir, "data/tcga_thca")
results_dir <- file.path(base_dir, "results")

cat("=== Loading data ===\n")
elderly_expr <- readRDS(file.path(data_dir, "rnaseq", "elderly_log2tpm_filtered.rds"))
elderly_info <- readRDS(file.path(data_dir, "clinical", "elderly_sample_info_with_subtypes.rds"))
maf_data <- readRDS(file.path(data_dir, "mutation", "TCGA_THCA_maf.rds"))

subtypes <- elderly_info$subtype
gene_names <- rownames(elderly_expr)
cat("Samples:", ncol(elderly_expr), "Genes:", length(gene_names), "\n")
cat("Subtypes:", table(subtypes), "\n")

# ============================================================================
# 2A: Mutation Landscape
# ============================================================================
cat("\n=== 2A: Mutation landscape ===\n")

elderly_patients <- substr(colnames(elderly_expr), 1, 12)
# Match MAF barcodes
maf_barcodes <- unique(substr(maf_data$Tumor_Sample_Barcode, 1, 12))
overlap <- intersect(elderly_patients, maf_barcodes)
cat("Patients with mutation data:", length(overlap), "\n")

if (length(overlap) > 0) {
  # Key driver genes
  drivers <- c("BRAF", "NRAS", "HRAS", "KRAS", "RET", "TP53", "PIK3CA", "PTEN", "AKT1", "TERT")
  driver_muts <- maf_data[maf_data$Hugo_Symbol %in% drivers &
                           substr(maf_data$Tumor_Sample_Barcode, 1, 12) %in% overlap, ]

  # Mutation frequency per subtype
  mut_freq <- data.frame()
  for (sub in unique(subtypes)) {
    sub_patients <- elderly_patients[subtypes == sub]
    sub_overlap <- intersect(sub_patients, maf_barcodes)
    for (gene in drivers) {
      n_mut <- length(unique(driver_muts$Tumor_Sample_Barcode[
        driver_muts$Hugo_Symbol == gene & substr(driver_muts$Tumor_Sample_Barcode, 1, 12) %in% sub_overlap]))
      mut_freq <- rbind(mut_freq, data.frame(
        Subtype = sub, Gene = gene,
        Mutated = n_mut, Total = length(sub_overlap),
        Frequency = round(n_mut / max(length(sub_overlap), 1) * 100, 1)
      ))
    }
  }
  write.csv(mut_freq, file.path(results_dir, "tables", "driver_mutation_frequency.csv"), row.names = FALSE)

  # TMB per subtype
  tmb_df <- maf_data %>%
    filter(substr(Tumor_Sample_Barcode, 1, 12) %in% overlap) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(n_mutations = n(), .groups = "drop") %>%
    mutate(patient = substr(Tumor_Sample_Barcode, 1, 12),
           TMB = n_mutations / 30) # ~30 Mb exome

  tmb_df$Subtype <- subtypes[match(tmb_df$patient, elderly_patients)]
  tmb_df <- tmb_df[!is.na(tmb_df$Subtype), ]

  p_tmb <- ggboxplot(tmb_df, x = "Subtype", y = "TMB", fill = "Subtype", palette = "jco") +
    stat_compare_means() +
    labs(y = "TMB (mutations/Mb)", title = "Tumor Mutational Burden by Subtype") +
    theme_minimal()
  ggsave(file.path(results_dir, "figures", "Fig2_TMB_by_subtype.pdf"), p_tmb, width = 7, height = 5)
  cat("TMB analysis done.\n")
}

# ============================================================================
# 2C: GSVA Pathway Scoring
# ============================================================================
cat("\n=== 2C: GSVA pathway scoring ===\n")

gene_sets <- list(
  EMT = c("VIM", "FN1", "CDH2", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", "MMP2", "MMP9", "COL1A1", "ACTA2"),
  Proliferation = c("MKI67", "PCNA", "TOP2A", "CDK1", "CCNB1", "CCND1", "CCNE1", "MCM2", "AURKA", "PLK1"),
  Angiogenesis = c("VEGFA", "VEGFB", "VEGFC", "FLT1", "KDR", "PECAM1", "CDH5", "ANGPT1", "ANGPT2", "HIF1A"),
  Adhesion = c("CDH1", "ITGB1", "ITGB3", "ITGA5", "ICAM1", "VCAM1", "CD44", "EPCAM", "CLDN1"),
  Glycolysis = c("HK1", "HK2", "PFKFB3", "ALDOA", "GAPDH", "PGK1", "ENO1", "PKM", "LDHA", "SLC2A1"),
  OXPHOS = c("NDUFA1", "SDHA", "SDHB", "UQCRC1", "COX4I1", "ATP5F1A", "CS", "IDH2"),
  SASP = c("IL6", "CXCL1", "CXCL2", "CCL2", "CCL5", "MMP1", "MMP3", "SERPINE1", "IGFBP3", "VEGFA"),
  Senescence = c("CDKN1A", "CDKN2A", "CDKN2B", "TP53", "RB1", "SERPINE1", "LMNB1", "HMGA1"),
  Immunosenescence = c("KLRG1", "TIGIT", "LAG3", "HAVCR2", "PDCD1", "TOX", "CTLA4", "CD244"),
  DNA_Repair = c("BRCA1", "BRCA2", "RAD51", "ATM", "ATR", "CHEK1", "CHEK2", "MLH1", "MSH2"),
  Apoptosis = c("BCL2", "BCL2L1", "BAX", "BAK1", "CASP3", "CASP8", "CASP9", "FAS")
)

# Filter to present genes
gs_filtered <- lapply(gene_sets, function(g) intersect(g, gene_names))
gs_filtered <- gs_filtered[sapply(gs_filtered, length) >= 3]

cat("Running GSVA on", length(gs_filtered), "gene sets...\n")
gsva_params <- gsvaParam(as.matrix(elderly_expr), gs_filtered, maxDiff = TRUE)
gsva_scores <- gsva(gsva_params)
saveRDS(gsva_scores, file.path(results_dir, "gsva_scores.rds"))

# Heatmap
ord <- order(subtypes)
sub_colors <- c("C1" = "#E41A1C", "C2" = "#377EB8", "C3" = "#4DAF4A")[1:length(unique(subtypes))]
ha <- HeatmapAnnotation(Subtype = subtypes[ord], col = list(Subtype = sub_colors))

pdf(file.path(results_dir, "figures", "Fig4_GSVA_heatmap.pdf"), width = 14, height = 8)
Heatmap(gsva_scores[, ord], name = "GSVA", top_annotation = ha,
        cluster_columns = FALSE, column_split = factor(subtypes[ord]),
        show_column_names = FALSE, col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")))
dev.off()

# Boxplots
gsva_df <- as.data.frame(t(gsva_scores))
gsva_df$Subtype <- subtypes
for (pw in rownames(gsva_scores)) {
  p <- ggboxplot(gsva_df, x = "Subtype", y = pw, fill = "Subtype", palette = "jco") +
    stat_compare_means() + labs(title = paste(pw, "by Subtype")) + theme_minimal()
  ggsave(file.path(results_dir, "figures", paste0("GSVA_", pw, "_boxplot.pdf")), p, width = 7, height = 5)
}
cat("GSVA done.\n")

# ============================================================================
# 2D: Thyroid Differentiation Score
# ============================================================================
cat("\n=== 2D: TDS ===\n")
tds_genes <- c("TG", "TPO", "SLC5A5", "SLC26A4", "DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", "PAX8", "TSHR")
tds_present <- intersect(tds_genes, gene_names)
if (length(tds_present) >= 3) {
  tds <- colMeans(elderly_expr[tds_present, ])
  p_tds <- ggboxplot(data.frame(Subtype = subtypes, TDS = tds), x = "Subtype", y = "TDS",
                      fill = "Subtype", palette = "jco") +
    stat_compare_means() + labs(title = "Thyroid Differentiation Score", y = "TDS") + theme_minimal()
  ggsave(file.path(results_dir, "figures", "Fig4_TDS_by_subtype.pdf"), p_tds, width = 7, height = 5)
  cat("TDS genes used:", paste(tds_present, collapse = ", "), "\n")
}

# ============================================================================
# 2E: Immune Deconvolution
# ============================================================================
cat("\n=== 2E: Immune scoring ===\n")
immune_markers <- list(
  CD8_T = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG"),
  Treg = c("FOXP3", "IL2RA", "CTLA4", "IKZF2"),
  M2_Macrophage = c("CD163", "MRC1", "MSR1", "CD68", "IL10"),
  M1_Macrophage = c("NOS2", "CD80", "CD86", "TNF", "IL1B"),
  NK = c("NCR1", "KLRK1", "KLRD1", "NKG7"),
  B_cell = c("CD19", "MS4A1", "CD79A", "PAX5"),
  CAF = c("FAP", "ACTA2", "PDGFRA", "PDGFRB", "COL1A1", "THY1"),
  MDSC = c("S100A8", "S100A9", "ARG1")
)

immune_df <- data.frame(row.names = colnames(elderly_expr))
for (ct in names(immune_markers)) {
  present <- intersect(immune_markers[[ct]], gene_names)
  if (length(present) >= 2) immune_df[[ct]] <- colMeans(elderly_expr[present, , drop = FALSE])
}
immune_df$Subtype <- subtypes
saveRDS(immune_df, file.path(results_dir, "immune_scores.rds"))

# Immune heatmap
immune_mat <- as.matrix(immune_df[, -ncol(immune_df)])
immune_means <- aggregate(immune_mat, by = list(Subtype = immune_df$Subtype), FUN = mean)
rownames(immune_means) <- immune_means$Subtype
immune_means$Subtype <- NULL

pdf(file.path(results_dir, "figures", "Fig3_immune_landscape.pdf"), width = 10, height = 6)
Heatmap(as.matrix(scale(t(immune_means))), name = "Z-score",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
dev.off()

# Immune boxplots
for (ct in names(immune_markers)) {
  if (ct %in% colnames(immune_df)) {
    p <- ggboxplot(immune_df, x = "Subtype", y = ct, fill = "Subtype", palette = "jco") +
      stat_compare_means() + labs(title = paste(ct, "Score by Subtype")) + theme_minimal()
    ggsave(file.path(results_dir, "figures", paste0("Immune_", ct, "_boxplot.pdf")), p, width = 7, height = 5)
  }
}

# Checkpoint expression
checkpoints <- c("CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "CD276")
cp_present <- intersect(checkpoints, gene_names)
if (length(cp_present) >= 2) {
  cp_df <- as.data.frame(t(elderly_expr[cp_present, ]))
  cp_df$Subtype <- subtypes
  cp_means <- aggregate(. ~ Subtype, cp_df, mean)
  rownames(cp_means) <- cp_means$Subtype
  cp_means$Subtype <- NULL

  pdf(file.path(results_dir, "figures", "Fig3_checkpoint_heatmap.pdf"), width = 8, height = 5)
  Heatmap(as.matrix(t(cp_means)), name = "log2TPM",
          col = colorRamp2(c(0, 3, 6), c("white", "orange", "red")))
  dev.off()
}
cat("Immune analysis done.\n")

# ============================================================================
# 2G: Clinical Association
# ============================================================================
cat("\n=== 2G: Clinical association ===\n")

# Stage by subtype
if ("ajcc_pathologic_stage" %in% colnames(elderly_info)) {
  stage_tab <- table(elderly_info$subtype, elderly_info$ajcc_pathologic_stage)
  cat("Stage distribution:\n")
  print(stage_tab)
  test <- fisher.test(stage_tab, simulate.p.value = TRUE)
  cat("Fisher p-value:", test$p.value, "\n")
}

# LN metastasis
if ("ajcc_pathologic_n" %in% colnames(elderly_info)) {
  elderly_info$LN_met <- ifelse(grepl("N1", elderly_info$ajcc_pathologic_n), "N1", "N0")
  ln_tab <- table(elderly_info$subtype, elderly_info$LN_met)
  cat("LN metastasis:\n")
  print(ln_tab)

  p_ln <- ggplot(elderly_info, aes(x = subtype, fill = LN_met)) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = c("N0" = "#4DAF4A", "N1" = "#E41A1C")) +
    labs(title = "LN Metastasis by Subtype", y = "Proportion") + theme_minimal()
  ggsave(file.path(results_dir, "figures", "Fig6_LN_metastasis.pdf"), p_ln, width = 7, height = 5)
}

# Survival
elderly_info$status <- ifelse(elderly_info$vital_status == "Dead", 1, 0)
elderly_info$time <- as.numeric(ifelse(elderly_info$status == 1,
                                        elderly_info$days_to_death,
                                        elderly_info$days_to_last_follow_up))
surv_data <- elderly_info[!is.na(elderly_info$time) & elderly_info$time > 0, ]

if (nrow(surv_data) > 10) {
  library(survival)
  library(survminer)
  fit <- survfit(Surv(time, status) ~ subtype, data = surv_data)
  pdf(file.path(results_dir, "figures", "Fig6_survival.pdf"), width = 10, height = 8)
  print(ggsurvplot(fit, data = surv_data, pval = TRUE, risk.table = TRUE, palette = "jco",
                    title = "OS by Subtype (Elderly)"))
  dev.off()
  cat("Survival plot done.\n")
}

# Save clinical summary
clin_summary <- elderly_info %>%
  group_by(subtype) %>%
  summarise(
    n = n(),
    age_mean = round(mean(age_years, na.rm = TRUE), 1),
    female_pct = round(100 * sum(gender == "female", na.rm = TRUE) / n(), 1),
    N1_pct = round(100 * sum(grepl("N1", ajcc_pathologic_n), na.rm = TRUE) / n(), 1),
    .groups = "drop"
  )
write.csv(clin_summary, file.path(results_dir, "tables", "subtype_clinical_summary.csv"), row.names = FALSE)
print(clin_summary)

cat("\n=== Characterization COMPLETE ===\n")
