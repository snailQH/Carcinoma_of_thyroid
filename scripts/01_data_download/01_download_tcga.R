#!/usr/bin/env Rscript
# =============================================================================
# 01_download_tcga.R — Download TCGA-THCA multi-omics data
# =============================================================================

library(TCGAbiolinks)
library(SummarizedExperiment)

base_dir <- "/data2/projects/Carcinoma_of_thyroid"
data_dir <- file.path(base_dir, "data/tcga_thca")

# ---------------------------------------------------------------------------
# 1. RNA-seq (HTSeq counts + TPM)
# ---------------------------------------------------------------------------
cat("=== Downloading RNA-seq data ===\n")
query_rnaseq <- GDCquery(
  project = "TCGA-THCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query_rnaseq, directory = file.path(data_dir, "rnaseq"))
se_rnaseq <- GDCprepare(query_rnaseq, directory = file.path(data_dir, "rnaseq"))

# Extract counts and TPM
counts_mat <- assay(se_rnaseq, "unstranded")
tpm_mat <- assay(se_rnaseq, "tpm_unstrand")
gene_info <- as.data.frame(rowData(se_rnaseq))
sample_info <- as.data.frame(colData(se_rnaseq))

saveRDS(se_rnaseq, file.path(data_dir, "rnaseq", "TCGA_THCA_rnaseq_SE.rds"))
saveRDS(counts_mat, file.path(data_dir, "rnaseq", "TCGA_THCA_counts.rds"))
saveRDS(tpm_mat, file.path(data_dir, "rnaseq", "TCGA_THCA_tpm.rds"))
saveRDS(gene_info, file.path(data_dir, "rnaseq", "gene_info.rds"))
cat("RNA-seq done. Samples:", ncol(counts_mat), "Genes:", nrow(counts_mat), "\n")

# ---------------------------------------------------------------------------
# 2. Somatic mutations (MAF)
# ---------------------------------------------------------------------------
cat("=== Downloading mutation data ===\n")
query_maf <- GDCquery(
  project = "TCGA-THCA",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)
GDCdownload(query_maf, directory = file.path(data_dir, "mutation"))
maf_data <- GDCprepare(query_maf, directory = file.path(data_dir, "mutation"))
saveRDS(maf_data, file.path(data_dir, "mutation", "TCGA_THCA_maf.rds"))
cat("Mutation data done. Rows:", nrow(maf_data), "\n")

# ---------------------------------------------------------------------------
# 3. Copy Number Variation
# ---------------------------------------------------------------------------
cat("=== Downloading CNV data ===\n")
query_cnv <- GDCquery(
  project = "TCGA-THCA",
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number"
)
GDCdownload(query_cnv, directory = file.path(data_dir, "cnv"))
cnv_data <- GDCprepare(query_cnv, directory = file.path(data_dir, "cnv"))
saveRDS(cnv_data, file.path(data_dir, "cnv", "TCGA_THCA_cnv.rds"))
cat("CNV data done.\n")

# ---------------------------------------------------------------------------
# 4. Clinical data
# ---------------------------------------------------------------------------
cat("=== Downloading clinical data ===\n")
clinical <- GDCquery_clinic("TCGA-THCA", type = "clinical")
saveRDS(clinical, file.path(data_dir, "clinical", "TCGA_THCA_clinical.rds"))
write.csv(clinical, file.path(data_dir, "clinical", "TCGA_THCA_clinical.csv"), row.names = FALSE)
cat("Clinical data done. Patients:", nrow(clinical), "\n")

# Save sample info from SE object (has age, stage, etc.)
saveRDS(sample_info, file.path(data_dir, "clinical", "TCGA_THCA_sample_info.rds"))
write.csv(sample_info, file.path(data_dir, "clinical", "TCGA_THCA_sample_info.csv"), row.names = FALSE)

cat("\n=== All TCGA-THCA data downloaded successfully ===\n")
