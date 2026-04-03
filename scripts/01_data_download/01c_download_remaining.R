library(TCGAbiolinks)
library(SummarizedExperiment)
base_dir <- "/data2/projects/Carcinoma_of_thyroid"
data_dir <- file.path(base_dir, "data/tcga_thca")

# Mutation
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
cat("Mutation done. Rows:", nrow(maf_data), "\n")

# Clinical
cat("=== Downloading clinical data ===\n")
clinical <- GDCquery_clinic("TCGA-THCA", type = "clinical")
saveRDS(clinical, file.path(data_dir, "clinical", "TCGA_THCA_clinical.rds"))
write.csv(clinical, file.path(data_dir, "clinical", "TCGA_THCA_clinical.csv"), row.names = FALSE)
cat("Clinical done. Patients:", nrow(clinical), "\n")

# Sample info from SE
se <- readRDS(file.path(data_dir, "rnaseq", "TCGA_THCA_rnaseq_SE.rds"))
sample_info <- as.data.frame(colData(se))
saveRDS(sample_info, file.path(data_dir, "clinical", "TCGA_THCA_sample_info.rds"))
write.csv(sample_info, file.path(data_dir, "clinical", "TCGA_THCA_sample_info.csv"), row.names = FALSE)
cat("Sample info saved. Samples:", nrow(sample_info), "\n")

# CNV
cat("=== Downloading CNV data ===\n")
query_cnv <- GDCquery(
  project = "TCGA-THCA",
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number"
)
GDCdownload(query_cnv, directory = file.path(data_dir, "cnv"))
cnv_data <- GDCprepare(query_cnv, directory = file.path(data_dir, "cnv"))
saveRDS(cnv_data, file.path(data_dir, "cnv", "TCGA_THCA_cnv.rds"))
cat("CNV done.\n")
cat("\n=== All downloads complete ===\n")
