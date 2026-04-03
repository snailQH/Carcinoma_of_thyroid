library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

base_dir <- "/data2/projects/Carcinoma_of_thyroid"
data_dir <- file.path(base_dir, "data/tcga_thca")

# ---- 1. Fix mutation: read MAF files manually ----
cat("=== Processing mutation data ===\n")
maf_dir <- file.path(data_dir, "mutation")
maf_files <- list.files(maf_dir, pattern = "\\.maf\\.gz$", recursive = TRUE, full.names = TRUE)

if (length(maf_files) == 0) {
  # Re-download mutation data
  query_maf <- GDCquery(
    project = "TCGA-THCA",
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
  )
  GDCdownload(query_maf, directory = maf_dir)
  maf_files <- list.files(maf_dir, pattern = "\\.maf\\.gz$", recursive = TRUE, full.names = TRUE)
}

cat("Found", length(maf_files), "MAF files\n")

# Read all MAF files, coercing all columns to character first
all_mafs <- lapply(maf_files, function(f) {
  tryCatch({
    df <- read.delim(gzfile(f), comment.char = "#", stringsAsFactors = FALSE)
    # Force all columns to character to avoid type conflicts
    df[] <- lapply(df, as.character)
    df
  }, error = function(e) {
    cat("  Error reading", basename(f), ":", e$message, "\n")
    NULL
  })
})
all_mafs <- all_mafs[!sapply(all_mafs, is.null)]
maf_combined <- bind_rows(all_mafs)
cat("Combined MAF rows:", nrow(maf_combined), "\n")

saveRDS(maf_combined, file.path(maf_dir, "TCGA_THCA_maf.rds"))
cat("Mutation data saved.\n")

# ---- 2. Clinical data ----
cat("\n=== Downloading clinical data ===\n")
clinical <- GDCquery_clinic("TCGA-THCA", type = "clinical")
saveRDS(clinical, file.path(data_dir, "clinical", "TCGA_THCA_clinical.rds"))
write.csv(clinical, file.path(data_dir, "clinical", "TCGA_THCA_clinical.csv"), row.names = FALSE)
cat("Clinical done. Patients:", nrow(clinical), "\n")

# ---- 3. Sample info from SE ----
cat("\n=== Extracting sample info ===\n")
se <- readRDS(file.path(data_dir, "rnaseq", "TCGA_THCA_rnaseq_SE.rds"))
sample_info <- as.data.frame(colData(se))
saveRDS(sample_info, file.path(data_dir, "clinical", "TCGA_THCA_sample_info.rds"))
write.csv(sample_info, file.path(data_dir, "clinical", "TCGA_THCA_sample_info.csv"), row.names = FALSE)
cat("Sample info saved. Samples:", nrow(sample_info), "\n")

# ---- 4. CNV ----
cat("\n=== Downloading CNV data ===\n")
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
