library(SummarizedExperiment)

base_dir <- "/data2/projects/Carcinoma_of_thyroid"
data_dir <- file.path(base_dir, "data/tcga_thca")

# Extract clinical from SE object (most reliable)
cat("=== Extracting clinical from SummarizedExperiment ===\n")
se <- readRDS(file.path(data_dir, "rnaseq", "TCGA_THCA_rnaseq_SE.rds"))
sample_info <- as.data.frame(colData(se))
saveRDS(sample_info, file.path(data_dir, "clinical", "TCGA_THCA_sample_info.rds"))
write.csv(sample_info, file.path(data_dir, "clinical", "TCGA_THCA_sample_info.csv"), row.names = FALSE)
cat("Sample info saved. Samples:", nrow(sample_info), "\n")
cat("Columns:", paste(colnames(sample_info)[1:30], collapse = ", "), "\n")

# Use sample_info as clinical (it contains age, gender, stage, etc.)
# Create a deduplicated clinical table by patient
clinical <- sample_info[!duplicated(sample_info$patient), ]
saveRDS(clinical, file.path(data_dir, "clinical", "TCGA_THCA_clinical.rds"))
write.csv(clinical, file.path(data_dir, "clinical", "TCGA_THCA_clinical.csv"), row.names = FALSE)
cat("Clinical saved. Patients:", nrow(clinical), "\n")

# Check age field
cat("Age field:", head(sample_info$age_at_index, 10), "\n")
cat("Gender:", table(sample_info$gender), "\n")
cat("Stage:", table(sample_info$ajcc_pathologic_stage), "\n")

cat("\n=== Clinical data ready ===\n")
