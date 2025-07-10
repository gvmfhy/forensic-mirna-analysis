#!/usr/bin/env Rscript

# Minimal CEL processor - assumes packages are installed
# If not installed, run:
# BiocManager::install(c("oligo", "limma"))

cat("=== Minimal CEL Processing ===\n")
cat("Loading libraries...\n")

library(oligo)
library(limma)

# Set paths
cel_dir <- "data/raw/GSE49630_CEL"
output_dir <- "data/processed/cel"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# List CEL files
cel_files <- list.files(cel_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE)
cat(sprintf("Found %d CEL files\n", length(cel_files)))

# Extract sample info from filenames
sample_names <- gsub(".CEL.gz", "", basename(cel_files))
fluid_types <- character(length(sample_names))
fluid_types[grep("Blood_", sample_names)] <- "blood"
fluid_types[grep("Semen_", sample_names)] <- "semen"
fluid_types[grep("Vaginal_", sample_names)] <- "vaginal"
fluid_types[grep("Saliva_", sample_names)] <- "saliva"

cat("\nSample distribution:\n")
print(table(fluid_types))

# Decompress files for oligo
cat("\nDecompressing CEL files...\n")
temp_dir <- tempdir()
temp_files <- character(length(cel_files))

for (i in 1:length(cel_files)) {
  temp_files[i] <- file.path(temp_dir, gsub(".gz", "", basename(cel_files[i])))
  if (!file.exists(temp_files[i])) {
    system2("gunzip", args = c("-c", cel_files[i]), stdout = temp_files[i])
  }
}

# Read CEL files
cat("\nReading CEL files with oligo...\n")
raw_data <- read.celfiles(temp_files)

# RMA normalization
cat("\nPerforming RMA normalization...\n")
eset <- rma(raw_data)

# Extract expression matrix
expr_matrix <- exprs(eset)
colnames(expr_matrix) <- sample_names

cat(sprintf("\nExpression matrix: %d probes x %d samples\n", 
            nrow(expr_matrix), ncol(expr_matrix)))

# Save expression matrix
write.csv(expr_matrix, file.path(output_dir, "cel_expression_matrix.csv"))

# Save metadata
metadata <- data.frame(
  sample = sample_names,
  fluid_type = fluid_types,
  stringsAsFactors = FALSE
)
write.csv(metadata, file.path(output_dir, "cel_metadata.csv"), row.names = FALSE)

# Clean up temp files
unlink(temp_files)

cat("\nProcessing complete!\n")
cat(sprintf("Output saved to: %s\n", output_dir))