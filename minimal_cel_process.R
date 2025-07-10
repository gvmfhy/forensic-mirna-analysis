#!/usr/bin/env Rscript

# Ultra-minimal CEL processing using only affy (which we have installed)
cat("=== Minimal CEL Processing with affy ===\n")

# Load only what we have
library(affy)
library(limma)

# Paths
cel_dir <- "data/raw/GSE49630_CEL"
output_dir <- "data/processed/cel"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# List CEL files
cel_files <- list.files(cel_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE)
cat(sprintf("Found %d CEL files\n", length(cel_files)))

# Since affy can't read compressed CEL files, decompress them
cat("\nDecompressing CEL files...\n")
temp_dir <- tempdir()
temp_files <- character(length(cel_files))

for (i in 1:length(cel_files)) {
  temp_files[i] <- file.path(temp_dir, gsub(".gz", "", basename(cel_files[i])))
  if (!file.exists(temp_files[i])) {
    system2("gunzip", args = c("-c", cel_files[i]), stdout = temp_files[i])
  }
}

# Try to read with affy
cat("\nAttempting to read CEL files with affy...\n")
tryCatch({
  # Read CEL files
  raw_data <- ReadAffy(filenames = temp_files)
  
  # RMA normalization
  cat("\nPerforming RMA normalization...\n")
  eset <- rma(raw_data)
  
  # Extract expression matrix
  expr_matrix <- exprs(eset)
  
  # Fix column names
  sample_names <- gsub(".CEL", "", basename(temp_files))
  colnames(expr_matrix) <- sample_names
  
  cat(sprintf("\nExpression matrix: %d probes x %d samples\n", 
              nrow(expr_matrix), ncol(expr_matrix)))
  
  # Save minimal output
  write.csv(expr_matrix, file.path(output_dir, "cel_expression_matrix_affy.csv"))
  
  # Create metadata
  fluid_types <- character(length(sample_names))
  fluid_types[grep("Blood_", sample_names)] <- "blood"
  fluid_types[grep("Semen_", sample_names)] <- "semen"
  fluid_types[grep("Vaginal_", sample_names)] <- "vaginal"
  fluid_types[grep("Saliva_", sample_names)] <- "saliva"
  
  metadata <- data.frame(
    sample = sample_names,
    fluid_type = fluid_types,
    stringsAsFactors = FALSE
  )
  
  write.csv(metadata, file.path(output_dir, "cel_metadata.csv"), row.names = FALSE)
  
  cat("\nSUCCESS! Basic processing complete.\n")
  cat(sprintf("Output saved to: %s\n", output_dir))
  
}, error = function(e) {
  cat("\nERROR: affy package cannot read these CEL files\n")
  cat("These appear to be newer Affymetrix ST arrays that require oligo\n")
  cat("Error details:", conditionMessage(e), "\n")
  
  # Alternative: create a file listing what needs to be done
  cat("\nCreating instructions for manual processing...\n")
  
  instructions <- c(
    "# Manual CEL Processing Instructions",
    "",
    "These Affymetrix miRNA ST arrays require the 'oligo' package.",
    "Since automatic installation is failing, here are manual steps:",
    "",
    "## Option 1: Use a different R environment",
    "1. Install R from CRAN (not homebrew)",
    "2. Install BiocManager: install.packages('BiocManager')",
    "3. Install oligo: BiocManager::install('oligo')",
    "4. Run: Rscript scripts/preprocessing/cel_processor_minimal.R",
    "",
    "## Option 2: Use Docker",
    "```bash",
    "docker run -it -v $(pwd):/data bioconductor/bioconductor_docker:latest R",
    "# Inside container:",
    "BiocManager::install('oligo')",
    "source('/data/scripts/preprocessing/cel_processor_minimal.R')",
    "```",
    "",
    "## Option 3: Process on another machine",
    "Transfer the CEL files to a machine with working Bioconductor",
    "",
    "## What the data contains:",
    paste("- ", length(cel_files), "CEL files (", 
          length(grep("Blood", cel_files)), "blood,",
          length(grep("Saliva", cel_files)), "saliva,",
          length(grep("Semen", cel_files)), "semen,",
          length(grep("Vaginal", cel_files)), "vaginal)"),
    "- Platform: Affymetrix miRNA 3.0/4.0 ST arrays",
    "- These are individual samples (not pooled like GPR)",
    ""
  )
  
  writeLines(instructions, file.path(output_dir, "CEL_PROCESSING_INSTRUCTIONS.md"))
  cat("Instructions saved to:", file.path(output_dir, "CEL_PROCESSING_INSTRUCTIONS.md"), "\n")
})

# Clean up
unlink(temp_files)