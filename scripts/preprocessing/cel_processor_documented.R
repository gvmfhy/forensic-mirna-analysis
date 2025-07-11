#!/usr/bin/env Rscript

#' Affymetrix CEL File Processor for Forensic miRNA Analysis
#'
#' This script processes Affymetrix Human Gene 1.0 ST Array CEL files from 
#' GSE49630 dataset containing individual body fluid samples for forensic
#' miRNA identification. It performs RMA normalization and prepares expression
#' data for downstream statistical analysis.
#'
#' CRITICAL: This script uses the oligo package, NOT the affy package.
#' The affy package does not support ST (Sense Target) arrays and will
#' produce invalid results. The oligo package is specifically designed
#' for newer Affymetrix platforms including ST arrays.
#'
#' WHY oligo vs affy:
#' - ST arrays use different probe design (25-mer vs 11-mer)
#' - ST arrays target exons rather than 3' UTR regions
#' - Different normalization algorithms required
#' - affy package will fail with "Could not obtain CDF environment" error
#'
#' Dataset Details (GSE49630):
#' - Platform: Affymetrix Human Gene 1.0 ST Array
#' - Samples: 20 individual donor samples (5 per body fluid)
#' - Body fluids: Blood, Saliva, Semen, Vaginal secretion
#' - Technology: Single-channel microarray (unlike two-color GPR)
#' - Probe content: ~28,000 well-annotated genes including miRNAs
#'
#' Processing Steps:
#' 1. Decompress gzipped CEL files to temporary location
#' 2. Read CEL files using oligo::read.celfiles()
#' 3. Perform RMA (Robust Multiarray Average) normalization
#' 4. Extract expression matrix and sample metadata
#' 5. Save processed data for downstream analysis
#'
#' RMA Normalization:
#' RMA performs three steps:
#' 1. Background correction using PM (Perfect Match) probes only
#' 2. Quantile normalization across arrays
#' 3. Median polish summarization across probe sets
#'
#' This produces log2-transformed, normalized expression values suitable
#' for parametric statistical analysis.
#'
#' Historical Context:
#' This processor was developed after multiple failed attempts using:
#' - affy package (incompatible with ST arrays)
#' - Custom Python implementations (overly complex)
#' - R/Python bridges (unreliable on macOS ARM64)
#'
#' The current implementation represents the most robust solution for
#' processing ST array data in the forensic miRNA analysis pipeline.
#'
#' Author: Forensic miRNA Analysis Pipeline
#' Date: 2024
#' R Version: 4.5+ required
#' Bioconductor: 3.21+ required

# Load required libraries with error handling
cat("=== Affymetrix CEL File Processor ===\n")
cat("Loading required Bioconductor packages...\n")

# Check for required packages and provide helpful error messages
required_packages <- c("oligo", "limma")
missing_packages <- character(0)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) > 0) {
  cat("ERROR: Missing required packages:\n")
  for (pkg in missing_packages) {
    cat(sprintf("  - %s\n", pkg))
  }
  cat("\nTo install missing packages, run:\n")
  cat("if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n")
  for (pkg in missing_packages) {
    cat(sprintf("BiocManager::install('%s')\n", pkg))
  }
  stop("Please install missing packages and try again")
}

# Load libraries with detailed logging
library(oligo)  # For ST array processing - CRITICAL: NOT affy
library(limma)  # For additional microarray utilities

cat(sprintf("✓ oligo version: %s\n", packageVersion("oligo")))
cat(sprintf("✓ limma version: %s\n", packageVersion("limma")))

#' Define file paths and parameters
#'
#' These paths are configured for the standard project structure.
#' Modify if your data is located elsewhere.

# Input and output directories
cel_dir <- "data/raw/GSE49630_CEL"  # Directory containing compressed CEL files
output_dir <- "data/processed/cel"   # Directory for processed output

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Validate input directory exists
if (!dir.exists(cel_dir)) {
  stop(sprintf("Input directory not found: %s\n", cel_dir))
}

#' Discover and validate CEL files
#'
#' GSE49630 CEL files are gzip-compressed with .CEL.gz extension.
#' The files follow a naming convention that includes body fluid information.

cat("\n=== File Discovery ===\n")

# Find all compressed CEL files
cel_files <- list.files(cel_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE)

if (length(cel_files) == 0) {
  stop(sprintf("No CEL.gz files found in %s\n", cel_dir))
}

cat(sprintf("Found %d CEL files:\n", length(cel_files)))

# Display file listing for verification
for (i in 1:min(5, length(cel_files))) {
  cat(sprintf("  %d. %s\n", i, basename(cel_files[i])))
}
if (length(cel_files) > 5) {
  cat(sprintf("  ... and %d more files\n", length(cel_files) - 5))
}

#' Extract sample metadata from filenames
#'
#' GSE49630 files follow this naming pattern:
#' GSM[ID]_[BodyFluid]_[Replicate]_[Platform].CEL.gz
#'
#' Example: GSM1204456_Blood_1_HuGene-1_0-st-v1.CEL.gz
#'
#' This function extracts:
#' - Sample names (without .CEL.gz extension)
#' - Body fluid types from filename patterns
#' - Replicate numbers where available

cat("\n=== Sample Metadata Extraction ===\n")

# Extract base sample names
sample_names <- gsub("\\.CEL\\.gz$", "", basename(cel_files))

# Initialize fluid type vector
fluid_types <- character(length(sample_names))

# Pattern matching for body fluid identification
# These patterns are specific to GSE49630 naming convention
fluid_patterns <- list(
  blood = c("Blood", "blood"),           # Blood samples
  semen = c("Semen", "semen"),           # Semen samples
  vaginal = c("Vaginal", "vaginal"),     # Vaginal secretion samples
  saliva = c("Saliva", "saliva")         # Saliva samples
)

# Assign fluid types based on filename patterns
for (fluid in names(fluid_patterns)) {
  for (pattern in fluid_patterns[[fluid]]) {
    matches <- grep(pattern, sample_names)
    fluid_types[matches] <- fluid
  }
}

# Check for unassigned samples
unassigned <- which(fluid_types == "")
if (length(unassigned) > 0) {
  cat("WARNING: Could not assign fluid type to samples:\n")
  for (idx in unassigned) {
    cat(sprintf("  %s\n", sample_names[idx]))
  }
  # Assign 'unknown' to unassigned samples
  fluid_types[unassigned] <- "unknown"
}

# Display sample distribution
cat("\nSample distribution by body fluid:\n")
fluid_table <- table(fluid_types)
print(fluid_table)

# Validate expected sample structure for GSE49630
expected_fluids <- c("blood", "semen", "vaginal", "saliva")
expected_per_fluid <- 5  # GSE49630 has 5 samples per fluid

for (fluid in expected_fluids) {
  count <- sum(fluid_types == fluid)
  if (count != expected_per_fluid) {
    cat(sprintf("WARNING: Expected %d %s samples, found %d\n", 
                expected_per_fluid, fluid, count))
  }
}

#' Decompress CEL files for processing
#'
#' The oligo package requires uncompressed CEL files. We decompress them
#' to a temporary directory to avoid cluttering the original data directory.
#'
#' Note: CEL files can be large (10-50MB each), so ensure sufficient
#' temporary disk space is available.

cat("\n=== CEL File Decompression ===\n")

# Create temporary directory for decompressed files
temp_dir <- tempdir()
cat(sprintf("Using temporary directory: %s\n", temp_dir))

# Initialize vector to store temporary file paths
temp_files <- character(length(cel_files))

# Decompress each file
cat("Decompressing files:\n")
for (i in 1:length(cel_files)) {
  # Create temporary filename (remove .gz extension)
  temp_files[i] <- file.path(temp_dir, gsub("\\.gz$", "", basename(cel_files[i])))
  
  # Skip if already decompressed
  if (file.exists(temp_files[i])) {
    cat(sprintf("  %d/%d: %s (already exists)\n", i, length(cel_files), basename(temp_files[i])))
  } else {
    cat(sprintf("  %d/%d: %s\n", i, length(cel_files), basename(temp_files[i])))
    
    # Decompress using system gunzip command
    # This is more reliable than R's built-in gzip functions for large files
    result <- system2("gunzip", 
                     args = c("-c", shQuote(cel_files[i])), 
                     stdout = temp_files[i],
                     stderr = FALSE)
    
    if (result != 0) {
      stop(sprintf("Failed to decompress %s", cel_files[i]))
    }
  }
}

# Verify all files were decompressed successfully
missing_temp <- !file.exists(temp_files)
if (any(missing_temp)) {
  cat("ERROR: Failed to decompress files:\n")
  for (idx in which(missing_temp)) {
    cat(sprintf("  %s\n", basename(cel_files[idx])))
  }
  stop("Decompression failed for some files")
}

cat(sprintf("✓ Successfully decompressed %d files\n", length(temp_files)))

#' Read CEL files using oligo package
#'
#' This is the critical step that reads the binary CEL file data and creates
#' an ExpressionFeatureSet object. The oligo package automatically detects
#' the array platform and applies appropriate processing.
#'
#' Memory usage: This step requires significant RAM (1-4GB depending on
#' number of arrays). Monitor system memory if processing large batches.

cat("\n=== CEL File Reading ===\n")
cat("Reading CEL files with oligo package...\n")
cat("This may take several minutes depending on array size...\n")

# Record start time for performance monitoring
start_time <- Sys.time()

# Read all CEL files into ExpressionFeatureSet object
# oligo automatically detects platform and probe design
tryCatch({
  raw_data <- read.celfiles(temp_files)
}, error = function(e) {
  cat("ERROR during CEL file reading:\n")
  cat(sprintf("  %s\n", e$message))
  cat("\nPossible causes:\n")
  cat("  - Corrupted CEL files\n")
  cat("  - Insufficient memory\n") 
  cat("  - Incompatible array platform\n")
  cat("  - Missing platform annotation packages\n")
  stop("CEL file reading failed")
})

# Calculate reading time
read_time <- Sys.time() - start_time
cat(sprintf("✓ CEL files read successfully in %.1f seconds\n", as.numeric(read_time)))

# Display information about the raw data object
cat("\nRaw data summary:\n")
cat(sprintf("  Platform: %s\n", annotation(raw_data)))
cat(sprintf("  Samples: %d\n", ncol(raw_data)))
cat(sprintf("  Probes: %d\n", nrow(raw_data)))

# Check sample names were preserved
if (length(sampleNames(raw_data)) != length(sample_names)) {
  cat("WARNING: Sample count mismatch after reading\n")
} else {
  cat("✓ All samples loaded successfully\n")
}

#' Perform RMA normalization
#'
#' RMA (Robust Multiarray Average) is the gold standard for Affymetrix
#' array normalization. It performs:
#'
#' 1. Background Correction:
#'    - Uses PM (Perfect Match) probes only
#'    - Applies convolution background correction
#'    - Removes systematic background noise
#'
#' 2. Quantile Normalization:
#'    - Makes expression distributions identical across arrays
#'    - Removes technical variation between arrays
#'    - Preserves biological differences
#'
#' 3. Probe Set Summarization:
#'    - Combines multiple probes per gene using median polish
#'    - Provides one expression value per gene per sample
#'    - Robust to outlier probes
#'
#' The result is log2-transformed expression values ready for statistical analysis.

cat("\n=== RMA Normalization ===\n")
cat("Performing RMA normalization (background correction + quantile normalization + summarization)...\n")
cat("This step may take 5-15 minutes depending on array size...\n")

# Record start time
norm_start_time <- Sys.time()

# Perform RMA normalization
# This creates an ExpressionSet object with normalized expression values
tryCatch({
  eset <- rma(raw_data)
}, error = function(e) {
  cat("ERROR during RMA normalization:\n")
  cat(sprintf("  %s\n", e$message))
  cat("\nPossible causes:\n")
  cat("  - Insufficient memory for normalization\n")
  cat("  - Corrupted probe annotation data\n")
  cat("  - Platform-specific normalization issues\n")
  stop("RMA normalization failed")
})

# Calculate normalization time
norm_time <- Sys.time() - norm_start_time
cat(sprintf("✓ RMA normalization completed in %.1f minutes\n", as.numeric(norm_time, units="mins")))

#' Extract and prepare expression matrix
#'
#' The RMA-normalized data is stored in an ExpressionSet object. We extract
#' the expression matrix and apply meaningful sample names for downstream analysis.

cat("\n=== Expression Matrix Preparation ===\n")

# Extract expression matrix (genes x samples)
expr_matrix <- exprs(eset)

# Apply meaningful sample names
colnames(expr_matrix) <- sample_names

# Display matrix dimensions and summary statistics
cat(sprintf("Expression matrix dimensions: %d probes × %d samples\n", 
            nrow(expr_matrix), ncol(expr_matrix)))

# Calculate summary statistics
cat("\nExpression value summary:\n")
cat(sprintf("  Range: %.2f to %.2f (log2 scale)\n", 
            min(expr_matrix, na.rm=TRUE), max(expr_matrix, na.rm=TRUE)))
cat(sprintf("  Median: %.2f\n", median(expr_matrix, na.rm=TRUE)))
cat(sprintf("  Missing values: %d (%.2f%%)\n", 
            sum(is.na(expr_matrix)), 
            100 * sum(is.na(expr_matrix)) / length(expr_matrix)))

# Check for potential normalization issues
if (any(is.infinite(expr_matrix))) {
  cat("WARNING: Expression matrix contains infinite values\n")
}

if (sd(colMeans(expr_matrix, na.rm=TRUE)) > 0.1) {
  cat("WARNING: High variation in sample means after normalization\n")
  cat("This may indicate normalization problems\n")
}

#' Save processed data
#'
#' Save both the expression matrix and sample metadata in CSV format
#' for easy import into downstream analysis tools (Python, R, etc.).

cat("\n=== Saving Results ===\n")

# Save expression matrix
expr_matrix_path <- file.path(output_dir, "cel_expression_matrix.csv")
write.csv(expr_matrix, expr_matrix_path, row.names = TRUE)
cat(sprintf("✓ Expression matrix saved: %s\n", expr_matrix_path))

# Create comprehensive metadata DataFrame
metadata <- data.frame(
  sample = sample_names,
  fluid_type = fluid_types,
  cel_file = basename(cel_files),
  processing_date = Sys.Date(),
  normalization_method = "RMA",
  platform = annotation(raw_data),
  stringsAsFactors = FALSE
)

# Add quality metrics if available
if (exists("raw_data")) {
  # Calculate basic quality metrics per sample
  # Note: More sophisticated QC can be added here
  metadata$mean_expression <- colMeans(expr_matrix, na.rm = TRUE)
  metadata$expression_sd <- apply(expr_matrix, 2, sd, na.rm = TRUE)
  metadata$probe_count <- rep(nrow(expr_matrix), nrow(metadata))
}

# Save metadata
metadata_path <- file.path(output_dir, "cel_metadata.csv")
write.csv(metadata, metadata_path, row.names = FALSE)
cat(sprintf("✓ Sample metadata saved: %s\n", metadata_path))

#' Quality control summary
#'
#' Generate a brief quality control report to assess the success of processing.

cat("\n=== Quality Control Summary ===\n")

# Sample-wise QC
cat("Per-sample expression summary:\n")
sample_stats <- data.frame(
  Sample = sample_names,
  Fluid = fluid_types,
  Mean = round(colMeans(expr_matrix, na.rm = TRUE), 2),
  SD = round(apply(expr_matrix, 2, sd, na.rm = TRUE), 2),
  Missing = colSums(is.na(expr_matrix))
)

# Display first few samples
print(head(sample_stats))
if (nrow(sample_stats) > 6) {
  cat(sprintf("... and %d more samples\n", nrow(sample_stats) - 6))
}

# Overall QC metrics
cat("\nOverall quality metrics:\n")
cat(sprintf("  ✓ Processing completed successfully\n"))
cat(sprintf("  ✓ %d samples processed\n", ncol(expr_matrix)))
cat(sprintf("  ✓ %d probes per sample\n", nrow(expr_matrix)))
cat(sprintf("  ✓ Mean sample correlation: %.3f\n", 
            mean(cor(expr_matrix, use="pairwise.complete.obs")[upper.tri(cor(expr_matrix, use="pairwise.complete.obs"))])))

# Check for potential batch effects by fluid type
fluid_means <- tapply(colMeans(expr_matrix, na.rm=TRUE), fluid_types, mean)
cat("\nMean expression by fluid type:\n")
for (fluid in names(fluid_means)) {
  cat(sprintf("  %s: %.2f\n", fluid, fluid_means[fluid]))
}

#' Cleanup temporary files
#'
#' Remove decompressed CEL files from temporary directory to free disk space.

cat("\n=== Cleanup ===\n")

# Remove temporary CEL files
files_removed <- 0
for (temp_file in temp_files) {
  if (file.exists(temp_file)) {
    unlink(temp_file)
    files_removed <- files_removed + 1
  }
}

if (files_removed > 0) {
  cat(sprintf("✓ Removed %d temporary files\n", files_removed))
}

#' Processing completion summary

cat("\n=== Processing Complete ===\n")
total_time <- Sys.time() - start_time
cat(sprintf("Total processing time: %.1f minutes\n", as.numeric(total_time, units="mins")))
cat(sprintf("Output directory: %s\n", normalizePath(output_dir)))
cat("\nFiles created:\n")
cat(sprintf("  - %s (expression matrix)\n", basename(expr_matrix_path)))
cat(sprintf("  - %s (sample metadata)\n", basename(metadata_path)))

# Final validation
if (file.exists(expr_matrix_path) && file.exists(metadata_path)) {
  cat("\n✓ All output files created successfully\n")
  cat("Ready for downstream analysis!\n")
} else {
  cat("\n✗ Some output files are missing\n")
  cat("Please check for errors above\n")
}

cat("\n=== Next Steps ===\n")
cat("The processed data can now be used for:\n")
cat("  1. Differential expression analysis\n")
cat("  2. Forensic marker identification\n") 
cat("  3. Statistical modeling\n")
cat("  4. Cross-platform validation\n")
cat("\nLoad the data in your analysis environment:\n")
cat("# In R:\n")
cat("expr_data <- read.csv('data/processed/cel/cel_expression_matrix.csv', row.names=1)\n")
cat("metadata <- read.csv('data/processed/cel/cel_metadata.csv')\n")
cat("\n# In Python:\n")
cat("import pandas as pd\n")
cat("expr_data = pd.read_csv('data/processed/cel/cel_expression_matrix.csv', index_col=0)\n")
cat("metadata = pd.read_csv('data/processed/cel/cel_metadata.csv')\n")