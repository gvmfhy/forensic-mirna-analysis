#!/usr/bin/env Rscript

# Platform Comparison Script
# Assesses compatibility between CEL and GPR platforms

library(tidyverse)

# Function to extract GPR data
parse_gpr_file <- function(gpr_path) {
  # Read the file
  lines <- readLines(gpr_path)
  
  # Find where the data starts (after headers)
  data_start <- which(grepl("^\"Block\"", lines))
  
  # Read the data portion
  gpr_data <- read.table(
    text = lines[(data_start):length(lines)],
    header = TRUE,
    sep = "\t",
    quote = "\"",
    stringsAsFactors = FALSE
  )
  
  return(gpr_data)
}

# Function to summarize platform characteristics
summarize_platform <- function() {
  cat("=== Platform Comparison Summary ===\n\n")
  
  # CEL Platform Info
  cat("GSE49630 (Affymetrix miRNA-3_0):\n")
  cat("- Binary CEL format\n")
  cat("- Single-channel intensity data\n")
  cat("- Requires affy/oligo package for reading\n")
  cat("- Probe-level data with PM/MM probes\n\n")
  
  # GPR Platform Info
  cat("GSE153135 (miRCURY LNA):\n")
  cat("- Text-based GPR format\n")
  cat("- Two-channel (Cy3/Cy5) intensities\n")
  cat("- Direct spot intensities\n")
  cat("- Background subtraction included\n\n")
  
  # Key differences
  cat("Key Integration Challenges:\n")
  cat("1. Single vs dual-channel data\n")
  cat("2. Probe sets vs direct spots\n")
  cat("3. Different miRNA coverage\n")
  cat("4. Platform-specific biases\n\n")
  
  # Integration strategy
  cat("Integration Strategy:\n")
  cat("1. Extract single intensity measure from each platform\n")
  cat("2. Map to common miRNA identifiers\n")
  cat("3. Quantile normalize separately\n")
  cat("4. Apply ComBat for batch correction\n")
}

# Test GPR parsing
test_gpr_parsing <- function() {
  gpr_files <- list.files(
    "data/raw/GSE153135_GPR/",
    pattern = "\\.gpr\\.gz$",
    full.names = TRUE
  )
  
  if (length(gpr_files) > 0) {
    # Decompress and parse first file
    temp_file <- tempfile()
    system2("gunzip", args = c("-c", gpr_files[1]), stdout = temp_file)
    
    gpr_data <- parse_gpr_file(temp_file)
    unlink(temp_file)
    
    cat("\nGPR File Structure:\n")
    cat("Dimensions:", nrow(gpr_data), "spots x", ncol(gpr_data), "columns\n")
    cat("\nKey columns:\n")
    print(names(gpr_data)[c(1:7, 9:10, 21:22)])
    
    # Check for miRNA names
    mirna_spots <- gpr_data %>%
      filter(grepl("miR", Name))
    
    cat("\nmiRNA spots found:", nrow(mirna_spots), "\n")
    cat("Example miRNAs:\n")
    print(head(unique(mirna_spots$Name), 10))
  }
}

# Main execution
if (!interactive()) {
  summarize_platform()
  
  # Create utils directory if needed
  dir.create("scripts/utils", recursive = TRUE, showWarnings = FALSE)
  
  # Test GPR parsing if files exist
  if (file.exists("data/raw/GSE153135_GPR/")) {
    test_gpr_parsing()
  }
}