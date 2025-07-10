#!/usr/bin/env Rscript

# Platform Comparison Script (Base R version)
# Assesses compatibility between CEL and GPR platforms

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

# Main execution
cat("=== Platform Comparison Summary ===\n\n")

# CEL Platform Info
cat("GSE49630 (Affymetrix miRNA-3_0):\n")
cat("- Binary CEL format\n")
cat("- Single-channel intensity data\n")
cat("- 20 samples: 5 each of blood, saliva, semen, vaginal\n")
cat("- Missing: menstrual blood\n\n")

# Count CEL files
cel_files <- list.files("data/raw/GSE49630_CEL/", pattern = "\\.CEL\\.gz$")
cat("CEL files found:", length(cel_files), "\n\n")

# GPR Platform Info
cat("GSE153135 (miRCURY LNA):\n")
cat("- Text-based GPR format\n")
cat("- Two-channel (Cy3/Cy5) intensities\n")
cat("- 10 samples including menstrual blood\n\n")

# Count GPR files
gpr_files <- list.files("data/raw/GSE153135_GPR/", pattern = "\\.gpr\\.gz$")
cat("GPR files found:", length(gpr_files), "\n")

# Parse first GPR file to check structure
if (length(gpr_files) > 0) {
  cat("\n--- Parsing first GPR file ---\n")
  temp_file <- tempfile()
  system2("gunzip", args = c("-c", paste0("data/raw/GSE153135_GPR/", gpr_files[1])), stdout = temp_file)
  
  gpr_data <- parse_gpr_file(temp_file)
  unlink(temp_file)
  
  cat("File:", gpr_files[1], "\n")
  cat("Total spots:", nrow(gpr_data), "\n")
  
  # Count miRNA spots
  mirna_spots <- grep("miR", gpr_data$Name)
  cat("miRNA spots:", length(mirna_spots), "\n")
  
  # Show some example miRNAs
  mirna_names <- unique(gpr_data$Name[mirna_spots])
  cat("\nExample miRNAs:\n")
  print(head(mirna_names[mirna_names != ""], 10))
  
  # Check intensity ranges
  f635 <- gpr_data$F635.Median[mirna_spots]
  f532 <- gpr_data$F532.Median[mirna_spots]
  
  cat("\nIntensity ranges:\n")
  cat("Channel 635nm (Red):", range(f635, na.rm = TRUE), "\n")
  cat("Channel 532nm (Green):", range(f532, na.rm = TRUE), "\n")
}

cat("\n=== Integration Challenges ===\n")
cat("1. Single-channel (CEL) vs dual-channel (GPR) data\n")
cat("2. Different array designs and probe sets\n")
cat("3. Platform-specific background correction\n")
cat("4. Potentially different miRNA coverage\n")

cat("\n=== Recommended Integration Approach ===\n")
cat("1. For GPR: Use Cy5/Cy3 ratio or single channel\n")
cat("2. Map both platforms to common miRNA names\n")
cat("3. Quantile normalize within platform\n")
cat("4. Apply ComBat batch correction\n")
cat("5. Validate with PCA and correlation analysis\n")