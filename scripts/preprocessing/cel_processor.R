#!/usr/bin/env Rscript

# CEL File Processor for Forensic miRNA Analysis
# Processes Affymetrix miRNA array data

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Check and install required packages
packages_needed <- c("affy", "oligo", "limma")
packages_to_install <- packages_needed[!packages_needed %in% installed.packages()[, "Package"]]

if (length(packages_to_install) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(packages_to_install)
}

# Load libraries
library(affy)
library(oligo)
library(limma)

# Set working directory to project root
# Assuming script is run from project root
setwd(".")

# Function to process CEL files
process_cel_files <- function(cel_dir, output_dir) {
  
  cat("=== CEL File Processing Pipeline ===\n")
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # List CEL files
  cel_files <- list.files(cel_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE)
  cat(sprintf("Found %d CEL files\n", length(cel_files)))
  
  # Extract sample information from filenames
  sample_info <- data.frame(
    file = cel_files,
    sample = gsub(".CEL.gz", "", basename(cel_files)),
    stringsAsFactors = FALSE
  )
  
  # Extract body fluid type from sample names
  sample_info$body_fluid <- sapply(sample_info$sample, function(x) {
    if (grepl("Blood_", x) && !grepl("Menstrual", x)) {
      return("blood")
    } else if (grepl("Semen_", x)) {
      return("semen")
    } else if (grepl("Vaginal_", x)) {
      return("vaginal")
    } else if (grepl("Saliva_", x)) {
      return("saliva")
    } else {
      return("unknown")
    }
  })
  
  cat("\nSample distribution:\n")
  print(table(sample_info$body_fluid))
  
  # Read CEL files
  cat("\nReading CEL files...\n")
  
  # Try using oligo package (preferred for newer arrays)
  tryCatch({
    # Decompress files temporarily
    temp_dir <- tempdir()
    temp_files <- c()
    
    for (i in 1:length(cel_files)) {
      temp_file <- file.path(temp_dir, gsub(".gz", "", basename(cel_files[i])))
      system2("gunzip", args = c("-c", cel_files[i]), stdout = temp_file)
      temp_files <- c(temp_files, temp_file)
    }
    
    # Read with oligo
    raw_data <- read.celfiles(temp_files)
    
    # Clean up temp files
    unlink(temp_files)
    
    cat("Successfully read CEL files with oligo\n")
    
  }, error = function(e) {
    cat("Error with oligo, trying affy package\n")
    stop("CEL reading failed. May need platform-specific CDF.")
  })
  
  # Get expression values - RMA normalization
  cat("\nPerforming RMA normalization...\n")
  eset <- rma(raw_data)
  
  # Extract expression matrix
  expr_matrix <- exprs(eset)
  
  # Update column names to match our sample names
  colnames(expr_matrix) <- sample_info$sample
  
  cat(sprintf("Expression matrix dimensions: %d probes x %d samples\n", 
              nrow(expr_matrix), ncol(expr_matrix)))
  
  # Load annotation to map probes to miRNAs
  annotation_file <- file.path(cel_dir, "GPL16384_miRNA-3_1-st-v1.annotations.20140513.csv.gz")
  
  if (file.exists(annotation_file)) {
    cat("\nLoading probe annotations...\n")
    
    # Read annotation file
    annot <- read.csv(gzfile(annotation_file), 
                      skip = 4,  # Skip header lines
                      stringsAsFactors = FALSE)
    
    # Match probes to miRNA names
    probe_ids <- rownames(expr_matrix)
    
    # Create mapping
    probe_to_mirna <- annot[match(probe_ids, annot$Probe.Set.ID), 
                            c("Probe.Set.ID", "Probe.Set.Name")]
    
    # Filter for miRNA probes
    mirna_mask <- grepl("miR", probe_to_mirna$Probe.Set.Name)
    mirna_expr <- expr_matrix[mirna_mask, ]
    mirna_names <- probe_to_mirna$Probe.Set.Name[mirna_mask]
    
    # Clean miRNA names (remove _st suffix)
    mirna_names <- gsub("_st$", "", mirna_names)
    rownames(mirna_expr) <- mirna_names
    
    cat(sprintf("Filtered to %d miRNA probes\n", nrow(mirna_expr)))
    
  } else {
    cat("Warning: Annotation file not found. Using probe IDs.\n")
    mirna_expr <- expr_matrix
  }
  
  # Save expression matrix
  output_file <- file.path(output_dir, "cel_expression_matrix.csv")
  write.csv(mirna_expr, output_file)
  cat(sprintf("\nSaved expression matrix to %s\n", output_file))
  
  # Save sample metadata
  metadata_file <- file.path(output_dir, "cel_sample_metadata.csv")
  write.csv(sample_info[, c("sample", "body_fluid")], 
            metadata_file, row.names = FALSE)
  cat(sprintf("Saved metadata to %s\n", metadata_file))
  
  # Create long format for compatibility with Python pipeline
  long_data <- data.frame()
  for (i in 1:ncol(mirna_expr)) {
    sample_data <- data.frame(
      miRNA = rownames(mirna_expr),
      intensity = mirna_expr[, i],
      sample = colnames(mirna_expr)[i],
      body_fluid = sample_info$body_fluid[i],
      stringsAsFactors = FALSE
    )
    long_data <- rbind(long_data, sample_data)
  }
  
  long_file <- file.path(output_dir, "cel_combined_raw.csv")
  write.csv(long_data, long_file, row.names = FALSE)
  cat(sprintf("Saved long format to %s\n", long_file))
  
  # Quality metrics
  quality_report <- data.frame(
    sample = colnames(mirna_expr),
    body_fluid = sample_info$body_fluid,
    mean_intensity = colMeans(mirna_expr),
    median_intensity = apply(mirna_expr, 2, median),
    sd_intensity = apply(mirna_expr, 2, sd)
  )
  
  quality_file <- file.path(output_dir, "cel_quality_report.csv")
  write.csv(quality_report, quality_file, row.names = FALSE)
  cat(sprintf("Saved quality report to %s\n", quality_file))
  
  # Summary statistics
  cat("\n=== Processing Summary ===\n")
  cat(sprintf("Total samples: %d\n", ncol(mirna_expr)))
  cat(sprintf("Total miRNA probes: %d\n", nrow(mirna_expr)))
  cat("\nIntensity range:\n")
  print(summary(as.vector(mirna_expr)))
  
  return(list(
    expression = mirna_expr,
    metadata = sample_info,
    quality = quality_report
  ))
}

# Main execution
if (!interactive()) {
  cel_dir <- "data/raw/GSE49630_CEL"
  output_dir <- "data/processed/cel"
  
  # Process files
  results <- process_cel_files(cel_dir, output_dir)
  
  cat("\nCEL processing complete!\n")
}