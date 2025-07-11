#!/usr/bin/env python3
"""
Final CEL Processor - Creates and executes R script for CEL processing
"""

import subprocess
import pandas as pd
from pathlib import Path
import logging
import sys

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def create_r_script(cel_dir, output_dir):
    """Create R script for CEL processing"""
    r_script = f"""
# CEL Processing Script for Forensic miRNA Analysis
# Auto-generated by cel_processor_final.py

# Set up
suppressPackageStartupMessages({{
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos="https://cloud.r-project.org")
    
    packages_needed <- c("oligo", "Biobase")
    packages_installed <- installed.packages()[, "Package"]
    missing <- packages_needed[!packages_needed %in% packages_installed]
    
    if (length(missing) > 0) {{
        cat("Installing missing packages:", missing, "\\n")
        BiocManager::install(missing, ask = FALSE, update = FALSE)
    }}
    
    library(oligo)
    library(Biobase)
}})

# Set paths
cel_dir <- "{cel_dir}"
output_dir <- "{output_dir}"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# List CEL files
cel_files <- list.files(cel_dir, pattern = "\\\\.CEL\\\\.gz$", full.names = TRUE)
cat("Found", length(cel_files), "CEL files\\n")

# Decompress files
temp_dir <- tempdir()
temp_files <- character(length(cel_files))

for (i in seq_along(cel_files)) {{
    temp_file <- file.path(temp_dir, gsub("\\\\.gz$", "", basename(cel_files[i])))
    system2("gunzip", args = c("-c", cel_files[i]), stdout = temp_file)
    temp_files[i] <- temp_file
}}

# Read CEL files
cat("Reading CEL files...\\n")
raw_data <- read.celfiles(temp_files, verbose = FALSE)

# RMA normalization
cat("Performing RMA normalization...\\n")
eset <- rma(raw_data, background = TRUE, normalize = TRUE)

# Extract expression matrix
expr_matrix <- exprs(eset)

# Clean sample names
sample_names <- gsub("\\\\.CEL$", "", basename(temp_files))
colnames(expr_matrix) <- sample_names

# Save expression matrix
output_file <- file.path(output_dir, "cel_raw_expression.csv")
write.csv(expr_matrix, output_file)
cat("Saved expression matrix to", output_file, "\\n")

# Print summary
cat("\\nProcessing complete!\\n")
cat("Dimensions:", nrow(expr_matrix), "probes x", ncol(expr_matrix), "samples\\n")
cat("Expression range:", range(expr_matrix), "\\n")

# Clean up temp files
unlink(temp_files)
"""
    
    # Save R script
    r_script_path = Path(output_dir) / "cel_processing.R"
    with open(r_script_path, 'w') as f:
        f.write(r_script)
    
    logger.info(f"Created R script: {r_script_path}")
    return r_script_path

def run_r_script(r_script_path):
    """Execute R script"""
    logger.info("Executing R script...")
    
    try:
        # Run R script
        result = subprocess.run(
            ["Rscript", str(r_script_path)],
            capture_output=True,
            text=True,
            check=True
        )
        
        logger.info("R script output:")
        print(result.stdout)
        
        if result.stderr:
            logger.warning("R script warnings:")
            print(result.stderr)
        
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"R script failed with error code {e.returncode}")
        logger.error("Error output:")
        print(e.stderr)
        return False
    except FileNotFoundError:
        logger.error("Rscript not found. Please ensure R is installed and in PATH")
        return False

def process_cel_output(output_dir):
    """Process the R output to create final miRNA expression matrix"""
    output_path = Path(output_dir)
    
    # Check if R output exists
    raw_expr_file = output_path / "cel_raw_expression.csv"
    if not raw_expr_file.exists():
        logger.error(f"Raw expression file not found: {raw_expr_file}")
        return False
    
    # Load raw expression
    logger.info("Loading raw expression data...")
    raw_expr = pd.read_csv(raw_expr_file, index_col=0)
    logger.info(f"Loaded expression matrix: {raw_expr.shape}")
    
    # Load probe mapping
    probe_mapping_file = output_path / "probe_to_mirna_mapping.csv"
    if not probe_mapping_file.exists():
        logger.error("Probe mapping file not found. Run cel_processor_simple.py first")
        return False
    
    probe_mapping = pd.read_csv(probe_mapping_file)
    probe_to_mirna = dict(zip(probe_mapping['probe_id'], probe_mapping['mirna_name']))
    
    # Map probes to miRNAs
    logger.info("Mapping probes to miRNAs...")
    mirna_data = []
    mapped_count = 0
    
    for probe_id in raw_expr.index:
        if probe_id in probe_to_mirna:
            mirna_name = probe_to_mirna[probe_id]
            mirna_data.append({
                'miRNA': mirna_name,
                **{sample: raw_expr.loc[probe_id, sample] for sample in raw_expr.columns}
            })
            mapped_count += 1
    
    logger.info(f"Mapped {mapped_count} out of {len(raw_expr)} probes to miRNAs")
    
    # Create miRNA expression matrix
    mirna_expr = pd.DataFrame(mirna_data).set_index('miRNA').T
    mirna_expr.to_csv(output_path / "cel_expression_matrix.csv")
    logger.info(f"Saved miRNA expression matrix: {mirna_expr.shape}")
    
    # Load metadata
    metadata = pd.read_csv(output_path / "cel_metadata.csv")
    
    # Create long format
    long_data = []
    for sample in mirna_expr.index:
        sample_meta = metadata[metadata['sample'] == sample].iloc[0]
        for mirna in mirna_expr.columns:
            long_data.append({
                'miRNA': mirna,
                'intensity': mirna_expr.loc[sample, mirna],
                'sample': sample,
                'body_fluid': sample_meta['body_fluid']
            })
    
    long_df = pd.DataFrame(long_data)
    long_df.to_csv(output_path / "cel_combined_raw.csv", index=False)
    
    # Quality report
    quality_data = []
    for sample in mirna_expr.index:
        sample_meta = metadata[metadata['sample'] == sample].iloc[0]
        quality_data.append({
            'sample': sample,
            'body_fluid': sample_meta['body_fluid'],
            'mean_intensity': mirna_expr.loc[sample].mean(),
            'median_intensity': mirna_expr.loc[sample].median(),
            'sd_intensity': mirna_expr.loc[sample].std(),
            'mirna_count': mirna_expr.shape[1]
        })
    
    quality_df = pd.DataFrame(quality_data)
    quality_df.to_csv(output_path / "cel_quality_report.csv", index=False)
    
    # Remove pending indicator
    pending_file = output_path / "cel_processing_pending.txt"
    if pending_file.exists():
        pending_file.unlink()
    
    logger.info("\n=== CEL Processing Complete ===")
    logger.info(f"Samples: {mirna_expr.shape[0]}")
    logger.info(f"miRNAs: {mirna_expr.shape[1]}")
    logger.info(f"Body fluids: {metadata['body_fluid'].value_counts().to_dict()}")
    
    return True

def main():
    """Main processing pipeline"""
    cel_dir = Path("data/raw/GSE49630_CEL")
    output_dir = Path("data/processed/cel")
    
    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if probe mapping exists
    if not (output_dir / "probe_to_mirna_mapping.csv").exists():
        logger.info("Running cel_processor_simple.py to create probe mapping...")
        subprocess.run([sys.executable, "scripts/preprocessing/cel_processor_simple.py"])
    
    # Create R script
    r_script_path = create_r_script(cel_dir, output_dir)
    
    # Run R script
    if run_r_script(r_script_path):
        # Process output
        if process_cel_output(output_dir):
            logger.info("CEL data successfully processed and ready for integration!")
        else:
            logger.error("Failed to process CEL output")
            sys.exit(1)
    else:
        logger.error("R script execution failed")
        logger.info("Please check that R and Bioconductor packages are installed")
        logger.info("You can manually run: Rscript " + str(r_script_path))
        sys.exit(1)

if __name__ == "__main__":
    main()