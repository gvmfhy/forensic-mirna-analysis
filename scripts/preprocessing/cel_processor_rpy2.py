#!/usr/bin/env python3
"""
CEL File Processor using rpy2
Processes Affymetrix CEL files using R's Bioconductor packages through Python
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import os
import gzip
import tempfile
import shutil

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def check_r_packages():
    """Check if required R packages are installed"""
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()
        
        # Check for required packages
        r_code = """
        packages_needed <- c("affy", "oligo", "Biobase")
        packages_installed <- installed.packages()[, "Package"]
        missing <- packages_needed[!packages_needed %in% packages_installed]
        
        if (length(missing) > 0) {
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager", repos="https://cloud.r-project.org")
            }
            BiocManager::install(missing, ask = FALSE, update = FALSE)
        }
        
        # Load libraries
        suppressPackageStartupMessages({
            library(oligo)
            library(Biobase)
        })
        
        TRUE
        """
        
        result = ro.r(r_code)
        logger.info("R packages checked and loaded successfully")
        return True
        
    except Exception as e:
        logger.error(f"Error checking R packages: {e}")
        logger.info("Trying alternative approach...")
        return False

def process_cel_files_rpy2(cel_dir, output_dir):
    """Process CEL files using rpy2"""
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()
        
        logger.info("Processing CEL files with rpy2...")
        
        # Create output directory
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # List CEL files
        cel_files = list(Path(cel_dir).glob("*.CEL.gz"))
        logger.info(f"Found {len(cel_files)} CEL files")
        
        # Decompress CEL files to temporary directory
        temp_dir = tempfile.mkdtemp()
        temp_files = []
        
        logger.info("Decompressing CEL files...")
        for cel_file in cel_files:
            temp_file = os.path.join(temp_dir, cel_file.stem)
            with gzip.open(cel_file, 'rb') as f_in:
                with open(temp_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            temp_files.append(temp_file)
        
        # Convert file paths to R vector
        ro.globalenv['cel_files'] = ro.StrVector(temp_files)
        
        # Process with R
        logger.info("Reading CEL files with oligo package...")
        r_processing = """
        # Read CEL files
        raw_data <- read.celfiles(cel_files, verbose = FALSE)
        
        # Perform RMA normalization
        eset <- rma(raw_data, background = TRUE, normalize = TRUE)
        
        # Extract expression matrix
        expr_matrix <- exprs(eset)
        
        # Get sample names
        sample_names <- sampleNames(eset)
        
        list(expr_matrix = expr_matrix, sample_names = sample_names)
        """
        
        result = ro.r(r_processing)
        
        # Convert R matrix to pandas DataFrame
        expr_matrix = np.array(result[0])
        sample_names = list(result[1])
        
        # Clean sample names
        sample_names = [os.path.basename(name).replace('.CEL', '') for name in sample_names]
        
        # Create DataFrame
        expr_df = pd.DataFrame(expr_matrix.T, columns=range(expr_matrix.shape[0]))
        expr_df.index = sample_names
        
        # Save raw expression matrix
        expr_df.to_csv(output_path / 'cel_raw_expression.csv')
        logger.info(f"Saved raw expression matrix: {expr_df.shape}")
        
        # Clean up temporary files
        shutil.rmtree(temp_dir)
        
        return expr_df
        
    except Exception as e:
        logger.error(f"Error processing CEL files with rpy2: {e}")
        return None

def map_probes_to_mirnas(expr_df, annotation_file, output_dir):
    """Map probe IDs to miRNA names"""
    logger.info("Mapping probes to miRNAs...")
    
    # Load probe mapping
    probe_mapping = pd.read_csv(Path(output_dir) / 'probe_to_mirna_mapping.csv')
    
    # Create mapping dictionary
    probe_to_mirna = dict(zip(probe_mapping['probe_id'], probe_mapping['mirna_name']))
    
    # Map column names
    mirna_names = []
    mirna_indices = []
    
    for i, probe in enumerate(expr_df.columns):
        probe_str = f"{probe}_st"  # Affymetrix naming convention
        if probe_str in probe_to_mirna:
            mirna_names.append(probe_to_mirna[probe_str])
            mirna_indices.append(i)
    
    # Filter to miRNA probes only
    mirna_expr = expr_df.iloc[:, mirna_indices]
    mirna_expr.columns = mirna_names
    
    logger.info(f"Mapped {len(mirna_names)} miRNA probes")
    
    # Save miRNA expression matrix
    output_path = Path(output_dir)
    mirna_expr.to_csv(output_path / 'cel_expression_matrix.csv')
    
    # Create long format for compatibility
    long_data = []
    metadata = pd.read_csv(output_path / 'cel_metadata.csv')
    
    for sample in mirna_expr.index:
        sample_metadata = metadata[metadata['sample'] == sample].iloc[0]
        for mirna in mirna_expr.columns:
            long_data.append({
                'miRNA': mirna,
                'intensity': mirna_expr.loc[sample, mirna],
                'sample': sample,
                'body_fluid': sample_metadata['body_fluid']
            })
    
    long_df = pd.DataFrame(long_data)
    long_df.to_csv(output_path / 'cel_combined_raw.csv', index=False)
    
    # Quality report
    quality_report = []
    for sample in mirna_expr.index:
        sample_metadata = metadata[metadata['sample'] == sample].iloc[0]
        quality_report.append({
            'sample': sample,
            'body_fluid': sample_metadata['body_fluid'],
            'mean_intensity': mirna_expr.loc[sample].mean(),
            'median_intensity': mirna_expr.loc[sample].median(),
            'sd_intensity': mirna_expr.loc[sample].std()
        })
    
    quality_df = pd.DataFrame(quality_report)
    quality_df.to_csv(output_path / 'cel_quality_report.csv', index=False)
    
    return mirna_expr

def main():
    """Main processing pipeline"""
    cel_dir = "data/raw/GSE49630_CEL"
    output_dir = "data/processed/cel"
    
    # Check if we can use rpy2
    if check_r_packages():
        # Process with rpy2
        expr_df = process_cel_files_rpy2(cel_dir, output_dir)
        
        if expr_df is not None:
            # Map probes to miRNAs
            annotation_file = Path(cel_dir) / "GPL16384_miRNA-3_1-st-v1.annotations.20140513.csv.gz"
            mirna_expr = map_probes_to_mirnas(expr_df, annotation_file, output_dir)
            
            logger.info("\n=== CEL Processing Complete ===")
            logger.info(f"Processed {mirna_expr.shape[0]} samples")
            logger.info(f"Extracted {mirna_expr.shape[1]} miRNAs")
            logger.info("Data ready for integration!")
            
            # Remove pending indicator
            pending_file = Path(output_dir) / 'cel_processing_pending.txt'
            if pending_file.exists():
                pending_file.unlink()
        else:
            logger.error("CEL processing failed")
            logger.info("Please ensure R and Bioconductor packages are installed")
            logger.info("Try: R -e 'BiocManager::install(c(\"affy\", \"oligo\"))'")
    else:
        logger.error("Could not load required R packages")
        logger.info("Creating R script for manual processing...")
        
        # Re-run cel_processor_simple to create R script
        os.system("python scripts/preprocessing/cel_processor_simple.py")

if __name__ == "__main__":
    main()