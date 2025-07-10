#!/usr/bin/env python3
"""
Alternative CEL processor using Python
Since R setup can be complex, we'll parse the annotation file
and prepare for R processing later
"""

import pandas as pd
import numpy as np
import gzip
import os
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_cel_annotation(annotation_file):
    """Parse Affymetrix annotation file to get probe-to-miRNA mapping"""
    logger.info("Parsing CEL annotation file")
    
    probe_mapping = {}
    mirna_probes = []
    
    with gzip.open(annotation_file, 'rt', encoding='latin-1') as f:
        # Skip header lines
        for line in f:
            if line.startswith('"Probe Set ID"'):
                break
        
        # Parse data
        for line in f:
            fields = line.strip().split('","')
            if len(fields) >= 2:
                probe_id = fields[0].strip('"')
                probe_name = fields[1].strip('"')
                
                if 'miR' in probe_name:
                    # Clean probe name
                    clean_name = probe_name.replace('_st', '').replace('_x_st', '')
                    probe_mapping[probe_id] = clean_name
                    mirna_probes.append(probe_id)
    
    logger.info(f"Found {len(mirna_probes)} miRNA probes in annotation")
    return probe_mapping, mirna_probes

def prepare_cel_metadata(cel_dir):
    """Prepare metadata for CEL files"""
    cel_files = list(Path(cel_dir).glob("*.CEL.gz"))
    logger.info(f"Found {len(cel_files)} CEL files")
    
    metadata = []
    for cel_file in cel_files:
        sample_name = cel_file.stem.replace('.CEL', '')
        
        # Extract body fluid type
        if 'Blood_' in sample_name and 'Menstrual' not in sample_name:
            body_fluid = 'blood'
        elif 'Semen_' in sample_name:
            body_fluid = 'semen'
        elif 'Vaginal_' in sample_name:
            body_fluid = 'vaginal_secretion'
        elif 'Saliva_' in sample_name:
            body_fluid = 'saliva'
        else:
            body_fluid = 'unknown'
        
        metadata.append({
            'file': str(cel_file),
            'sample': sample_name,
            'body_fluid': body_fluid
        })
    
    metadata_df = pd.DataFrame(metadata)
    
    # Summary
    logger.info("\nSample distribution:")
    print(metadata_df.groupby('body_fluid').size())
    
    return metadata_df

def create_r_processing_script(metadata_df, output_path):
    """Create a simplified R script for CEL processing"""
    r_script = """
# Simplified CEL processing script
# Run this in R after installing: BiocManager::install(c("affy", "oligo"))

library(oligo)

# File paths
cel_files <- c({files})
sample_names <- c({samples})

# Read CEL files
cat("Reading CEL files...\\n")
raw_data <- read.celfiles(cel_files)

# RMA normalization
cat("Performing RMA normalization...\\n")
eset <- rma(raw_data)

# Extract expression matrix
expr_matrix <- exprs(eset)
colnames(expr_matrix) <- sample_names

# Save as CSV
write.csv(expr_matrix, "{output_file}")
cat("Saved expression matrix to {output_file}\\n")

# Print summary
cat(sprintf("Processed %d probes x %d samples\\n", nrow(expr_matrix), ncol(expr_matrix)))
""".format(
        files=',\n  '.join([f'"{f}"' for f in metadata_df['file'].values]),
        samples=',\n  '.join([f'"{s}"' for s in metadata_df['sample'].values]),
        output_file=output_path
    )
    
    return r_script

def main():
    """Main processing function"""
    cel_dir = Path("data/raw/GSE49630_CEL")
    output_dir = Path("data/processed/cel")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Parse annotation
    annotation_file = cel_dir / "GPL16384_miRNA-3_1-st-v1.annotations.20140513.csv.gz"
    if annotation_file.exists():
        probe_mapping, mirna_probes = parse_cel_annotation(annotation_file)
        
        # Save probe mapping
        mapping_df = pd.DataFrame(list(probe_mapping.items()), 
                                columns=['probe_id', 'mirna_name'])
        mapping_df.to_csv(output_dir / 'probe_to_mirna_mapping.csv', index=False)
        logger.info(f"Saved probe mapping to {output_dir / 'probe_to_mirna_mapping.csv'}")
    
    # Prepare metadata
    metadata_df = prepare_cel_metadata(cel_dir)
    metadata_df.to_csv(output_dir / 'cel_metadata.csv', index=False)
    logger.info(f"Saved metadata to {output_dir / 'cel_metadata.csv'}")
    
    # Create R script for later processing
    r_script = create_r_processing_script(
        metadata_df, 
        str(output_dir / 'cel_raw_expression.csv')
    )
    
    with open(output_dir / 'process_cel_files.R', 'w') as f:
        f.write(r_script)
    logger.info(f"Created R script at {output_dir / 'process_cel_files.R'}")
    
    # Create sample preprocessing info
    logger.info("\n=== CEL Preprocessing Summary ===")
    logger.info(f"Total CEL files: {len(metadata_df)}")
    logger.info(f"Total miRNA probes: {len(probe_mapping)}")
    logger.info("\nNext steps:")
    logger.info("1. Run the generated R script to process CEL files")
    logger.info("2. Use probe_to_mirna_mapping.csv to map probes to miRNAs")
    logger.info("3. Merge with GPR data using common miRNA names")
    
    # Create instructions for next steps
    logger.info("\nNext steps for CEL processing:")
    logger.info("1. Install R packages: BiocManager::install(c('affy', 'oligo'))")
    logger.info("2. Run the generated R script: Rscript data/processed/cel/process_cel_files.R")
    logger.info("3. This will create cel_raw_expression.csv")
    logger.info("4. Run cel_postprocessor.py to complete the preprocessing")
    
    # Create a file to indicate CEL processing is pending
    with open(output_dir / 'cel_processing_pending.txt', 'w') as f:
        f.write("CEL files are ready for R processing.\n")
        f.write("Run: Rscript data/processed/cel/process_cel_files.R\n")
    
    logger.info("CEL preprocessing preparation complete!")

if __name__ == "__main__":
    main()