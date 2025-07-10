#!/usr/bin/env python3
"""
Examine raw data structure without processing
Understanding what we have before applying methods
"""

import gzip
from pathlib import Path
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def examine_gpr_structure():
    """Look at GPR file structure in detail"""
    logger.info("\n=== GPR Data Structure ===")
    
    # Read one GPR file
    gpr_file = Path('data/raw/GSE153135_GPR/GSM4634783_saliva-array1.gpr.gz')
    
    with gzip.open(gpr_file, 'rt') as f:
        lines = f.readlines()
    
    # Find header end
    header_end = 0
    for i, line in enumerate(lines):
        if line.startswith('"Block"'):
            header_end = i
            break
    
    logger.info(f"Header lines: {header_end}")
    logger.info(f"Data rows: {len(lines) - header_end - 1}")
    
    # Parse column names
    columns = lines[header_end].strip().split('\t')
    logger.info(f"Total columns: {len(columns)}")
    
    # Key columns for forensics
    key_cols = ['Name', 'F635 Median', 'B635 Median', 'F532 Median', 'B532 Median']
    logger.info("\nKey columns for analysis:")
    for col in key_cols:
        if col in columns:
            logger.info(f"  - {col}: position {columns.index(col)}")
    
    # Sample first few data rows
    logger.info("\nFirst 3 miRNA entries:")
    for i in range(header_end + 1, min(header_end + 4, len(lines))):
        parts = lines[i].strip().split('\t')
        name_idx = columns.index('"Name"')
        logger.info(f"  {parts[name_idx]}")
        

def examine_cel_structure():
    """Look at CEL dataset structure"""
    logger.info("\n=== CEL Data Structure ===")
    
    # List CEL files by fluid type
    cel_dir = Path('data/raw/GSE49630_CEL')
    cel_files = list(cel_dir.glob('*.CEL.gz'))
    
    # Group by fluid
    fluids = {}
    for f in cel_files:
        if 'Blood' in f.name and 'GSM' in f.name:
            fluids.setdefault('Blood', []).append(f.name)
        elif 'Semen' in f.name:
            fluids.setdefault('Semen', []).append(f.name)
        elif 'Vaginal' in f.name:
            fluids.setdefault('Vaginal', []).append(f.name)
        elif 'Saliva' in f.name:
            fluids.setdefault('Saliva', []).append(f.name)
    
    logger.info("Samples per fluid:")
    for fluid, files in fluids.items():
        logger.info(f"  {fluid}: {len(files)} samples")
        for f in sorted(files)[:2]:  # Show first 2
            logger.info(f"    - {f}")
            
    # Check annotation file
    annot_file = cel_dir / 'GPL16384_miRNA-3_1-st-v1.annotations.20140513.csv.gz'
    if annot_file.exists():
        logger.info(f"\nAnnotation file found: {annot_file.name}")
        # Peek at structure - skip comment lines
        try:
            df = pd.read_csv(annot_file, comment='#', nrows=5)
            logger.info(f"Annotation columns: {list(df.columns)[:5]}...")
        except:
            logger.info("Complex annotation format - needs special parsing")
        

def compare_data_characteristics():
    """Compare key characteristics between datasets"""
    logger.info("\n=== Dataset Comparison ===")
    
    comparison = """
    | Characteristic | GPR (GSE153135) | CEL (GSE49630) |
    |----------------|-----------------|----------------|
    | Platform       | GenePix 2-color | Affymetrix     |
    | Samples/fluid  | 2               | 5              |
    | Total samples  | 10              | 20             |
    | Fluids         | 5 (incl. menst) | 4 (no menst)   |
    | Data type      | Ratios (Cy5/Cy3)| Absolute       |
    | Normalization  | Loess typical   | RMA typical    |
    """
    logger.info(comparison)
    
    logger.info("\nForensic Implications:")
    logger.info("- GPR: Relative expression makes direct comparison harder")
    logger.info("- CEL: Absolute values better for presence/absence")
    logger.info("- Both: Need platform-specific quality thresholds")
    logger.info("- Integration: Rank-based methods most appropriate")
    

def suggest_analysis_priorities():
    """Based on data structure, what should we prioritize?"""
    logger.info("\n=== Analysis Priorities ===")
    
    priorities = """
    1. **Quality Control First**
       - Detection call rates per platform
       - Background vs signal distributions
       - Sample clustering by fluid type
       
    2. **Platform-Specific Analysis**
       - GPR: Focus on high-ratio miRNAs
       - CEL: Focus on consistently detected miRNAs
       - Both: Binary presence/absence may be most robust
       
    3. **Forensic Marker Criteria**
       - Specificity > Sensitivity for court admissibility
       - Biological plausibility check against literature
       - Stability/degradation resistance important
       
    4. **Statistical Approach**
       - With n=5 (CEL): Wilcoxon test appropriate
       - No parametric assumptions needed
       - Multiple testing correction essential
       
    5. **Validation Strategy**
       - If marker found in both platforms = strong candidate
       - Platform-specific markers need extra scrutiny
       - Consider technical replicates vs biological replicates
    """
    logger.info(priorities)


def main():
    """Run data structure examination"""
    examine_gpr_structure()
    examine_cel_structure()
    compare_data_characteristics()
    suggest_analysis_priorities()
    
    logger.info("\n=== Next Steps ===")
    logger.info("1. Process CEL files with R script")
    logger.info("2. Run quality assessment on processed data")
    logger.info("3. Apply forensic marker criteria")
    logger.info("4. Compare findings between platforms")


if __name__ == "__main__":
    main()