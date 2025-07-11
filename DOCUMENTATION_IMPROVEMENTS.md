# Documentation Improvements Guide

This guide provides specific improvements for code documentation.

## 1. Add Module-Level Constants


### gpr_continuous_analysis.py
```python
# Add at top of file after imports:

# Analysis thresholds and constants
ZERO_OFFSET = 0.01  # Small constant to avoid log(0) in calculations
DONOR_MULTIPLIER = 2.5  # Average donors per pooled sample
HIGH_CONFIDENCE_THRESHOLD = -2.0  # Log2 ratio threshold for high confidence detection
MODERATE_CONFIDENCE_THRESHOLD = -6.0  # Log2 ratio threshold for moderate confidence
LOW_CONFIDENCE_THRESHOLD = -10.0  # Log2 ratio threshold for low confidence
```


### cel_practical_forensic.py
```python
# Add at top of file after imports:

# Analysis thresholds and constants
P_VALUE_THRESHOLD = 0.01  # Nominal p-value threshold for significance
FOLD_CHANGE_THRESHOLD = 2  # Log2 fold change threshold (4x linear)
EFFECT_SIZE_THRESHOLD = 1.5  # Cohen's d threshold for very large effect
DETECTION_RATE_THRESHOLD = 1.0  # Minimum detection rate in target fluid
```


## 2. Improve Function Docstrings


### cel_practical_forensic.py

#### Function: `practical_forensic_analysis`
```python
def practical_forensic_analysis(...):
    """
    Perform practical forensic analysis on CEL data with relaxed statistical criteria

    Identifies miRNA markers suitable for forensic body fluid identification

    Args:
    results_dir (Path, optional): Directory containing preprocessed results
        Default: Path('results/cel_continuous_expression')

    Returns:
    dict: Dictionary mapping body fluids to their marker DataFrames
        Keys: 'blood', 'saliva', 'semen', 'vaginal'
        Values: DataFrames with columns for miRNA, statistics, and rankings

    Note:
        With small sample sizes (n=5), strict FDR correction is too conservative.
This function applies practical criteria focusing on effect size and
biological relevance rather than strict statistical significance.
    """
```

#### Function: `create_practical_visualizations`
```python
def create_practical_visualizations(...):
    """
    Create comprehensive visualizations for forensic marker analysis

    Args:
    markers_by_fluid (dict): Dictionary of DataFrames with markers per fluid
    output_dir (Path): Directory to save visualization files

    Returns:
    None: Saves PNG files to output directory

    """
```

### gpr_continuous_analysis.py

#### Function: `calculate_expression_tiers`
```python
def calculate_expression_tiers(...):
    """
    Categorize miRNA expression into confidence tiers based on intensity

    Returns:
    pd.DataFrame: Tier assignments with columns:
        - mirna: miRNA identifier
        - fluid: body fluid type
        - expression_tier: assigned confidence level
        - log2_intensity: original log2 ratio value

    """
```

#### Function: `wilcoxon_tests`
```python
def wilcoxon_tests(...):
    """
    Perform Wilcoxon rank-sum tests comparing each fluid vs all others

    Args:
    multiple_correction (str): Correction method for multiple testing
        Options: 'bonferroni' (conservative), 'fdr' (recommended)
        Default: 'bonferroni' 

    """
```

## 3. Add Inline Comments


### gpr_continuous_analysis.py

line_156:
```python
# Add small constant to prevent division by zero while maintaining ratio relationships
specificity = (min_in_fluid + 0.01) / (max_in_others + 0.01)
```

line_337:
```python
# Estimate total donors: ~2.5 donors pooled per array based on experimental design
donor_count = int(count * 2.5)
```

line_65:
```python
# Combined score weights significance, effect size, and specificity equally
score = -np.log10(row["padj"] + 1e-10) * abs(row["log2fc"]) * spec_score
```

### cel_practical_forensic.py

line_65:
```python
# Scoring formula: significance × fold change × effect size
        # Higher scores indicate better forensic markers
candidates['score'] = (
```

line_115:
```python
# Mathematical limit: With n=5, minimum p-value ≈ 0.001
    # Multiple testing burden: 0.05 / 13,564 = 0.0000037
report.append("- Testing 13,564 miRNAs requires p < 0.0000037 for FDR < 0.05")
```

## 4. Add Type Hints


### All Python files
```python
from typing import Dict, List, Tuple, Optional, Union
import pandas as pd
from pathlib import Path

# Example function signatures with type hints:
def load_data(self) -> None:
def calculate_expression_tiers(self) -> pd.DataFrame:
def wilcoxon_tests(self, multiple_correction: str = 'bonferroni') -> pd.DataFrame:
def create_visualizations(self, output_dir: Path) -> None:
```

## 5. Add R Documentation


### cel_processor_minimal.R
```r
#' Process Affymetrix miRNA ST Array CEL Files
#'
#' @description
#' Reads CEL files from GSE49630 and performs RMA normalization
#' using the oligo package. The oligo package is required for
#' ST arrays as the older affy package does not support them.
#'
#' @param cel_dir Character. Directory containing .CEL.gz files
#' @param output_dir Character. Directory for output files
#'
#' @details
#' This function processes Affymetrix Human Gene 1.0 ST arrays
#' which require specialized handling compared to 3' IVT arrays.
#'
#' @return
#' Invisible NULL. Writes two files:
#' \itemize{
#'   \item cel_expression_matrix.csv: Normalized expression values
#'   \item cel_metadata.csv: Sample annotations
#' }
#'
#' @examples
#' \dontrun{
#' process_cel_files('data/GSE49630_RAW', 'data/processed/cel')
#' }
#'
#' @export
process_cel_files <- function(cel_dir, output_dir) {
```

## 6. Files to Remove or Archive


Move these abandoned/test files to an `archive/` directory:
- `scripts/preprocessing/cel_processor.R` (replaced by minimal version)
- `scripts/preprocessing/cel_processor_simple.py` (incomplete)
- `scripts/preprocessing/cel_processor_rpy2.py` (abandoned approach)
- `scripts/preprocessing/cel_processor_final.py` (another attempt)
- `minimal_cel_process.R` (root directory duplicate)
- `process_affymetrix_mirna.R` (alternative approach)

## 7. Create Preprocessing README


### scripts/preprocessing/README.md
```markdown
# Preprocessing Scripts

## Overview
This directory contains scripts for processing raw microarray data from two platforms:
- GenePix Result (GPR) files from two-color Agilent arrays
- CEL files from Affymetrix miRNA ST arrays

## Production Scripts

### gpr_parser.py
Processes GPR files from GSE153135 (pooled body fluid samples)
- Input: .gpr or .gpr.gz files
- Output: Normalized log2 ratios (Cy5/Cy3)
- Key feature: Background correction

### cel_processor_minimal.R
Processes CEL files from GSE49630 (individual samples)
- Input: .CEL.gz files
- Output: RMA-normalized expression matrix
- Requirement: oligo package (NOT affy)

## Requirements

### Python
- Python 3.8+
- pandas, numpy
- See requirements.txt

### R
- R 4.5+
- BiocManager::install('oligo')
- BiocManager::install('limma')

## Usage Examples

```bash
# Process GPR files
python gpr_parser.py \
    --input data/GSE153135_RAW \
    --output data/processed/gpr

# Process CEL files
Rscript cel_processor_minimal.R \
    data/GSE49630_RAW \
    data/processed/cel
```
```