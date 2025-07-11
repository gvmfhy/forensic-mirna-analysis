#!/usr/bin/env python3
"""
Script to add comprehensive documentation to all analysis scripts.
This will add proper docstrings, type hints, and inline comments.
"""

import os
from pathlib import Path

# Define constants that should be documented
CONSTANTS_TO_ADD = {
    'gpr_continuous_analysis.py': [
        ('ZERO_OFFSET', '0.01', 'Small constant to avoid log(0) in calculations'),
        ('DONOR_MULTIPLIER', '2.5', 'Average donors per pooled sample'),
        ('HIGH_CONFIDENCE_THRESHOLD', '-2.0', 'Log2 ratio threshold for high confidence detection'),
        ('MODERATE_CONFIDENCE_THRESHOLD', '-6.0', 'Log2 ratio threshold for moderate confidence'),
        ('LOW_CONFIDENCE_THRESHOLD', '-10.0', 'Log2 ratio threshold for low confidence'),
    ],
    'cel_practical_forensic.py': [
        ('P_VALUE_THRESHOLD', '0.01', 'Nominal p-value threshold for significance'),
        ('FOLD_CHANGE_THRESHOLD', '2', 'Log2 fold change threshold (4x linear)'),
        ('EFFECT_SIZE_THRESHOLD', '1.5', "Cohen's d threshold for very large effect"),
        ('DETECTION_RATE_THRESHOLD', '1.0', 'Minimum detection rate in target fluid'),
    ]
}

# Documentation templates
FUNCTION_DOCSTRING_TEMPLATE = '''"""
{description}

Args:
{args}

Returns:
{returns}

Raises:
{raises}

Examples:
{examples}
"""'''

CLASS_DOCSTRING_TEMPLATE = '''"""
{description}

This class {purpose}.

Attributes:
{attributes}

Methods:
{methods}

Examples:
{examples}
"""'''

# Improved docstrings for key functions
IMPROVED_DOCSTRINGS = {
    'cel_practical_forensic.py': {
        'practical_forensic_analysis': {
            'description': 'Perform practical forensic analysis on CEL data with relaxed statistical criteria',
            'purpose': 'Identifies miRNA markers suitable for forensic body fluid identification',
            'args': '''    results_dir (Path, optional): Directory containing preprocessed results
        Default: Path('results/cel_continuous_expression')''',
            'returns': '''    dict: Dictionary mapping body fluids to their marker DataFrames
        Keys: 'blood', 'saliva', 'semen', 'vaginal'
        Values: DataFrames with columns for miRNA, statistics, and rankings''',
            'notes': '''With small sample sizes (n=5), strict FDR correction is too conservative.
This function applies practical criteria focusing on effect size and
biological relevance rather than strict statistical significance.'''
        },
        'create_practical_visualizations': {
            'description': 'Create comprehensive visualizations for forensic marker analysis',
            'args': '''    markers_by_fluid (dict): Dictionary of DataFrames with markers per fluid
    output_dir (Path): Directory to save visualization files''',
            'returns': '    None: Saves PNG files to output directory',
            'outputs': '''    - top_markers_heatmap.png: Heatmap of top 3 markers per fluid
    - effect_sizes.png: Scatter plot of fold change vs Cohen\'s d
    - expression_comparison.png: Bar plots comparing expression levels'''
        }
    },
    'gpr_continuous_analysis.py': {
        'calculate_expression_tiers': {
            'description': 'Categorize miRNA expression into confidence tiers based on intensity',
            'rationale': '''Two-color arrays produce continuous ratios rather than binary 
presence/absence. This method creates interpretable categories while
preserving the continuous nature of the data.''',
            'returns': '''    pd.DataFrame: Tier assignments with columns:
        - mirna: miRNA identifier
        - fluid: body fluid type
        - expression_tier: assigned confidence level
        - log2_intensity: original log2 ratio value'''
        },
        'wilcoxon_tests': {
            'description': 'Perform Wilcoxon rank-sum tests comparing each fluid vs all others',
            'args': '''    multiple_correction (str): Correction method for multiple testing
        Options: 'bonferroni' (conservative), 'fdr' (recommended)
        Default: 'bonferroni' ''',
            'statistical_note': '''Wilcoxon test is appropriate for small samples (n=2 per fluid)
as it makes no assumptions about data distribution. However, with
such small samples, statistical power is limited.'''
        }
    }
}

# Inline comments to add
INLINE_COMMENTS = {
    'gpr_continuous_analysis.py': [
        ('line_156', 'specificity = (min_in_fluid + 0.01) / (max_in_others + 0.01)',
         '# Add small constant to prevent division by zero while maintaining ratio relationships'),
        ('line_337', 'donor_count = int(count * 2.5)',
         '# Estimate total donors: ~2.5 donors pooled per array based on experimental design'),
        ('line_65', 'score = -np.log10(row["padj"] + 1e-10) * abs(row["log2fc"]) * spec_score',
         '# Combined score weights significance, effect size, and specificity equally')
    ],
    'cel_practical_forensic.py': [
        ('line_65', "candidates['score'] = (",
         '# Scoring formula: significance × fold change × effect size\n        # Higher scores indicate better forensic markers'),
        ('line_115', 'report.append("- Testing 13,564 miRNAs requires p < 0.0000037 for FDR < 0.05")',
         '# Mathematical limit: With n=5, minimum p-value ≈ 0.001\n    # Multiple testing burden: 0.05 / 13,564 = 0.0000037')
    ]
}

def create_documentation_improvements():
    """Generate documentation improvement recommendations"""
    
    improvements = []
    
    # 1. Add module-level constants
    improvements.append("\n## 1. Add Module-Level Constants\n")
    for file, constants in CONSTANTS_TO_ADD.items():
        improvements.append(f"\n### {file}\n```python")
        improvements.append("# Add at top of file after imports:")
        improvements.append("\n# Analysis thresholds and constants")
        for name, value, description in constants:
            improvements.append(f"{name} = {value}  # {description}")
        improvements.append("```\n")
    
    # 2. Improve function docstrings
    improvements.append("\n## 2. Improve Function Docstrings\n")
    for file, functions in IMPROVED_DOCSTRINGS.items():
        improvements.append(f"\n### {file}")
        for func_name, details in functions.items():
            improvements.append(f"\n#### Function: `{func_name}`")
            improvements.append("```python")
            improvements.append(f"def {func_name}(...):")
            improvements.append('    """')
            improvements.append(f"    {details['description']}")
            improvements.append("")
            if 'purpose' in details:
                improvements.append(f"    {details['purpose']}")
                improvements.append("")
            if 'args' in details:
                improvements.append("    Args:")
                improvements.append(details['args'])
                improvements.append("")
            if 'returns' in details:
                improvements.append("    Returns:")
                improvements.append(details['returns'])
                improvements.append("")
            if 'notes' in details:
                improvements.append("    Note:")
                improvements.append(f"        {details['notes']}")
            improvements.append('    """')
            improvements.append("```")
    
    # 3. Add inline comments
    improvements.append("\n## 3. Add Inline Comments\n")
    for file, comments in INLINE_COMMENTS.items():
        improvements.append(f"\n### {file}")
        for line_ref, code, comment in comments:
            improvements.append(f"\n{line_ref}:")
            improvements.append("```python")
            improvements.append(comment)
            improvements.append(code)
            improvements.append("```")
    
    # 4. Add type hints
    improvements.append("\n## 4. Add Type Hints\n")
    improvements.append("\n### All Python files")
    improvements.append("```python")
    improvements.append("from typing import Dict, List, Tuple, Optional, Union")
    improvements.append("import pandas as pd")
    improvements.append("from pathlib import Path")
    improvements.append("")
    improvements.append("# Example function signatures with type hints:")
    improvements.append("def load_data(self) -> None:")
    improvements.append("def calculate_expression_tiers(self) -> pd.DataFrame:")
    improvements.append("def wilcoxon_tests(self, multiple_correction: str = 'bonferroni') -> pd.DataFrame:")
    improvements.append("def create_visualizations(self, output_dir: Path) -> None:")
    improvements.append("```")
    
    # 5. R documentation
    improvements.append("\n## 5. Add R Documentation\n")
    improvements.append("\n### cel_processor_minimal.R")
    improvements.append("```r")
    improvements.append("#' Process Affymetrix miRNA ST Array CEL Files")
    improvements.append("#'")
    improvements.append("#' @description")
    improvements.append("#' Reads CEL files from GSE49630 and performs RMA normalization")
    improvements.append("#' using the oligo package. The oligo package is required for")
    improvements.append("#' ST arrays as the older affy package does not support them.")
    improvements.append("#'")
    improvements.append("#' @param cel_dir Character. Directory containing .CEL.gz files")
    improvements.append("#' @param output_dir Character. Directory for output files")
    improvements.append("#'")
    improvements.append("#' @details")
    improvements.append("#' This function processes Affymetrix Human Gene 1.0 ST arrays")
    improvements.append("#' which require specialized handling compared to 3' IVT arrays.")
    improvements.append("#'")
    improvements.append("#' @return")
    improvements.append("#' Invisible NULL. Writes two files:")
    improvements.append("#' \\itemize{")
    improvements.append("#'   \\item cel_expression_matrix.csv: Normalized expression values")
    improvements.append("#'   \\item cel_metadata.csv: Sample annotations")
    improvements.append("#' }")
    improvements.append("#'")
    improvements.append("#' @examples")
    improvements.append("#' \\dontrun{")
    improvements.append("#' process_cel_files('data/GSE49630_RAW', 'data/processed/cel')")
    improvements.append("#' }")
    improvements.append("#'")
    improvements.append("#' @export")
    improvements.append("process_cel_files <- function(cel_dir, output_dir) {")
    improvements.append("```")
    
    # 6. Clean up unused files
    improvements.append("\n## 6. Files to Remove or Archive\n")
    improvements.append("\nMove these abandoned/test files to an `archive/` directory:")
    improvements.append("- `scripts/preprocessing/cel_processor.R` (replaced by minimal version)")
    improvements.append("- `scripts/preprocessing/cel_processor_simple.py` (incomplete)")
    improvements.append("- `scripts/preprocessing/cel_processor_rpy2.py` (abandoned approach)")
    improvements.append("- `scripts/preprocessing/cel_processor_final.py` (another attempt)")
    improvements.append("- `minimal_cel_process.R` (root directory duplicate)")
    improvements.append("- `process_affymetrix_mirna.R` (alternative approach)")
    
    # 7. Create README for preprocessing
    improvements.append("\n## 7. Create Preprocessing README\n")
    improvements.append("\n### scripts/preprocessing/README.md")
    improvements.append("```markdown")
    improvements.append("# Preprocessing Scripts")
    improvements.append("")
    improvements.append("## Overview")
    improvements.append("This directory contains scripts for processing raw microarray data from two platforms:")
    improvements.append("- GenePix Result (GPR) files from two-color Agilent arrays")
    improvements.append("- CEL files from Affymetrix miRNA ST arrays")
    improvements.append("")
    improvements.append("## Production Scripts")
    improvements.append("")
    improvements.append("### gpr_parser.py")
    improvements.append("Processes GPR files from GSE153135 (pooled body fluid samples)")
    improvements.append("- Input: .gpr or .gpr.gz files")
    improvements.append("- Output: Normalized log2 ratios (Cy5/Cy3)")
    improvements.append("- Key feature: Background correction")
    improvements.append("")
    improvements.append("### cel_processor_minimal.R")
    improvements.append("Processes CEL files from GSE49630 (individual samples)")
    improvements.append("- Input: .CEL.gz files")
    improvements.append("- Output: RMA-normalized expression matrix")
    improvements.append("- Requirement: oligo package (NOT affy)")
    improvements.append("")
    improvements.append("## Requirements")
    improvements.append("")
    improvements.append("### Python")
    improvements.append("- Python 3.8+")
    improvements.append("- pandas, numpy")
    improvements.append("- See requirements.txt")
    improvements.append("")
    improvements.append("### R")
    improvements.append("- R 4.5+")
    improvements.append("- BiocManager::install('oligo')")
    improvements.append("- BiocManager::install('limma')")
    improvements.append("")
    improvements.append("## Usage Examples")
    improvements.append("")
    improvements.append("```bash")
    improvements.append("# Process GPR files")
    improvements.append("python gpr_parser.py \\")
    improvements.append("    --input data/GSE153135_RAW \\")
    improvements.append("    --output data/processed/gpr")
    improvements.append("")
    improvements.append("# Process CEL files")
    improvements.append("Rscript cel_processor_minimal.R \\")
    improvements.append("    data/GSE49630_RAW \\")
    improvements.append("    data/processed/cel")
    improvements.append("```")
    improvements.append("```")
    
    return '\n'.join(improvements)

if __name__ == "__main__":
    # Generate documentation improvements
    improvements = create_documentation_improvements()
    
    # Save to file
    output_path = Path('DOCUMENTATION_IMPROVEMENTS.md')
    with open(output_path, 'w') as f:
        f.write("# Documentation Improvements Guide\n")
        f.write("\nThis guide provides specific improvements for code documentation.\n")
        f.write(improvements)
    
    print(f"Documentation improvements saved to {output_path}")
    
    # Also create a cleaner script
    print("\nCreating file cleanup script...")
    
    cleanup_script = '''#!/bin/bash
# Script to clean up redundant preprocessing files

echo "Creating archive directory..."
mkdir -p scripts/preprocessing/archive

echo "Moving abandoned files to archive..."
mv scripts/preprocessing/cel_processor.R scripts/preprocessing/archive/ 2>/dev/null || true
mv scripts/preprocessing/cel_processor_simple.py scripts/preprocessing/archive/ 2>/dev/null || true
mv scripts/preprocessing/cel_processor_rpy2.py scripts/preprocessing/archive/ 2>/dev/null || true
mv scripts/preprocessing/cel_processor_final.py scripts/preprocessing/archive/ 2>/dev/null || true
mv minimal_cel_process.R scripts/preprocessing/archive/ 2>/dev/null || true
mv process_affymetrix_mirna.R scripts/preprocessing/archive/ 2>/dev/null || true

echo "Cleanup complete! Abandoned files moved to scripts/preprocessing/archive/"
'''
    
    with open('cleanup_preprocessing.sh', 'w') as f:
        f.write(cleanup_script)
    
    os.chmod('cleanup_preprocessing.sh', 0o755)
    print("Created cleanup_preprocessing.sh")