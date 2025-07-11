#!/usr/bin/env python3
"""
Continuous Expression Analysis for CEL Data
Applies the same framework we developed for GPR to CEL data
Now with proper sample sizes (n=5 per fluid)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import sys

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent.parent))

# Import our existing analyzer
from scripts.analysis.gpr_continuous_analysis import ContinuousExpressionAnalyzer

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class CELContinuousAnalyzer(ContinuousExpressionAnalyzer):
    """Extends GPR analyzer for CEL data specifics"""
    
    def load_data(self):
        """Load CEL expression data"""
        logger.info("Loading CEL expression data...")
        
        # Read expression matrix
        expr_df = pd.read_csv(self.expression_path, index_col=0)
        
        # Transpose to have samples as rows
        self.expr_data = expr_df.T
        
        # Extract metadata from sample names
        self.metadata = pd.DataFrame({
            'sample': self.expr_data.index,
            'fluid_type': [self._extract_fluid_type_cel(s) for s in self.expr_data.index]
        })
        
        logger.info(f"Loaded {len(self.expr_data)} samples, {len(self.expr_data.columns)} probes")
        logger.info(f"Fluid distribution: {self.metadata['fluid_type'].value_counts().to_dict()}")
        
    def _extract_fluid_type_cel(self, sample_name):
        """Extract fluid type from CEL sample names"""
        if 'Blood' in sample_name:
            return 'blood'
        elif 'Semen' in sample_name:
            return 'semen'
        elif 'Vaginal' in sample_name:
            return 'vaginal'
        elif 'Saliva' in sample_name:
            return 'saliva'
        else:
            return 'unknown'
            
    def filter_to_mirnas(self):
        """Filter to only human miRNA probes"""
        # Keep only probes that look like miRNAs
        mirna_cols = [col for col in self.expr_data.columns 
                     if ('mir' in col.lower() or 'let-' in col.lower()) 
                     and 'st' in col.lower()]
        
        logger.info(f"Filtering from {len(self.expr_data.columns)} to {len(mirna_cols)} miRNA probes")
        self.expr_data = self.expr_data[mirna_cols]
        
    def generate_report(self, output_dir: Path):
        """Generate comprehensive analysis report with CEL-specific details"""
        # Filter to miRNAs first
        self.filter_to_mirnas()
        
        # Call parent method
        forensic_markers = super().generate_report(output_dir)
        
        # Add CEL-specific summary
        with open(output_dir / 'continuous_analysis_report.md', 'a') as f:
            f.write("\n## CEL-Specific Notes\n")
            f.write("- Platform: Affymetrix miRNA 3.0 ST arrays\n")
            f.write("- Processing: oligo::rma() normalization\n")
            f.write("- Sample type: Individual donors (not pooled)\n")
            f.write("- Statistical power: Improved with n=5 per fluid\n")
            
        return forensic_markers


def main():
    """Run continuous expression analysis on CEL data"""
    # Paths
    expression_path = Path('data/processed/cel/cel_expression_matrix.csv')
    output_dir = Path('results/cel_continuous_expression')
    
    # Run analysis
    analyzer = CELContinuousAnalyzer(expression_path)
    analyzer.load_data()
    forensic_markers = analyzer.generate_report(output_dir)
    
    if len(forensic_markers) > 0:
        logger.info(f"\nFound {len(forensic_markers)} forensic marker candidates!")
        logger.info("With n=5 per fluid, these results are more statistically reliable than GPR.")
    else:
        logger.info("\nNo markers met all forensic criteria.")
        logger.info("This is unexpected with n=5. Check relaxed criteria results.")
    
    logger.info(f"\nAnalysis complete! Results in {output_dir}")
    

if __name__ == "__main__":
    main()