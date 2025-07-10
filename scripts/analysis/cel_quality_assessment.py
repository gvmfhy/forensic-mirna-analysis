#!/usr/bin/env python3
"""
CEL Data Quality Assessment
Focus on understanding data structure before analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class CELQualityAssessment:
    """Assess CEL data quality and structure"""
    
    def __init__(self, cel_expr_path: Path, annotation_path: Path):
        self.cel_expr_path = cel_expr_path
        self.annotation_path = annotation_path
        self.expr_data = None
        self.metadata = None
        self.annotation = None
        
    def load_data(self):
        """Load processed CEL expression data if available"""
        if not self.cel_expr_path.exists():
            logger.warning(f"CEL expression file not found at {self.cel_expr_path}")
            logger.info("Please run R processing script first:")
            logger.info("cd data/processed/cel && R --vanilla < process_cel_files.R")
            return False
            
        # Load expression data
        self.expr_data = pd.read_csv(self.cel_expr_path, index_col=0)
        logger.info(f"Loaded expression data: {self.expr_data.shape}")
        
        # Load metadata
        metadata_path = self.cel_expr_path.parent / 'cel_metadata.csv'
        if metadata_path.exists():
            self.metadata = pd.read_csv(metadata_path)
            logger.info(f"Loaded metadata: {self.metadata.shape}")
        
        return True
        
    def assess_detection_rates(self):
        """Calculate detection rates for each miRNA"""
        if self.expr_data is None:
            return None
            
        # Define detection threshold (log2 expression > 0 means above background)
        detection_threshold = 0
        
        # Calculate detection rate per miRNA
        detection_rates = (self.expr_data > detection_threshold).sum(axis=1) / self.expr_data.shape[1]
        
        # Summarize
        summary = {
            'total_mirnas': len(detection_rates),
            'always_detected': (detection_rates == 1.0).sum(),
            'never_detected': (detection_rates == 0.0).sum(),
            'sometimes_detected': ((detection_rates > 0) & (detection_rates < 1)).sum(),
            'detected_in_majority': (detection_rates > 0.5).sum()
        }
        
        return detection_rates, summary
        
    def assess_expression_ranges(self):
        """Understand the dynamic range of expression"""
        if self.expr_data is None:
            return None
            
        # Get expression statistics
        expr_stats = {
            'min': self.expr_data.min().min(),
            'max': self.expr_data.max().max(),
            'mean': self.expr_data.mean().mean(),
            'median': self.expr_data.median().median(),
            'std': self.expr_data.std().mean()
        }
        
        # Find miRNAs with highest dynamic range
        mirna_ranges = self.expr_data.max(axis=1) - self.expr_data.min(axis=1)
        top_dynamic = mirna_ranges.nlargest(20)
        
        return expr_stats, top_dynamic
        
    def assess_sample_quality(self):
        """Check quality metrics per sample"""
        if self.expr_data is None:
            return None
            
        sample_stats = pd.DataFrame({
            'mean_expression': self.expr_data.mean(axis=0),
            'detected_mirnas': (self.expr_data > 0).sum(axis=0),
            'strong_signals': (self.expr_data > 5).sum(axis=0),  # Log2 > 5 is strong
            'sample': self.expr_data.columns
        })
        
        # Add fluid type if metadata available
        if self.metadata is not None:
            sample_stats = sample_stats.merge(self.metadata[['sample', 'fluid_type']], on='sample')
            
        return sample_stats
        
    def find_fluid_specific_candidates(self):
        """Initial search for fluid-specific patterns without ML"""
        if self.expr_data is None or self.metadata is None:
            return None
            
        results = {}
        
        for fluid in self.metadata['fluid_type'].unique():
            # Get samples for this fluid
            fluid_samples = self.metadata[self.metadata['fluid_type'] == fluid]['sample'].tolist()
            other_samples = self.metadata[self.metadata['fluid_type'] != fluid]['sample'].tolist()
            
            # Calculate mean expression
            fluid_expr = self.expr_data[fluid_samples].mean(axis=1)
            other_expr = self.expr_data[other_samples].mean(axis=1)
            
            # Find miRNAs high in this fluid, low in others
            # Forensic criteria: >5-fold difference
            fold_change = fluid_expr - other_expr  # Log2 scale
            
            # Additional criteria: must be reliably detected in target fluid
            detection_in_fluid = (self.expr_data[fluid_samples] > 0).mean(axis=1)
            
            # Candidates: >5-fold (2.32 in log2) and detected in >80% of fluid samples
            candidates = (fold_change > 2.32) & (detection_in_fluid > 0.8)
            
            results[fluid] = {
                'candidates': fold_change[candidates].sort_values(ascending=False),
                'n_candidates': candidates.sum()
            }
            
        return results
        
    def create_quality_report(self, output_dir: Path):
        """Generate comprehensive quality report"""
        output_dir.mkdir(exist_ok=True)
        
        # Detection rates
        detection_rates, detection_summary = self.assess_detection_rates()
        
        # Expression ranges  
        expr_stats, top_dynamic = self.assess_expression_ranges()
        
        # Sample quality
        sample_stats = self.assess_sample_quality()
        
        # Fluid-specific candidates
        candidates = self.find_fluid_specific_candidates()
        
        # Create report
        report = []
        report.append("# CEL Data Quality Assessment Report\n")
        report.append(f"## Data Overview")
        report.append(f"- Total miRNAs: {self.expr_data.shape[0]}")
        report.append(f"- Total samples: {self.expr_data.shape[1]}")
        report.append(f"- Expression range: {expr_stats['min']:.2f} to {expr_stats['max']:.2f} (log2)\n")
        
        report.append(f"## Detection Summary")
        report.append(f"- Always detected: {detection_summary['always_detected']} miRNAs")
        report.append(f"- Never detected: {detection_summary['never_detected']} miRNAs") 
        report.append(f"- Sometimes detected: {detection_summary['sometimes_detected']} miRNAs")
        report.append(f"- Detected in majority: {detection_summary['detected_in_majority']} miRNAs\n")
        
        report.append(f"## Sample Quality")
        for fluid in sample_stats['fluid_type'].unique():
            fluid_data = sample_stats[sample_stats['fluid_type'] == fluid]
            report.append(f"\n### {fluid}")
            report.append(f"- Mean detected miRNAs: {fluid_data['detected_mirnas'].mean():.0f}")
            report.append(f"- Mean strong signals: {fluid_data['strong_signals'].mean():.0f}")
            
        report.append(f"\n## Forensic Marker Candidates")
        for fluid, data in candidates.items():
            report.append(f"\n### {fluid}")
            report.append(f"- Candidates found: {data['n_candidates']}")
            if data['n_candidates'] > 0:
                report.append(f"- Top 5:")
                for mirna, fc in data['candidates'].head().items():
                    report.append(f"  - {mirna}: {fc:.2f} log2 fold change")
                    
        # Save report
        with open(output_dir / 'cel_quality_report.md', 'w') as f:
            f.write('\n'.join(report))
            
        logger.info(f"Quality report saved to {output_dir / 'cel_quality_report.md'}")
        
        return detection_summary, expr_stats, sample_stats, candidates


def main():
    """Run CEL quality assessment"""
    # Paths
    cel_expr_path = Path('data/processed/cel/cel_expression_matrix.csv')
    annotation_path = Path('data/raw/GSE49630_CEL/GPL16384_miRNA-3_1-st-v1.annotations.20140513.csv.gz')
    output_dir = Path('results/cel_quality_assessment')
    
    # Run assessment
    qa = CELQualityAssessment(cel_expr_path, annotation_path)
    
    if qa.load_data():
        qa.create_quality_report(output_dir)
        logger.info("Quality assessment complete")
    else:
        logger.error("Could not load CEL data. Please process CEL files first.")
        

if __name__ == "__main__":
    main()