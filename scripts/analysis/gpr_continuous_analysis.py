#!/usr/bin/env python3
"""
Continuous Expression Analysis for GPR Data
Implements multi-tier detection framework and proper statistics
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from itertools import combinations

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ContinuousExpressionAnalyzer:
    """Analyze miRNA expression as continuous variables"""
    
    def __init__(self, expression_path: Path):
        self.expression_path = expression_path
        self.expr_data = None
        self.metadata = None
        
        # Multi-tier thresholds (log2 scale)
        self.thresholds = {
            'high_confidence': -2.0,
            'moderate_confidence': -6.0,
            'low_confidence': -10.0
        }
        
    def load_data(self):
        """Load expression data and extract metadata"""
        logger.info("Loading expression data...")
        self.expr_data = pd.read_csv(self.expression_path, index_col=0)
        
        # Extract metadata from sample names
        self.metadata = pd.DataFrame({
            'sample': self.expr_data.index,
            'fluid_type': [self._extract_fluid_type(s) for s in self.expr_data.index]
        })
        
        logger.info(f"Loaded {len(self.expr_data)} samples, {len(self.expr_data.columns)} miRNAs")
        logger.info(f"Fluid distribution: {self.metadata['fluid_type'].value_counts().to_dict()}")
        
    def _extract_fluid_type(self, sample_name):
        """Extract fluid type from sample name"""
        if 'saliva' in sample_name:
            return 'saliva'
        elif 'peripheral_blood' in sample_name:
            return 'blood'
        elif 'semen' in sample_name:
            return 'semen'
        elif 'vaginal' in sample_name:
            return 'vaginal'
        elif 'menstrual' in sample_name:
            return 'menstrual'
        else:
            return 'unknown'
            
    def calculate_expression_tiers(self):
        """Assign expression values to confidence tiers"""
        tiers = pd.DataFrame(index=self.expr_data.index, columns=self.expr_data.columns)
        
        for col in self.expr_data.columns:
            values = self.expr_data[col]
            tiers[col] = pd.cut(values, 
                               bins=[-np.inf, self.thresholds['low_confidence'], 
                                     self.thresholds['moderate_confidence'], 
                                     self.thresholds['high_confidence'], np.inf],
                               labels=['below_detection', 'low', 'moderate', 'high'])
        
        return tiers
        
    def wilcoxon_tests(self, multiple_correction='bonferroni'):
        """Perform Wilcoxon tests for each fluid vs all others"""
        results = []
        
        for fluid in self.metadata['fluid_type'].unique():
            if fluid == 'unknown':
                continue
                
            fluid_samples = self.metadata[self.metadata['fluid_type'] == fluid]['sample'].tolist()
            other_samples = self.metadata[self.metadata['fluid_type'] != fluid]['sample'].tolist()
            
            # Skip if too few samples
            if len(fluid_samples) < 2 or len(other_samples) < 2:
                logger.warning(f"Skipping {fluid} - too few samples")
                continue
                
            for mirna in self.expr_data.columns:
                fluid_expr = self.expr_data.loc[fluid_samples, mirna].values
                other_expr = self.expr_data.loc[other_samples, mirna].values
                
                # Wilcoxon test
                try:
                    stat, pval = stats.ranksums(fluid_expr, other_expr)
                    
                    # Calculate effect size (Cohen's d)
                    pooled_std = np.sqrt((np.var(fluid_expr) + np.var(other_expr)) / 2)
                    if pooled_std > 0:
                        cohens_d = (np.mean(fluid_expr) - np.mean(other_expr)) / pooled_std
                    else:
                        cohens_d = 0
                        
                    # Log2 fold change
                    log2fc = np.mean(fluid_expr) - np.mean(other_expr)
                    
                    results.append({
                        'fluid': fluid,
                        'mirna': mirna,
                        'pvalue': pval,
                        'log2fc': log2fc,
                        'cohens_d': cohens_d,
                        'fluid_mean': np.mean(fluid_expr),
                        'other_mean': np.mean(other_expr),
                        'fluid_std': np.std(fluid_expr),
                        'other_std': np.std(other_expr)
                    })
                    
                except Exception as e:
                    logger.warning(f"Test failed for {fluid} - {mirna}: {e}")
                    
        results_df = pd.DataFrame(results)
        
        # Multiple testing correction
        if multiple_correction == 'bonferroni':
            results_df['padj'] = results_df['pvalue'] * len(results_df)
            results_df['padj'] = results_df['padj'].clip(upper=1.0)
        elif multiple_correction == 'fdr':
            from statsmodels.stats.multitest import multipletests
            _, results_df['padj'], _, _ = multipletests(results_df['pvalue'], method='fdr_bh')
            
        return results_df
        
    def calculate_specificity_scores(self):
        """Calculate specificity scores for each miRNA/fluid combination"""
        scores = []
        
        for fluid in self.metadata['fluid_type'].unique():
            if fluid == 'unknown':
                continue
                
            fluid_samples = self.metadata[self.metadata['fluid_type'] == fluid]['sample'].tolist()
            other_samples = self.metadata[self.metadata['fluid_type'] != fluid]['sample'].tolist()
            
            for mirna in self.expr_data.columns:
                fluid_expr = self.expr_data.loc[fluid_samples, mirna].values
                other_expr = self.expr_data.loc[other_samples, mirna].values
                
                # Specificity score: min(fluid) / max(others)
                # Add small constant to avoid division by zero
                min_fluid = np.min(fluid_expr) + 0.01
                max_other = np.max(other_expr) + 0.01
                
                if min_fluid > 0 and max_other > 0:
                    spec_score = min_fluid / max_other
                else:
                    spec_score = 0
                    
                # Detection rate in target fluid
                detection_rate = np.mean(fluid_expr > self.thresholds['low_confidence'])
                
                scores.append({
                    'fluid': fluid,
                    'mirna': mirna,
                    'specificity_score': spec_score,
                    'detection_rate': detection_rate,
                    'fluid_min': np.min(fluid_expr),
                    'fluid_max': np.max(fluid_expr),
                    'other_max': np.max(other_expr)
                })
                
        return pd.DataFrame(scores)
        
    def find_forensic_markers(self, results_df, spec_df):
        """Identify forensically useful markers based on multiple criteria"""
        # Merge statistical and specificity results
        merged = results_df.merge(spec_df, on=['fluid', 'mirna'])
        
        # Apply forensic criteria
        forensic_markers = merged[
            (merged['padj'] < 0.05) &  # Statistically significant
            (np.abs(merged['log2fc']) > 2.32) &  # >5-fold change
            (np.abs(merged['cohens_d']) > 0.8) &  # Large effect size
            (merged['detection_rate'] > 0.8) &  # Detected in >80% of target
            (merged['specificity_score'] > 5)  # High specificity
        ].copy()
        
        # Sort by combined score
        forensic_markers['combined_score'] = (
            np.abs(forensic_markers['log2fc']) * 
            np.abs(forensic_markers['cohens_d']) * 
            forensic_markers['specificity_score']
        )
        
        forensic_markers = forensic_markers.sort_values('combined_score', ascending=False)
        
        return forensic_markers
        
    def create_visualizations(self, output_dir: Path):
        """Create comprehensive visualizations"""
        output_dir.mkdir(exist_ok=True)
        
        # 1. Expression distribution by confidence tier
        fig, ax = plt.subplots(figsize=(10, 6))
        
        all_values = self.expr_data.values.flatten()
        
        # Plot distribution with tier boundaries
        ax.hist(all_values, bins=50, alpha=0.7, color='blue', edgecolor='black')
        
        # Add tier boundaries
        colors = ['red', 'orange', 'green']
        labels = ['Low conf.', 'Moderate conf.', 'High conf.']
        
        for threshold, color, label in zip(self.thresholds.values(), colors, labels):
            ax.axvline(threshold, color=color, linestyle='--', linewidth=2, label=f'{label} (>{threshold})')
            
        ax.set_xlabel('Log2 Expression')
        ax.set_ylabel('Frequency')
        ax.set_title('Expression Distribution with Confidence Tiers')
        ax.legend()
        
        plt.savefig(output_dir / 'expression_tiers.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 2. Fluid-specific expression patterns
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()
        
        for i, fluid in enumerate(self.metadata['fluid_type'].unique()):
            if fluid == 'unknown' or i >= 6:
                continue
                
            ax = axes[i]
            
            fluid_samples = self.metadata[self.metadata['fluid_type'] == fluid]['sample'].tolist()
            
            # Get mean expression for this fluid
            fluid_means = self.expr_data.loc[fluid_samples].mean(axis=0)
            
            # Plot histogram
            ax.hist(fluid_means, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
            ax.set_xlabel('Mean Log2 Expression')
            ax.set_ylabel('Number of miRNAs')
            ax.set_title(f'{fluid.capitalize()} (n={len(fluid_samples)})')
            
            # Add threshold lines
            for threshold, color in zip(self.thresholds.values(), colors):
                ax.axvline(threshold, color=color, linestyle='--', alpha=0.5)
                
        # Remove extra subplots
        for i in range(len(self.metadata['fluid_type'].unique()), 6):
            fig.delaxes(axes[i])
            
        plt.suptitle('Expression Distributions by Body Fluid')
        plt.tight_layout()
        plt.savefig(output_dir / 'fluid_distributions.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 3. Sample quality assessment
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Mean expression per sample
        sample_means = self.expr_data.mean(axis=1)
        sample_stds = self.expr_data.std(axis=1)
        
        colors_map = {'saliva': 'blue', 'blood': 'red', 'semen': 'green', 
                     'vaginal': 'purple', 'menstrual': 'orange'}
        
        for fluid in self.metadata['fluid_type'].unique():
            if fluid == 'unknown':
                continue
            fluid_samples = self.metadata[self.metadata['fluid_type'] == fluid]['sample'].tolist()
            x = sample_means[fluid_samples]
            y = sample_stds[fluid_samples]
            ax1.scatter(x, y, label=fluid, color=colors_map.get(fluid, 'gray'), s=100, alpha=0.7)
            
        ax1.set_xlabel('Mean Expression')
        ax1.set_ylabel('Standard Deviation')
        ax1.set_title('Sample Quality: Mean vs Variability')
        ax1.legend()
        
        # Detection rates per sample
        detection_rates = (self.expr_data > self.thresholds['low_confidence']).sum(axis=1) / len(self.expr_data.columns)
        
        # Sort by fluid type for visualization
        sorted_samples = []
        for fluid in ['saliva', 'blood', 'semen', 'vaginal', 'menstrual']:
            fluid_samples = self.metadata[self.metadata['fluid_type'] == fluid]['sample'].tolist()
            sorted_samples.extend(fluid_samples)
            
        y_pos = range(len(sorted_samples))
        colors_list = []
        for sample in sorted_samples:
            fluid = self.metadata[self.metadata['sample'] == sample]['fluid_type'].values[0]
            colors_list.append(colors_map.get(fluid, 'gray'))
            
        ax2.barh(y_pos, detection_rates[sorted_samples], color=colors_list, alpha=0.7)
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels([s.split('_')[0] for s in sorted_samples], fontsize=8)
        ax2.set_xlabel('Detection Rate')
        ax2.set_title('miRNA Detection Rate by Sample')
        ax2.axvline(0.5, color='red', linestyle='--', alpha=0.5)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'sample_quality.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Visualizations saved to {output_dir}")
        
    def generate_report(self, output_dir: Path):
        """Generate comprehensive analysis report"""
        # Run analyses
        tiers = self.calculate_expression_tiers()
        wilcox_results = self.wilcoxon_tests()
        spec_scores = self.calculate_specificity_scores()
        forensic_markers = self.find_forensic_markers(wilcox_results, spec_scores)
        
        # Create visualizations
        self.create_visualizations(output_dir)
        
        # Generate report
        report = []
        report.append("# Continuous Expression Analysis Report (GPR Data)")
        report.append("\n## Dataset Overview")
        report.append(f"- Total samples: {len(self.expr_data)}")
        report.append(f"- Total miRNAs: {len(self.expr_data.columns)}")
        report.append(f"- Expression range: {self.expr_data.min().min():.2f} to {self.expr_data.max().max():.2f} (log2)")
        
        report.append("\n## Sample Distribution")
        for fluid, count in self.metadata['fluid_type'].value_counts().items():
            report.append(f"- {fluid}: {count} samples (pooled from {count*2.5:.0f} donors)")
            
        report.append("\n## Expression Tier Distribution")
        tier_counts = tiers.apply(lambda x: x.value_counts()).sum(axis=1)
        total_measurements = tier_counts.sum()
        report.append(f"- High confidence (>{self.thresholds['high_confidence']}): {tier_counts['high']:.0f} ({tier_counts['high']/total_measurements*100:.1f}%)")
        report.append(f"- Moderate confidence: {tier_counts['moderate']:.0f} ({tier_counts['moderate']/total_measurements*100:.1f}%)")
        report.append(f"- Low confidence: {tier_counts['low']:.0f} ({tier_counts['low']/total_measurements*100:.1f}%)")
        report.append(f"- Below detection: {tier_counts['below_detection']:.0f} ({tier_counts['below_detection']/total_measurements*100:.1f}%)")
        
        report.append("\n## Statistical Analysis Results")
        report.append(f"- Total comparisons: {len(wilcox_results)}")
        report.append(f"- Significant (p<0.05 before correction): {(wilcox_results['pvalue'] < 0.05).sum()}")
        report.append(f"- Significant (Bonferroni corrected): {(wilcox_results['padj'] < 0.05).sum()}")
        
        report.append("\n## Forensic Marker Candidates")
        if len(forensic_markers) > 0:
            report.append(f"\nFound {len(forensic_markers)} candidates meeting all criteria:")
            
            for fluid in forensic_markers['fluid'].unique():
                fluid_markers = forensic_markers[forensic_markers['fluid'] == fluid]
                report.append(f"\n### {fluid.capitalize()}")
                
                for _, marker in fluid_markers.head(3).iterrows():
                    report.append(f"- **{marker['mirna']}**")
                    report.append(f"  - Log2 fold change: {marker['log2fc']:.2f}")
                    report.append(f"  - Cohen's d: {marker['cohens_d']:.2f}")
                    report.append(f"  - Specificity score: {marker['specificity_score']:.2f}")
                    report.append(f"  - Detection rate: {marker['detection_rate']*100:.0f}%")
                    report.append(f"  - p-value (adj): {marker['padj']:.3e}")
        else:
            report.append("\n**WARNING**: No markers met all forensic criteria.")
            report.append("This is expected with n=2 pooled samples per fluid.")
            report.append("\nRelaxed criteria results (p<0.05, |log2FC|>1):")
            
            relaxed = wilcox_results[
                (wilcox_results['pvalue'] < 0.05) & 
                (np.abs(wilcox_results['log2fc']) > 1)
            ].sort_values('pvalue')
            
            for fluid in relaxed['fluid'].unique()[:3]:
                fluid_results = relaxed[relaxed['fluid'] == fluid]
                report.append(f"\n### {fluid.capitalize()}")
                for _, result in fluid_results.head(3).iterrows():
                    report.append(f"- {result['mirna']}: log2FC={result['log2fc']:.2f}, p={result['pvalue']:.3f}")
                    
        report.append("\n## Limitations")
        report.append("- Small sample size (n=2 arrays per fluid) limits statistical power")
        report.append("- Arrays represent pooled samples, reducing biological variability")
        report.append("- Results require validation in larger, individual sample cohorts")
        report.append("- CEL data analysis recommended for more robust findings")
        
        # Save report
        report_path = output_dir / 'continuous_analysis_report.md'
        with open(report_path, 'w') as f:
            f.write('\n'.join(report))
            
        # Save detailed results
        wilcox_results.to_csv(output_dir / 'wilcoxon_results.csv', index=False)
        spec_scores.to_csv(output_dir / 'specificity_scores.csv', index=False)
        if len(forensic_markers) > 0:
            forensic_markers.to_csv(output_dir / 'forensic_markers.csv', index=False)
            
        logger.info(f"Report saved to {report_path}")
        
        return forensic_markers


def main():
    """Run continuous expression analysis on GPR data"""
    # Paths
    expression_path = Path('data/processed/gpr/gpr_expression_matrix.csv')
    output_dir = Path('results/continuous_expression_analysis')
    
    # Run analysis
    analyzer = ContinuousExpressionAnalyzer(expression_path)
    analyzer.load_data()
    forensic_markers = analyzer.generate_report(output_dir)
    
    logger.info("Analysis complete!")
    

if __name__ == "__main__":
    main()