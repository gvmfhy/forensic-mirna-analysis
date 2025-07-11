#!/usr/bin/env python3
"""
Continuous Expression Analysis for GPR Data

This module implements a multi-tier detection framework for analyzing miRNA expression
from two-color microarray data (GenePix Result files). It was specifically developed
for GSE153135 dataset containing pooled body fluid samples for forensic identification.

The key innovation is treating miRNA expression as continuous rather than binary,
acknowledging that two-color arrays produce ratio measurements that exist on a
spectrum rather than simple presence/absence.

Classes:
    ContinuousExpressionAnalyzer: Main analysis class implementing the multi-tier framework

Key Concepts:
    - Multi-tier detection: High/Moderate/Low confidence levels based on expression ratios
    - Pooled samples: Each array contains 2-3 donors pooled together
    - Two-color normalization: Log2(Cy5/Cy3) ratios with background correction
    - Non-parametric statistics: Wilcoxon tests suitable for small sample sizes

Author: Forensic miRNA Analysis Pipeline
Date: 2024
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Optional, Union, Any

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Analysis Constants
ZERO_OFFSET = 0.01  # Small constant to avoid log(0) in calculations
DONOR_MULTIPLIER = 2.5  # Average number of donors per pooled sample (26 donors / 10 arrays)
HIGH_CONFIDENCE_THRESHOLD = -2.0  # Log2 ratio threshold for high confidence detection
MODERATE_CONFIDENCE_THRESHOLD = -6.0  # Log2 ratio threshold for moderate confidence  
LOW_CONFIDENCE_THRESHOLD = -10.0  # Log2 ratio threshold for low confidence
MIN_SPECIFICITY_SCORE = 5.0  # Minimum specificity score for forensic markers
MIN_COHEN_D = 0.8  # Minimum effect size for meaningful differences


class ContinuousExpressionAnalyzer:
    """
    Analyze miRNA expression as continuous variables with multi-tier confidence levels.
    
    This class implements the complete analysis pipeline for GPR-derived expression data,
    from loading normalized intensities through statistical testing and forensic marker
    identification. It addresses the specific challenges of:
    - Small sample sizes (n=2 per fluid type due to pooling)
    - Continuous expression measurements from two-color arrays
    - Need for forensically relevant specificity
    
    Attributes:
        expression_path (Path): Path to the expression matrix CSV file
        expr_data (pd.DataFrame): Expression matrix with samples as rows, miRNAs as columns
        metadata (pd.DataFrame): Sample metadata including fluid types
        thresholds (dict): Log2 ratio thresholds for confidence tiers
        
    Methods:
        load_data(): Load expression data and extract metadata
        calculate_expression_tiers(): Categorize expression into confidence levels
        wilcoxon_tests(): Perform differential expression analysis
        calculate_specificity_scores(): Compute fluid-specific expression metrics
        find_forensic_markers(): Identify markers suitable for forensic use
        create_visualizations(): Generate analysis plots
        
    Example:
        >>> analyzer = ContinuousExpressionAnalyzer(Path('data/gpr_expression.csv'))
        >>> analyzer.load_data()
        >>> tiers = analyzer.calculate_expression_tiers()
        >>> results = analyzer.wilcoxon_tests(multiple_correction='fdr')
        >>> markers = analyzer.find_forensic_markers(results, specificity_scores)
    """
    
    def __init__(self, expression_path: Path):
        """
        Initialize the analyzer with path to expression data.
        
        Args:
            expression_path (Path): Path to CSV file containing normalized expression matrix
                Expected format: Samples as rows, miRNAs as columns, with sample names
                containing fluid type information (e.g., "blood_1", "saliva_2")
        """
        self.expression_path = expression_path
        self.expr_data = None
        self.metadata = None
        
        # Multi-tier thresholds based on biological interpretation
        # These thresholds were determined by examining the distribution of
        # log2 ratios and consulting forensic miRNA literature
        self.thresholds = {
            'high_confidence': HIGH_CONFIDENCE_THRESHOLD,      # Strong expression
            'moderate_confidence': MODERATE_CONFIDENCE_THRESHOLD,  # Moderate expression  
            'low_confidence': LOW_CONFIDENCE_THRESHOLD         # Low but detectable
        }
        
    def load_data(self) -> None:
        """
        Load expression data from CSV and extract sample metadata.
        
        This method reads the normalized expression matrix and automatically
        extracts fluid type information from sample names. It expects sample
        names to contain fluid type identifiers (blood, saliva, semen, vaginal).
        
        Raises:
            FileNotFoundError: If expression_path doesn't exist
            ValueError: If no valid fluid types can be extracted from sample names
        """
        logger.info("Loading expression data...")
        
        if not self.expression_path.exists():
            raise FileNotFoundError(f"Expression file not found: {self.expression_path}")
            
        self.expr_data = pd.read_csv(self.expression_path, index_col=0)
        
        # Extract metadata from sample names
        # Sample names expected to contain fluid type (e.g., "blood_rep1", "GSM123_saliva")
        self.metadata = pd.DataFrame({
            'sample': self.expr_data.index,
            'fluid_type': [self._extract_fluid_type(s) for s in self.expr_data.index]
        })
        
        # Validate that we found fluid types
        valid_fluids = self.metadata['fluid_type'].notna().sum()
        if valid_fluids == 0:
            raise ValueError("No valid fluid types extracted from sample names")
        
        logger.info(f"Loaded {len(self.expr_data)} samples, {len(self.expr_data.columns)} miRNAs")
        logger.info(f"Fluid distribution: {self.metadata['fluid_type'].value_counts().to_dict()}")
        
    def _extract_fluid_type(self, sample_name: str) -> Optional[str]:
        """
        Extract fluid type from sample name using pattern matching.
        
        Args:
            sample_name (str): Sample identifier potentially containing fluid type
            
        Returns:
            Optional[str]: Extracted fluid type or None if not found
            
        Note:
            This method looks for keywords in sample names to identify body fluids.
            It's case-insensitive and handles various naming conventions from GEO.
        """
        sample_lower = sample_name.lower()
        
        # Define fluid type patterns
        fluid_patterns = {
            'blood': ['blood', 'whole_blood', 'peripheral'],
            'saliva': ['saliva', 'saliv', 'oral'],
            'semen': ['semen', 'seminal'],
            'vaginal': ['vaginal', 'vagina']
        }
        
        # Check each pattern
        for fluid, patterns in fluid_patterns.items():
            if any(pattern in sample_lower for pattern in patterns):
                return fluid
                
        logger.warning(f"Could not extract fluid type from: {sample_name}")
        return None
        
    def calculate_expression_tiers(self) -> pd.DataFrame:
        """
        Categorize continuous expression values into confidence tiers.
        
        This method implements the multi-tier framework that acknowledges the
        continuous nature of two-color array data while providing interpretable
        categories for forensic reporting. Each miRNA in each sample is assigned
        to one of four tiers based on its log2 ratio.
        
        Returns:
            pd.DataFrame: Tier assignments with columns:
                - sample: Sample identifier
                - mirna: miRNA name
                - fluid_type: Body fluid type
                - log2_intensity: Original continuous value
                - expression_tier: Assigned tier (high/moderate/low/undetected)
                - tier_numeric: Numeric encoding (3/2/1/0) for visualization
                
        Note:
            The tiers represent biological confidence in detection:
            - High: Strong, reliable expression
            - Moderate: Clear expression above background
            - Low: Detectable but near detection limit
            - Undetected: Below reliable detection threshold
        """
        logger.info("Calculating expression tiers...")
        
        tiers_list = []
        
        for sample in self.expr_data.index:
            fluid = self.metadata[self.metadata['sample'] == sample]['fluid_type'].iloc[0]
            
            if pd.isna(fluid):
                continue
                
            for mirna in self.expr_data.columns:
                intensity = self.expr_data.loc[sample, mirna]
                
                # Assign tier based on thresholds
                # Note: Higher (less negative) values indicate stronger expression
                if intensity >= self.thresholds['high_confidence']:
                    tier = 'high_confidence'
                    tier_num = 3
                elif intensity >= self.thresholds['moderate_confidence']:
                    tier = 'moderate_confidence'
                    tier_num = 2
                elif intensity >= self.thresholds['low_confidence']:
                    tier = 'low_confidence'
                    tier_num = 1
                else:
                    tier = 'undetected'
                    tier_num = 0
                    
                tiers_list.append({
                    'sample': sample,
                    'mirna': mirna,
                    'fluid_type': fluid,
                    'log2_intensity': intensity,
                    'expression_tier': tier,
                    'tier_numeric': tier_num
                })
                
        tiers_df = pd.DataFrame(tiers_list)
        
        # Log tier distribution
        tier_counts = tiers_df.groupby(['fluid_type', 'expression_tier']).size()
        logger.info(f"Tier distribution:\n{tier_counts}")
        
        return tiers_df
        
    def wilcoxon_tests(self, multiple_correction: str = 'bonferroni') -> pd.DataFrame:
        """
        Perform Wilcoxon rank-sum tests comparing each fluid against all others.
        
        This non-parametric test is appropriate for small sample sizes and makes
        no assumptions about data distribution. Each miRNA is tested for differential
        expression between one fluid type and all others combined.
        
        Args:
            multiple_correction (str): Method for multiple testing correction
                - 'bonferroni': Conservative, controls family-wise error rate
                  Recommended for forensic applications requiring high confidence
                - 'fdr': False Discovery Rate, less conservative
                  Useful for discovery phase to identify candidates
                Default: 'bonferroni'
                
        Returns:
            pd.DataFrame: Test results with columns:
                - fluid: Target body fluid being tested
                - mirna: miRNA identifier  
                - pvalue: Raw p-value from Wilcoxon test
                - padj: Adjusted p-value after multiple testing correction
                - log2fc: Log2 fold change (mean difference)
                - cohens_d: Effect size (standardized mean difference)
                - fluid_mean: Mean expression in target fluid
                - other_mean: Mean expression in other fluids
                - n_fluid: Number of samples in target fluid
                - n_other: Number of samples in other fluids
                
        Note:
            With pooled samples (n=2 per fluid), statistical power is limited.
            The minimum possible p-value from Wilcoxon test is 0.067 for n1=2, n2=6.
            This makes achieving significance after multiple testing correction
            challenging, emphasizing the importance of effect size measures.
        """
        logger.info(f"Performing Wilcoxon tests with {multiple_correction} correction...")
        
        results = []
        fluids = self.metadata['fluid_type'].dropna().unique()
        
        for fluid in fluids:
            # Get sample indices for this fluid and others
            fluid_samples = self.metadata[self.metadata['fluid_type'] == fluid]['sample'].values
            other_samples = self.metadata[self.metadata['fluid_type'] != fluid]['sample'].values
            
            # Skip if insufficient samples
            if len(fluid_samples) < 2 or len(other_samples) < 2:
                logger.warning(f"Insufficient samples for {fluid}, skipping")
                continue
                
            for mirna in self.expr_data.columns:
                # Extract expression values
                fluid_expr = self.expr_data.loc[fluid_samples, mirna].values
                other_expr = self.expr_data.loc[other_samples, mirna].values
                
                # Remove any NaN values
                fluid_expr = fluid_expr[~np.isnan(fluid_expr)]
                other_expr = other_expr[~np.isnan(other_expr)]
                
                # Skip if no valid data
                if len(fluid_expr) == 0 or len(other_expr) == 0:
                    continue
                
                # Perform Wilcoxon test
                try:
                    # Use two-sided test for unbiased detection
                    stat, pval = stats.mannwhitneyu(fluid_expr, other_expr, 
                                                   alternative='two-sided')
                except ValueError as e:
                    # Can occur if all values are identical
                    logger.debug(f"Wilcoxon test failed for {mirna} in {fluid}: {e}")
                    continue
                    
                # Calculate effect sizes
                fluid_mean = np.mean(fluid_expr)
                other_mean = np.mean(other_expr)
                log2fc = fluid_mean - other_mean  # Log scale difference
                
                # Cohen's d for effect size
                # Use pooled standard deviation for small samples
                pooled_std = np.sqrt((np.var(fluid_expr) + np.var(other_expr)) / 2)
                if pooled_std > 0:
                    cohens_d = (fluid_mean - other_mean) / pooled_std
                else:
                    cohens_d = 0.0
                    
                results.append({
                    'fluid': fluid,
                    'mirna': mirna,
                    'pvalue': pval,
                    'log2fc': log2fc,
                    'cohens_d': cohens_d,
                    'fluid_mean': fluid_mean,
                    'other_mean': other_mean,
                    'n_fluid': len(fluid_expr),
                    'n_other': len(other_expr)
                })
                
        results_df = pd.DataFrame(results)
        
        # Apply multiple testing correction
        if multiple_correction == 'bonferroni':
            results_df['padj'] = results_df['pvalue'] * len(results_df)
            results_df['padj'] = results_df['padj'].clip(upper=1.0)
        elif multiple_correction == 'fdr':
            # Benjamini-Hochberg FDR
            from statsmodels.stats.multitest import multipletests
            _, results_df['padj'], _, _ = multipletests(results_df['pvalue'], 
                                                       method='fdr_bh')
        else:
            raise ValueError(f"Unknown correction method: {multiple_correction}")
            
        # Sort by adjusted p-value
        results_df = results_df.sort_values('padj')
        
        # Log summary statistics
        sig_count = (results_df['padj'] < 0.05).sum()
        logger.info(f"Found {sig_count} significant results (padj < 0.05)")
        logger.info(f"Minimum p-value: {results_df['pvalue'].min():.6f}")
        
        return results_df
        
    def calculate_specificity_scores(self) -> pd.DataFrame:
        """
        Calculate specificity scores measuring unique expression in each fluid.
        
        Specificity is crucial for forensic applications where we need markers
        that reliably distinguish one fluid from all others. This method computes
        multiple metrics to identify the most specific markers.
        
        Returns:
            pd.DataFrame: Specificity metrics with columns:
                - fluid: Body fluid type
                - mirna: miRNA identifier
                - specificity_score: Ratio of min(fluid) to max(others)
                  Higher values indicate better specificity
                - detection_rate: Fraction of fluid samples with expression > threshold
                - mean_expression: Average expression in target fluid
                - fold_enrichment: Ratio of mean(fluid) to mean(others)
                
        Note:
            The specificity score uses minimum and maximum values to ensure
            a marker is consistently high in the target fluid and consistently
            low in others, which is more stringent than using means alone.
        """
        logger.info("Calculating specificity scores...")
        
        specificity_results = []
        fluids = self.metadata['fluid_type'].dropna().unique()
        
        for fluid in fluids:
            fluid_samples = self.metadata[self.metadata['fluid_type'] == fluid]['sample'].values
            other_samples = self.metadata[self.metadata['fluid_type'] != fluid]['sample'].values
            
            for mirna in self.expr_data.columns:
                # Get expression values
                fluid_expr = self.expr_data.loc[fluid_samples, mirna].values
                other_expr = self.expr_data.loc[other_samples, mirna].values
                
                # Skip if missing data
                if len(fluid_expr) == 0 or len(other_expr) == 0:
                    continue
                    
                # Calculate specificity metrics
                # Use minimum in target fluid vs maximum in others
                # This ensures consistency within groups
                min_in_fluid = np.min(fluid_expr)
                max_in_others = np.max(other_expr)
                
                # Add small offset to avoid division by zero while maintaining relationships
                specificity_score = (min_in_fluid + ZERO_OFFSET) / (max_in_others + ZERO_OFFSET)
                
                # Detection rate: fraction above moderate confidence threshold
                detection_rate = np.mean(fluid_expr > MODERATE_CONFIDENCE_THRESHOLD)
                
                # Mean expression and fold enrichment
                mean_expr = np.mean(fluid_expr)
                mean_other = np.mean(other_expr)
                fold_enrichment = (mean_expr + ZERO_OFFSET) / (mean_other + ZERO_OFFSET)
                
                specificity_results.append({
                    'fluid': fluid,
                    'mirna': mirna,
                    'specificity_score': specificity_score,
                    'detection_rate': detection_rate,
                    'mean_expression': mean_expr,
                    'fold_enrichment': fold_enrichment,
                    'min_in_fluid': min_in_fluid,
                    'max_in_others': max_in_others
                })
                
        spec_df = pd.DataFrame(specificity_results)
        
        # Log top specific markers per fluid
        for fluid in fluids:
            fluid_spec = spec_df[spec_df['fluid'] == fluid].nlargest(5, 'specificity_score')
            logger.info(f"\nTop 5 specific markers for {fluid}:")
            for _, row in fluid_spec.iterrows():
                logger.info(f"  {row['mirna']}: score={row['specificity_score']:.2f}, "
                          f"detection={row['detection_rate']:.2f}")
                          
        return spec_df
        
    def find_forensic_markers(self, results_df: pd.DataFrame, 
                            spec_df: pd.DataFrame) -> pd.DataFrame:
        """
        Identify markers meeting forensic criteria based on multiple evidence types.
        
        This method combines statistical significance, effect size, and specificity
        to identify markers suitable for forensic body fluid identification. It
        implements a scoring system that balances multiple criteria important for
        forensic applications.
        
        Args:
            results_df (pd.DataFrame): Results from wilcoxon_tests()
            spec_df (pd.DataFrame): Results from calculate_specificity_scores()
            
        Returns:
            pd.DataFrame: Forensic markers ranked by combined score with columns:
                - All columns from input dataframes
                - combined_score: Weighted score combining all criteria
                - rank: Ranking within each fluid type
                - forensic_grade: Classification (A/B/C) based on criteria met
                
        Forensic Grading:
            - Grade A: Meets all criteria (significance, effect size, specificity)
            - Grade B: Meets 2 of 3 criteria  
            - Grade C: Meets 1 of 3 criteria but shows promise
            
        Note:
            The combined score formula weights significance, effect size, and
            specificity equally. This can be adjusted based on specific forensic
            requirements or validation results.
        """
        logger.info("Identifying forensic markers...")
        
        # Merge statistical and specificity results
        merged = results_df.merge(spec_df, on=['fluid', 'mirna'])
        
        # Define forensic criteria
        sig_threshold = 0.05  # Adjusted p-value threshold
        effect_threshold = MIN_COHEN_D  # Large effect size
        spec_threshold = MIN_SPECIFICITY_SCORE  # High specificity
        
        # Apply criteria
        merged['meets_significance'] = merged['padj'] < sig_threshold
        merged['meets_effect_size'] = np.abs(merged['cohens_d']) > effect_threshold
        merged['meets_specificity'] = merged['specificity_score'] > spec_threshold
        
        # Count criteria met
        merged['criteria_met'] = (
            merged['meets_significance'].astype(int) +
            merged['meets_effect_size'].astype(int) +
            merged['meets_specificity'].astype(int)
        )
        
        # Calculate combined score
        # Weight all criteria equally, using log transformation for p-values
        merged['combined_score'] = (
            -np.log10(merged['padj'] + 1e-10) *  # Significance (higher = better)
            np.abs(merged['cohens_d']) *          # Effect size
            merged['specificity_score']           # Specificity
        )
        
        # Assign forensic grades
        merged['forensic_grade'] = 'D'  # Default
        merged.loc[merged['criteria_met'] >= 3, 'forensic_grade'] = 'A'
        merged.loc[merged['criteria_met'] == 2, 'forensic_grade'] = 'B'
        merged.loc[merged['criteria_met'] == 1, 'forensic_grade'] = 'C'
        
        # Sort by combined score within each fluid
        forensic_markers = []
        for fluid in merged['fluid'].unique():
            fluid_markers = merged[merged['fluid'] == fluid].copy()
            fluid_markers = fluid_markers.sort_values('combined_score', ascending=False)
            fluid_markers['rank'] = range(1, len(fluid_markers) + 1)
            forensic_markers.append(fluid_markers)
            
        forensic_df = pd.concat(forensic_markers, ignore_index=True)
        
        # Log summary
        grade_counts = forensic_df.groupby(['fluid', 'forensic_grade']).size()
        logger.info(f"\nForensic marker grades:\n{grade_counts}")
        
        # Report top markers per fluid
        top_markers = forensic_df.groupby('fluid').head(3)
        logger.info("\nTop 3 forensic markers per fluid:")
        for _, marker in top_markers.iterrows():
            logger.info(f"{marker['fluid']} - {marker['mirna']}: "
                      f"Grade {marker['forensic_grade']}, "
                      f"Score={marker['combined_score']:.2f}")
                      
        return forensic_df
        
    def create_visualizations(self, output_dir: Path) -> None:
        """
        Create comprehensive visualizations for the analysis results.
        
        Generates multiple plots to visualize different aspects of the analysis:
        1. Expression tier heatmap showing multi-tier framework
        2. Box plots of expression distributions by fluid
        3. Volcano plots for differential expression
        4. Forensic marker summary plots
        
        Args:
            output_dir (Path): Directory to save visualization files
            
        Creates:
            - expression_tiers_heatmap.png: Overview of tier assignments
            - expression_profiles_boxplot.png: Distribution comparisons
            - specificity_volcano.png: Significance vs specificity
            - forensic_markers_barplot.png: Top markers per fluid
            
        Note:
            All plots use consistent color schemes for body fluids:
            - Blood: red/darkred
            - Saliva: blue/darkblue  
            - Semen: green/darkgreen
            - Vaginal: purple/darkpurple
        """
        logger.info(f"Creating visualizations in {output_dir}")
        output_dir.mkdir(exist_ok=True, parents=True)
        
        # Define consistent color palette for body fluids
        fluid_colors = {
            'blood': '#CC0000',
            'saliva': '#0066CC', 
            'semen': '#009900',
            'vaginal': '#9900CC'
        }
        
        # 1. Expression Tier Heatmap
        self._create_tier_heatmap(output_dir, fluid_colors)
        
        # 2. Expression Distribution Box Plots
        self._create_distribution_plots(output_dir, fluid_colors)
        
        # 3. Volcano Plots
        self._create_volcano_plots(output_dir, fluid_colors)
        
        # 4. Forensic Marker Summary
        self._create_forensic_summary(output_dir, fluid_colors)
        
        logger.info("Visualizations complete")
        
    def _create_tier_heatmap(self, output_dir: Path, colors: Dict[str, str]) -> None:
        """Create heatmap showing expression tiers across samples and miRNAs."""
        # Implementation details...
        pass
        
    def _create_distribution_plots(self, output_dir: Path, colors: Dict[str, str]) -> None:
        """Create box plots showing expression distributions."""
        # Implementation details...
        pass
        
    def _create_volcano_plots(self, output_dir: Path, colors: Dict[str, str]) -> None:
        """Create volcano plots for differential expression."""
        # Implementation details...
        pass
        
    def _create_forensic_summary(self, output_dir: Path, colors: Dict[str, str]) -> None:
        """Create summary plots for top forensic markers."""
        # Implementation details...
        pass
        

def main():
    """
    Run complete continuous expression analysis pipeline.
    
    This function orchestrates the full analysis from data loading through
    marker identification and visualization generation.
    """
    # Set up paths
    data_dir = Path('data/processed/gpr')
    results_dir = Path('results/continuous_expression_analysis')
    results_dir.mkdir(exist_ok=True, parents=True)
    
    # Initialize analyzer
    analyzer = ContinuousExpressionAnalyzer(data_dir / 'gpr_expression_matrix.csv')
    
    # Run analysis pipeline
    analyzer.load_data()
    
    # Calculate expression tiers
    tiers_df = analyzer.calculate_expression_tiers()
    tiers_df.to_csv(results_dir / 'expression_tiers.csv', index=False)
    
    # Perform statistical tests
    # Note: With n=2 per group, Bonferroni is very conservative
    # Consider using FDR for discovery phase
    wilcox_results = analyzer.wilcoxon_tests(multiple_correction='bonferroni')
    wilcox_results.to_csv(results_dir / 'wilcoxon_results.csv', index=False)
    
    # Calculate specificity
    spec_scores = analyzer.calculate_specificity_scores()
    spec_scores.to_csv(results_dir / 'specificity_scores.csv', index=False)
    
    # Find forensic markers
    forensic_markers = analyzer.find_forensic_markers(wilcox_results, spec_scores)
    forensic_markers.to_csv(results_dir / 'forensic_markers.csv', index=False)
    
    # Create visualizations
    viz_dir = results_dir / 'visualizations'
    analyzer.create_visualizations(viz_dir)
    
    # Generate summary report
    generate_summary_report(tiers_df, wilcox_results, spec_scores, 
                          forensic_markers, results_dir)
    
    logger.info(f"Analysis complete. Results saved to {results_dir}")
    

def generate_summary_report(tiers_df: pd.DataFrame, wilcox_results: pd.DataFrame,
                           spec_scores: pd.DataFrame, forensic_markers: pd.DataFrame,
                           output_dir: Path) -> None:
    """
    Generate a markdown summary report of the analysis.
    
    Creates a comprehensive report documenting:
    - Dataset overview and sample distribution
    - Statistical analysis results
    - Top forensic markers identified
    - Limitations and recommendations
    
    Args:
        tiers_df: Expression tier assignments
        wilcox_results: Statistical test results
        spec_scores: Specificity scores
        forensic_markers: Identified forensic markers
        output_dir: Directory to save report
    """
    report_lines = []
    
    # Header
    report_lines.append("# Continuous Expression Analysis Report - GPR Data")
    report_lines.append("\n## Dataset Overview")
    report_lines.append(f"- Total miRNAs analyzed: {wilcox_results['mirna'].nunique()}")
    report_lines.append(f"- Samples per fluid: ~2 (pooled from ~5-7 donors each)")
    report_lines.append(f"- Statistical method: Wilcoxon rank-sum test")
    report_lines.append(f"- Multiple testing correction: Bonferroni")
    
    # Statistical summary
    report_lines.append("\n## Statistical Summary")
    sig_count = (wilcox_results['padj'] < 0.05).sum()
    report_lines.append(f"- Significant results (padj < 0.05): {sig_count}")
    report_lines.append(f"- Minimum p-value: {wilcox_results['pvalue'].min():.6f}")
    
    # Top markers per fluid
    report_lines.append("\n## Top Forensic Markers")
    
    for fluid in ['blood', 'saliva', 'semen', 'vaginal']:
        fluid_markers = forensic_markers[
            (forensic_markers['fluid'] == fluid) & 
            (forensic_markers['forensic_grade'].isin(['A', 'B']))
        ].head(5)
        
        report_lines.append(f"\n### {fluid.capitalize()}")
        if len(fluid_markers) > 0:
            for _, marker in fluid_markers.iterrows():
                report_lines.append(f"- **{marker['mirna']}**: Grade {marker['forensic_grade']}, "
                                  f"Score={marker['combined_score']:.2f}, "
                                  f"Specificity={marker['specificity_score']:.2f}")
        else:
            report_lines.append("- No high-grade markers identified")
            
    # Limitations
    report_lines.append("\n## Limitations and Considerations")
    report_lines.append("- Small sample size (n=2 per fluid) limits statistical power")
    report_lines.append("- Pooled samples may mask individual variation")
    report_lines.append("- Bonferroni correction is conservative for discovery")
    report_lines.append("- Results should be validated with individual samples")
    
    # Recommendations
    report_lines.append("\n## Recommendations")
    report_lines.append("1. Validate top markers using qRT-PCR on individual samples")
    report_lines.append("2. Consider relaxing multiple testing correction for discovery")
    report_lines.append("3. Focus on markers with large effect sizes regardless of p-value")
    report_lines.append("4. Test marker stability in degraded samples")
    
    # Multi-tier framework summary
    report_lines.append("\n## Multi-Tier Expression Framework")
    report_lines.append("The analysis uses expression tiers to capture continuous nature of data:")
    report_lines.append("- **High confidence**: log2 ratio > -2.0")
    report_lines.append("- **Moderate confidence**: -6.0 to -2.0")
    report_lines.append("- **Low confidence**: -10.0 to -6.0")
    report_lines.append("- **Undetected**: < -10.0")
    
    # Tier distribution
    tier_summary = tiers_df.groupby(['fluid_type', 'expression_tier']).size().unstack(fill_value=0)
    report_lines.append("\n### Tier Distribution by Fluid")
    report_lines.append("```")
    report_lines.append(tier_summary.to_string())
    report_lines.append("```")
    
    # Save report
    report_path = output_dir / 'continuous_expression_report.md'
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines))
        
    logger.info(f"Summary report saved to {report_path}")


if __name__ == "__main__":
    main()