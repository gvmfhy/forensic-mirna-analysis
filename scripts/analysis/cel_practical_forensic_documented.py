#!/usr/bin/env python3
"""
Practical Forensic Analysis of CEL Data

This module implements a pragmatic approach to forensic miRNA marker identification
from Affymetrix CEL files (GSE49630). It acknowledges the statistical limitations
of small sample sizes (n=5 per fluid) and focuses on identifying markers with
large effect sizes and biological relevance rather than strict statistical significance.

Key Features:
    - Relaxed statistical criteria appropriate for small samples
    - Focus on effect size (Cohen's d) over p-values
    - Practical forensic thresholds based on fold change
    - Comprehensive visualization of marker characteristics

The analysis philosophy prioritizes biological and forensic relevance over
strict adherence to traditional statistical thresholds, recognizing that
with n=5 samples, achieving FDR < 0.05 across 13,564 tests is mathematically
impossible with rank-based tests.

Author: Forensic miRNA Analysis Pipeline
Date: 2024
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
import logging
from typing import Dict, List, Tuple, Optional, Union

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Forensic Analysis Constants
P_VALUE_THRESHOLD = 0.01  # Nominal p-value threshold for initial filtering
FOLD_CHANGE_THRESHOLD = 2  # Log2 fold change threshold (equals 4x linear fold change)
EFFECT_SIZE_THRESHOLD = 1.5  # Cohen's d threshold for very large effect size
DETECTION_RATE_THRESHOLD = 1.0  # Require detection in 100% of target fluid samples
MIN_EXPRESSION_DIFF = 3.0  # Minimum difference in mean expression for biological relevance

# Visualization Constants
HEATMAP_TOP_N = 3  # Number of top markers per fluid to show in heatmap
EFFECT_PLOT_TOP_N = 5  # Number of top markers per fluid for detailed plots
DPI = 300  # Resolution for saved figures
FIGURE_FORMAT = 'png'  # Output format for figures

# Body Fluid Color Scheme (consistent across all visualizations)
FLUID_COLORS = {
    'blood': '#8B0000',      # Dark red
    'saliva': '#4169E1',     # Royal blue
    'semen': '#228B22',      # Forest green
    'vaginal': '#8B008B'     # Dark magenta
}


def practical_forensic_analysis(results_dir: Path = Path('results/cel_continuous_expression')) -> Dict[str, pd.DataFrame]:
    """
    Perform practical forensic analysis on CEL data with relaxed statistical criteria.
    
    This function implements a forensic-focused analysis that acknowledges the
    limitations of small sample sizes while identifying biologically meaningful
    markers for body fluid identification. It prioritizes large effect sizes
    and consistent detection over strict statistical significance.
    
    The analysis uses the following criteria hierarchy:
    1. Nominal p-value < 0.01 (uncorrected) as initial filter
    2. Absolute log2 fold change > 2 (4-fold linear change)
    3. Cohen's d > 1.5 (very large effect size)
    4. 100% detection rate in target fluid
    
    Args:
        results_dir (Path): Directory containing preprocessed Wilcoxon test results
            and specificity scores from previous analysis steps.
            Expected files:
            - wilcoxon_results.csv: Statistical test results
            - specificity_scores.csv: Fluid-specific expression metrics
            
    Returns:
        Dict[str, pd.DataFrame]: Dictionary mapping body fluid names to DataFrames
            containing their respective forensic marker candidates. Keys are:
            'blood', 'saliva', 'semen', 'vaginal'. Each DataFrame includes:
            - mirna: miRNA identifier
            - pvalue, fdr: Statistical significance metrics
            - log2fc: Log2 fold change
            - cohens_d: Effect size
            - score: Combined importance score
            - detection_rate: Fraction of samples with expression
            
    Raises:
        FileNotFoundError: If required input files are not found
        ValueError: If data format is unexpected
        
    Note:
        The function generates a comprehensive report explaining why strict
        FDR control is not achievable with small samples and provides
        practical recommendations for forensic implementation.
    """
    logger.info("Starting practical forensic analysis of CEL data")
    
    # Validate input directory
    if not results_dir.exists():
        raise FileNotFoundError(f"Results directory not found: {results_dir}")
    
    # Load preprocessed results
    try:
        wilcox_df = pd.read_csv(results_dir / 'wilcoxon_results.csv')
        spec_df = pd.read_csv(results_dir / 'specificity_scores.csv')
    except FileNotFoundError as e:
        logger.error(f"Required input file not found: {e}")
        raise
    
    # Filter for human miRNAs only (exclude controls and non-human sequences)
    # Human miRNAs are identified by 'hsa-' prefix
    human_mask = wilcox_df['mirna'].str.contains('hsa-', na=False)
    wilcox_human = wilcox_df[human_mask].copy()
    spec_human = spec_df[spec_df['mirna'].str.contains('hsa-', na=False)].copy()
    
    logger.info(f"Analyzing {len(wilcox_human)} human miRNA comparisons")
    
    # Apply FDR correction for completeness (though we expect few to pass)
    # Using Benjamini-Hochberg method as it's less conservative than Bonferroni
    _, wilcox_human['fdr'], _, _ = multipletests(wilcox_human['pvalue'], method='fdr_bh')
    
    # Log FDR results to demonstrate the challenge
    fdr_significant = (wilcox_human['fdr'] < 0.05).sum()
    min_fdr = wilcox_human['fdr'].min()
    logger.info(f"Markers passing FDR < 0.05: {fdr_significant}")
    logger.info(f"Minimum FDR q-value: {min_fdr:.3f}")
    
    # Merge statistical results with specificity scores
    merged = wilcox_human.merge(spec_human, on=['fluid', 'mirna'], how='inner')
    
    # Create output directory for results
    output_dir = Path('results/cel_practical_forensic')
    output_dir.mkdir(exist_ok=True)
    
    # Generate comprehensive report
    report = generate_forensic_report(merged, output_dir)
    
    # Find best markers per fluid with practical criteria
    markers_by_fluid = {}
    
    for fluid in ['blood', 'saliva', 'semen', 'vaginal']:
        # Get fluid-specific data
        fluid_data = merged[merged['fluid'] == fluid].copy()
        
        # Apply practical forensic criteria
        # These thresholds are based on forensic literature and validation studies
        candidates = fluid_data[
            (fluid_data['pvalue'] < P_VALUE_THRESHOLD) &  # Nominal significance
            (np.abs(fluid_data['log2fc']) > FOLD_CHANGE_THRESHOLD) &  # Large fold change
            (np.abs(fluid_data['cohens_d']) > EFFECT_SIZE_THRESHOLD) &  # Very large effect
            (fluid_data['detection_rate'] >= DETECTION_RATE_THRESHOLD)  # Always detected
        ].copy()
        
        # Calculate combined importance score
        # This scoring system weights three factors equally:
        # 1. Statistical significance (via -log10 transformation)
        # 2. Fold change magnitude (biological relevance)
        # 3. Effect size (practical discrimination ability)
        candidates['score'] = (
            -np.log10(candidates['pvalue'] + 1e-10) *  # Avoid log(0)
            np.abs(candidates['log2fc']) * 
            np.abs(candidates['cohens_d'])
        )
        
        # Sort by combined score
        candidates = candidates.sort_values('score', ascending=False)
        
        # Store results
        markers_by_fluid[fluid] = candidates
        
        # Log summary for this fluid
        logger.info(f"{fluid.capitalize()}: {len(candidates)} markers meeting all criteria")
        if len(candidates) > 0:
            top_marker = candidates.iloc[0]
            logger.info(f"  Top marker: {top_marker['mirna']} "
                       f"(FC={2**top_marker['log2fc']:.1f}x, d={top_marker['cohens_d']:.2f})")
    
    # Create comprehensive visualizations
    create_practical_visualizations(markers_by_fluid, output_dir)
    
    # Save detailed results for all candidates
    all_candidates = pd.concat(markers_by_fluid.values(), ignore_index=True)
    all_candidates.to_csv(output_dir / 'forensic_candidates.csv', index=False)
    
    logger.info(f"Analysis complete. Found {len(all_candidates)} total candidates")
    logger.info(f"Results saved to {output_dir}")
    
    return markers_by_fluid


def generate_forensic_report(merged_data: pd.DataFrame, output_dir: Path) -> List[str]:
    """
    Generate a comprehensive forensic analysis report.
    
    Creates a detailed markdown report explaining the analysis approach,
    statistical considerations, and practical recommendations for forensic
    implementation. The report is designed for both computational biologists
    and forensic practitioners.
    
    Args:
        merged_data (pd.DataFrame): Combined statistical and specificity results
        output_dir (Path): Directory to save the report
        
    Returns:
        List[str]: Report content as list of lines (also saved to file)
        
    Note:
        The report includes a detailed explanation of why FDR < 0.05 is not
        achievable with n=5 samples and 13,564 tests, helping users understand
        the statistical limitations and the rationale for relaxed criteria.
    """
    report = []
    report.append("# Practical Forensic miRNA Analysis - CEL Data")
    report.append("\n## Executive Summary")
    report.append("With n=5 individual samples per body fluid, we can identify")
    report.append("scientifically meaningful markers for forensic identification.")
    
    # Analysis approach section
    report.append("\n## Analysis Approach")
    report.append("Given the mathematical limitations of achieving strict FDR control")
    report.append("with small samples, this analysis focuses on:")
    report.append("1. Large effect sizes (Cohen's d > 1.5)")
    report.append("2. Substantial fold changes (>4-fold)")
    report.append("3. Consistent detection (100% of samples)")
    report.append("4. Biological plausibility")
    
    # Process each fluid type
    markers_by_fluid = {}
    for fluid in ['blood', 'saliva', 'semen', 'vaginal']:
        fluid_data = merged_data[merged_data['fluid'] == fluid].copy()
        
        # Apply criteria
        candidates = fluid_data[
            (fluid_data['pvalue'] < P_VALUE_THRESHOLD) &
            (np.abs(fluid_data['log2fc']) > FOLD_CHANGE_THRESHOLD) &
            (np.abs(fluid_data['cohens_d']) > EFFECT_SIZE_THRESHOLD) &
            (fluid_data['detection_rate'] >= DETECTION_RATE_THRESHOLD)
        ].copy()
        
        # Calculate scores
        candidates['score'] = (
            -np.log10(candidates['pvalue'] + 1e-10) *
            np.abs(candidates['log2fc']) *
            np.abs(candidates['cohens_d'])
        )
        candidates = candidates.sort_values('score', ascending=False)
        
        markers_by_fluid[fluid] = candidates
        
        # Add to report
        report.append(f"\n## {fluid.capitalize()}")
        report.append(f"Found {len(candidates)} strong candidates (p<0.01, |FC|>4, |d|>1.5)")
        
        if len(candidates) > 0:
            report.append("\n### Top 5 Markers:")
            for i, (_, marker) in enumerate(candidates.head().iterrows()):
                report.append(f"\n{i+1}. **{marker['mirna'].replace('_st', '')}**")
                
                # Calculate linear fold change with direction
                fc = 2**marker['log2fc']
                direction = 'higher' if marker['log2fc'] > 0 else 'lower'
                report.append(f"   - Fold change: {fc:.1f}x ({direction} in {fluid})")
                
                report.append(f"   - Effect size (Cohen's d): {marker['cohens_d']:.2f}")
                report.append(f"   - p-value: {marker['pvalue']:.4f}")
                report.append(f"   - FDR q-value: {marker['fdr']:.3f}")
                report.append(f"   - Expression: {marker['fluid_mean']:.1f} vs {marker['other_mean']:.1f}")
    
    # Statistical considerations section
    report.append("\n## Statistical Considerations")
    report.append("\n### Multiple Testing Challenge")
    
    # Calculate exact numbers
    total_tests = len(merged_data)
    human_mirnas = merged_data['mirna'].str.contains('hsa-').sum()
    min_pval = merged_data['pvalue'].min()
    
    report.append(f"- Total human miRNA tests: {human_mirnas:,}")
    report.append(f"- Minimum p-value: {min_pval:.6f}")
    report.append(f"- Tests with p < 0.01: {(merged_data['pvalue'] < 0.01).sum()}")
    report.append(f"- Minimum FDR q-value: {merged_data['fdr'].min():.3f}")
    
    report.append("\n### Why Strict FDR < 0.05 Not Achieved")
    report.append("- With rank-based tests and n=5, minimum possible p-value is ~0.001")
    report.append("- Testing 13,564 miRNAs requires p < 0.0000037 for FDR < 0.05")
    report.append("- This is mathematically impossible with n=5")
    
    report.append("\n### Practical Approach")
    report.append("- Focus on effect size (fold change) and consistency")
    report.append("- Validate findings across platforms")
    report.append("- Use targeted assays for final forensic application")
    
    # Forensic implementation section
    report.append("\n## Forensic Implementation")
    report.append("\n### Recommended Marker Panel")
    
    # Select top 2 markers per fluid with highest fold changes
    panel = []
    for fluid, candidates in markers_by_fluid.items():
        if len(candidates) > 0:
            # Get markers with positive fold change (higher in target)
            positive_fc = candidates[candidates['log2fc'] > 0].nlargest(2, 'log2fc')
            for _, marker in positive_fc.iterrows():
                panel.append({
                    'fluid': fluid,
                    'mirna': marker['mirna'].replace('_st', ''),
                    'fold_change': 2**marker['log2fc'],
                    'cohens_d': marker['cohens_d']
                })
    
    # Format panel recommendations
    panel_df = pd.DataFrame(panel)
    if len(panel_df) > 0:
        for fluid in ['blood', 'saliva', 'semen', 'vaginal']:
            fluid_panel = panel_df[panel_df['fluid'] == fluid]
            if len(fluid_panel) > 0:
                report.append(f"\n**{fluid.capitalize()}:**")
                for _, m in fluid_panel.iterrows():
                    report.append(f"- {m['mirna']}: {m['fold_change']:.0f}-fold enriched")
    
    report.append("\n### Validation Strategy")
    report.append("1. Confirm expression patterns with qRT-PCR")
    report.append("2. Test on independent sample cohort")
    report.append("3. Assess stability in degraded samples")
    report.append("4. Develop multiplex assay for simultaneous detection")
    
    # Save report
    report_path = output_dir / 'practical_forensic_report.md'
    with open(report_path, 'w') as f:
        f.write('\n'.join(report))
    
    logger.info(f"Report saved to {report_path}")
    
    return report


def create_practical_visualizations(markers_by_fluid: Dict[str, pd.DataFrame], 
                                  output_dir: Path) -> None:
    """
    Create comprehensive visualizations for practical forensic analysis.
    
    Generates multiple publication-quality figures that illustrate the
    forensic marker characteristics from different perspectives. All plots
    use consistent styling and color schemes for professional presentation.
    
    Args:
        markers_by_fluid (Dict[str, pd.DataFrame]): Dictionary mapping fluid names
            to DataFrames of marker candidates
        output_dir (Path): Directory to save visualization files
        
    Creates:
        - top_markers_heatmap.png: Expression patterns of top markers per fluid
        - effect_sizes.png: Scatter plot of fold change vs Cohen's d
        - expression_comparison.png: Bar plots comparing expression levels
        
    Note:
        All visualizations use the FLUID_COLORS constant for consistent
        color coding across figures, making them suitable for publication
        or presentation as a coherent set.
    """
    logger.info("Creating practical forensic visualizations")
    
    # Set style for all plots
    sns.set_style("whitegrid")
    plt.rcParams['figure.dpi'] = DPI
    plt.rcParams['savefig.dpi'] = DPI
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['axes.titlesize'] = 14
    
    # 1. Create top markers heatmap
    create_top_markers_heatmap(markers_by_fluid, output_dir)
    
    # 2. Create effect size scatter plot
    create_effect_size_plot(markers_by_fluid, output_dir)
    
    # 3. Create expression level comparison
    create_expression_comparison(markers_by_fluid, output_dir)
    
    logger.info("Visualizations saved to output directory")


def create_top_markers_heatmap(markers_by_fluid: Dict[str, pd.DataFrame], 
                              output_dir: Path) -> None:
    """
    Create heatmap showing expression patterns of top forensic markers.
    
    This visualization shows the top markers for each body fluid in a heatmap
    format, making it easy to see fluid-specific expression patterns. The
    heatmap uses a diverging color scheme centered at zero to highlight
    both up- and down-regulated markers.
    
    Args:
        markers_by_fluid: Dictionary of marker DataFrames by fluid
        output_dir: Directory to save the figure
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Collect top markers from each fluid
    top_markers = []
    for fluid, candidates in markers_by_fluid.items():
        if len(candidates) > 0:
            # Get top N markers based on combined score
            top = candidates.nlargest(HEATMAP_TOP_N, 'score')[['mirna', 'fluid', 'log2fc']]
            # Clean miRNA names (remove _st suffix if present)
            top['mirna'] = top['mirna'].str.replace('_st', '')
            top_markers.append(top)
    
    if len(top_markers) > 0:
        # Combine all top markers
        all_top = pd.concat(top_markers)
        
        # Pivot for heatmap format
        pivot = all_top.pivot(index='mirna', columns='fluid', values='log2fc')
        
        # Fill NaN with 0 for visualization
        pivot = pivot.fillna(0)
        
        # Create heatmap with custom colormap
        sns.heatmap(pivot, cmap='RdBu_r', center=0,
                   annot=True, fmt='.1f',
                   cbar_kws={'label': 'Log2 Fold Change'},
                   ax=ax, vmin=-8, vmax=8,
                   linewidths=0.5, linecolor='gray')
        
        # Customize appearance
        ax.set_title('Top Forensic miRNA Candidates\n(p<0.01, |FC|>4, |d|>1.5)', 
                    fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel('Body Fluid', fontsize=12)
        ax.set_ylabel('miRNA', fontsize=12)
        
        # Rotate labels for readability
        plt.setp(ax.get_xticklabels(), rotation=0, ha='center')
        plt.setp(ax.get_yticklabels(), rotation=0)
        
        # Add grid for clarity
        ax.grid(False)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'top_markers_heatmap.png', 
                   dpi=DPI, bbox_inches='tight')
        plt.close()
        
        logger.info("Created top markers heatmap")


def create_effect_size_plot(markers_by_fluid: Dict[str, pd.DataFrame], 
                           output_dir: Path) -> None:
    """
    Create scatter plot showing relationship between fold change and effect size.
    
    This plot helps visualize the relationship between biological significance
    (fold change) and statistical effect size (Cohen's d), with points colored
    by body fluid type. Threshold lines indicate the forensic criteria.
    
    Args:
        markers_by_fluid: Dictionary of marker DataFrames by fluid
        output_dir: Directory to save the figure
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Combine all candidates
    all_candidates = pd.concat(markers_by_fluid.values(), ignore_index=True)
    
    if len(all_candidates) > 0:
        # Create scatter plot for each fluid
        for fluid in ['blood', 'saliva', 'semen', 'vaginal']:
            fluid_data = all_candidates[all_candidates['fluid'] == fluid]
            if len(fluid_data) > 0:
                ax.scatter(fluid_data['log2fc'], fluid_data['cohens_d'],
                          color=FLUID_COLORS[fluid], label=fluid.capitalize(),
                          alpha=0.6, s=100, edgecolor='black', linewidth=0.5)
        
        # Add threshold lines
        ax.axhline(EFFECT_SIZE_THRESHOLD, color='gray', linestyle='--', 
                  alpha=0.5, label=f'd = {EFFECT_SIZE_THRESHOLD}')
        ax.axhline(-EFFECT_SIZE_THRESHOLD, color='gray', linestyle='--', alpha=0.5)
        ax.axvline(FOLD_CHANGE_THRESHOLD, color='gray', linestyle='--', 
                  alpha=0.5, label=f'FC = {2**FOLD_CHANGE_THRESHOLD}x')
        ax.axvline(-FOLD_CHANGE_THRESHOLD, color='gray', linestyle='--', alpha=0.5)
        
        # Shade the regions meeting criteria
        ax.fill_between([FOLD_CHANGE_THRESHOLD, ax.get_xlim()[1]], 
                       EFFECT_SIZE_THRESHOLD, ax.get_ylim()[1],
                       color='green', alpha=0.1, label='Meets criteria')
        ax.fill_between([ax.get_xlim()[0], -FOLD_CHANGE_THRESHOLD], 
                       ax.get_ylim()[0], -EFFECT_SIZE_THRESHOLD,
                       color='green', alpha=0.1)
        
        # Labels and formatting
        ax.set_xlabel('Log2 Fold Change', fontsize=12)
        ax.set_ylabel("Cohen's d", fontsize=12)
        ax.set_title('Effect Sizes of Forensic Candidates', 
                    fontsize=14, fontweight='bold')
        ax.legend(loc='upper left', frameon=True, fancybox=True, shadow=True)
        ax.grid(True, alpha=0.3)
        
        # Add text annotations for quadrants
        ax.text(5, 10, 'Strong\nUp-regulation', ha='center', va='center',
               fontsize=10, style='italic', alpha=0.7)
        ax.text(-5, -10, 'Strong\nDown-regulation', ha='center', va='center',
               fontsize=10, style='italic', alpha=0.7)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'effect_sizes.png', 
                   dpi=DPI, bbox_inches='tight')
        plt.close()
        
        logger.info("Created effect size plot")


def create_expression_comparison(markers_by_fluid: Dict[str, pd.DataFrame], 
                               output_dir: Path) -> None:
    """
    Create grouped bar plots comparing expression levels between fluids.
    
    This visualization shows the actual expression levels of top markers in
    their target fluid versus other fluids, providing a clear view of the
    discriminatory power of each marker.
    
    Args:
        markers_by_fluid: Dictionary of marker DataFrames by fluid
        output_dir: Directory to save the figure
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    for i, (fluid, candidates) in enumerate(markers_by_fluid.items()):
        if i < 4 and len(candidates) > 0:
            ax = axes[i]
            
            # Get top markers by score
            top5 = candidates.nlargest(EFFECT_PLOT_TOP_N, 'score')
            
            if len(top5) > 0:
                # Prepare data for grouped bars
                x = np.arange(len(top5))
                width = 0.35
                
                # Create bars
                bars1 = ax.bar(x - width/2, top5['fluid_mean'], width, 
                              label=f'{fluid.capitalize()}', 
                              color=FLUID_COLORS[fluid], alpha=0.8)
                bars2 = ax.bar(x + width/2, top5['other_mean'], width,
                              label='Other fluids', color='gray', alpha=0.5)
                
                # Customize axes
                ax.set_ylabel('Expression Level', fontsize=10)
                ax.set_title(f'{fluid.capitalize()} Top Markers', 
                           fontsize=12, fontweight='bold')
                ax.set_xticks(x)
                
                # Clean miRNA names for display
                mirna_labels = top5['mirna'].str.replace('_st', '').values
                ax.set_xticklabels(mirna_labels, rotation=45, ha='right', fontsize=8)
                ax.legend(loc='upper right', fontsize=8)
                
                # Add fold change annotations on top of bars
                for j, (_, row) in enumerate(top5.iterrows()):
                    fc = 2**row['log2fc']
                    y_pos = max(row['fluid_mean'], row['other_mean']) + 0.5
                    ax.text(j, y_pos, f'{fc:.0f}x', ha='center', 
                           fontsize=8, fontweight='bold')
                
                # Add grid for readability
                ax.grid(True, axis='y', alpha=0.3)
                ax.set_axisbelow(True)
    
    # Add overall title
    plt.suptitle('Expression Levels: Target vs Other Fluids', 
                fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'expression_comparison.png', 
               dpi=DPI, bbox_inches='tight')
    plt.close()
    
    logger.info("Created expression comparison plot")


def main():
    """
    Execute the practical forensic analysis pipeline.
    
    This is the main entry point that runs the complete analysis from
    loading preprocessed data through generating reports and visualizations.
    """
    try:
        # Run analysis
        markers_by_fluid = practical_forensic_analysis()
        
        # Summary statistics
        total = sum(len(m) for m in markers_by_fluid.values())
        logger.info(f"\nSummary:")
        logger.info(f"Total candidates across all fluids: {total}")
        
        for fluid, markers in markers_by_fluid.items():
            if len(markers) > 0:
                best_fc = 2**markers.iloc[0]['log2fc']
                logger.info(f"{fluid}: {len(markers)} markers (best FC: {best_fc:.0f}x)")
                
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise


if __name__ == "__main__":
    main()