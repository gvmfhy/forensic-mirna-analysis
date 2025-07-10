#!/usr/bin/env python3
"""
Forensic Signature Analysis - Identifying body fluid-specific miRNA patterns
This is what we SHOULD have done from the start
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

def load_gpr_data():
    """Load the GPR expression data"""
    expr = pd.read_csv('data/processed/gpr/gpr_expression_matrix.csv', index_col=0)
    raw = pd.read_csv('data/processed/gpr/gpr_combined_raw.csv')
    
    # Get metadata
    metadata = raw[['sample', 'body_fluid']].drop_duplicates().set_index('sample')
    
    return expr, metadata

def calculate_fluid_signatures(expr, metadata):
    """Calculate mean expression and specificity for each miRNA in each fluid"""
    results = []
    
    fluids = metadata['body_fluid'].unique()
    
    for mirna in expr.columns:
        mirna_data = {'miRNA': mirna}
        
        # Calculate mean expression per fluid
        for fluid in fluids:
            fluid_samples = metadata[metadata['body_fluid'] == fluid].index
            fluid_expr = expr.loc[fluid_samples, mirna]
            
            mirna_data[f'{fluid}_mean'] = fluid_expr.mean()
            mirna_data[f'{fluid}_std'] = fluid_expr.std()
            mirna_data[f'{fluid}_n'] = len(fluid_expr)
        
        # Calculate specificity score (how specific is this miRNA to each fluid)
        for target_fluid in fluids:
            target_mean = mirna_data[f'{target_fluid}_mean']
            other_means = [mirna_data[f'{f}_mean'] for f in fluids if f != target_fluid]
            
            # Specificity score: target expression / mean of others
            if np.mean(other_means) != 0:
                specificity = target_mean / np.mean(other_means)
            else:
                specificity = np.inf if target_mean > 0 else 1
            
            mirna_data[f'{target_fluid}_specificity'] = specificity
        
        results.append(mirna_data)
    
    return pd.DataFrame(results)

def find_fluid_specific_markers(signatures, threshold=2.0):
    """Find miRNAs that are highly specific to each fluid"""
    fluids = ['menstrual_blood', 'peripheral_blood', 'saliva', 'semen', 'vaginal_secretion']
    
    fluid_markers = {}
    
    for fluid in fluids:
        # Find miRNAs with high specificity for this fluid
        fluid_specific = signatures[
            (signatures[f'{fluid}_specificity'] > threshold) & 
            (signatures[f'{fluid}_mean'] > signatures[[f'{f}_mean' for f in fluids]].mean(axis=1))
        ].copy()
        
        fluid_specific = fluid_specific.sort_values(f'{fluid}_specificity', ascending=False)
        fluid_markers[fluid] = fluid_specific.head(10)['miRNA'].tolist()
    
    return fluid_markers

def create_forensic_visualizations(expr, metadata, signatures):
    """Create proper forensic-focused visualizations"""
    
    fig = plt.figure(figsize=(20, 15))
    
    # 1. Expression distributions by fluid
    ax1 = plt.subplot(3, 2, 1)
    
    # Prepare data for box plot
    plot_data = []
    for fluid in metadata['body_fluid'].unique():
        fluid_samples = metadata[metadata['body_fluid'] == fluid].index
        for sample in fluid_samples:
            for mirna in expr.columns[:20]:  # Top 20 for visibility
                plot_data.append({
                    'Body Fluid': fluid,
                    'Expression': expr.loc[sample, mirna],
                    'miRNA': mirna
                })
    
    plot_df = pd.DataFrame(plot_data)
    sns.boxplot(data=plot_df, x='Body Fluid', y='Expression', ax=ax1)
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')
    ax1.set_title('Expression Distribution by Body Fluid')
    ax1.set_ylabel('Log2 Expression')
    
    # 2. Specificity scores heatmap
    ax2 = plt.subplot(3, 2, 2)
    
    # Get top specific miRNAs for each fluid
    fluids = ['menstrual_blood', 'peripheral_blood', 'saliva', 'semen', 'vaginal_secretion']
    top_specific = []
    
    for fluid in fluids:
        top = signatures.nlargest(5, f'{fluid}_specificity')
        top_specific.extend(top['miRNA'].tolist())
    
    top_specific = list(set(top_specific))[:25]  # Unique, limit to 25
    
    # Create specificity matrix
    spec_matrix = []
    for mirna in top_specific:
        row = signatures[signatures['miRNA'] == mirna].iloc[0]
        spec_matrix.append([row[f'{fluid}_specificity'] for fluid in fluids])
    
    spec_df = pd.DataFrame(spec_matrix, index=top_specific, columns=fluids)
    
    sns.heatmap(spec_df, cmap='YlOrRd', center=1, ax=ax2, cbar_kws={'label': 'Specificity Score'})
    ax2.set_title('miRNA Specificity Scores by Body Fluid')
    ax2.set_xlabel('Body Fluid')
    ax2.set_ylabel('miRNA')
    
    # 3. Mean expression patterns
    ax3 = plt.subplot(3, 2, 3)
    
    # Select diverse miRNAs
    mean_matrix = []
    mirna_list = []
    
    for fluid in fluids:
        # Get top 3 specific for each fluid
        top = signatures.nlargest(3, f'{fluid}_specificity')['miRNA'].tolist()
        mirna_list.extend(top)
    
    mirna_list = list(set(mirna_list))[:15]
    
    for mirna in mirna_list:
        row = signatures[signatures['miRNA'] == mirna].iloc[0]
        mean_matrix.append([row[f'{fluid}_mean'] for fluid in fluids])
    
    mean_df = pd.DataFrame(mean_matrix, index=mirna_list, columns=fluids)
    
    sns.heatmap(mean_df, cmap='RdBu_r', center=0, ax=ax3, 
                annot=True, fmt='.1f', cbar_kws={'label': 'Mean Expression'})
    ax3.set_title('Mean Expression of Fluid-Specific miRNAs')
    ax3.set_xlabel('Body Fluid')
    ax3.set_ylabel('miRNA')
    
    # 4. Fluid-specific markers
    ax4 = plt.subplot(3, 2, 4)
    
    fluid_markers = find_fluid_specific_markers(signatures, threshold=1.5)
    
    # Create presence/absence matrix
    all_markers = []
    for markers in fluid_markers.values():
        all_markers.extend(markers)
    all_markers = list(set(all_markers))[:20]
    
    presence_matrix = []
    for marker in all_markers:
        row = []
        for fluid in fluids:
            if marker in fluid_markers.get(fluid, []):
                row.append(1)
            else:
                row.append(0)
        presence_matrix.append(row)
    
    presence_df = pd.DataFrame(presence_matrix, index=all_markers, columns=fluids)
    
    sns.heatmap(presence_df, cmap='Blues', ax=ax4, cbar=False, 
                annot=True, fmt='d', linewidths=0.5)
    ax4.set_title('Fluid-Specific Marker Presence')
    ax4.set_xlabel('Body Fluid')
    ax4.set_ylabel('miRNA')
    
    # 5. Expression profiles of top markers
    ax5 = plt.subplot(3, 1, 3)
    
    # Get most specific markers overall
    max_specs = []
    for idx, row in signatures.iterrows():
        max_spec = max([row[f'{fluid}_specificity'] for fluid in fluids])
        max_specs.append((row['miRNA'], max_spec))
    
    top_markers = sorted(max_specs, key=lambda x: x[1], reverse=True)[:10]
    top_mirnas = [m[0] for m in top_markers]
    
    # Plot expression profiles
    x = np.arange(len(fluids))
    width = 0.08
    
    for i, mirna in enumerate(top_mirnas):
        row = signatures[signatures['miRNA'] == mirna].iloc[0]
        means = [row[f'{fluid}_mean'] for fluid in fluids]
        ax5.bar(x + i*width, means, width, label=mirna)
    
    ax5.set_xlabel('Body Fluid')
    ax5.set_ylabel('Mean Expression')
    ax5.set_title('Expression Profiles of Top Specific miRNAs')
    ax5.set_xticks(x + width * 4.5)
    ax5.set_xticklabels(fluids, rotation=45, ha='right')
    ax5.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig('results/forensic_signatures_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return fluid_markers

def generate_forensic_report(signatures, fluid_markers):
    """Generate a proper forensic analysis report"""
    
    report = []
    report.append("# Forensic miRNA Signature Analysis\n")
    report.append("## Overview")
    report.append("Analysis of body fluid-specific miRNA expression patterns\n")
    
    fluids = ['menstrual_blood', 'peripheral_blood', 'saliva', 'semen', 'vaginal_secretion']
    
    for fluid in fluids:
        report.append(f"\n## {fluid.replace('_', ' ').title()}")
        
        # Get top specific markers
        top_specific = signatures.nlargest(5, f'{fluid}_specificity')
        
        report.append(f"\n### Top Specific miRNAs:")
        for _, row in top_specific.iterrows():
            report.append(f"- **{row['miRNA']}**: Specificity={row[f'{fluid}_specificity']:.2f}, "
                         f"Mean Expression={row[f'{fluid}_mean']:.2f}")
        
        # Expression characteristics
        high_expr = signatures[signatures[f'{fluid}_mean'] > 0].nlargest(5, f'{fluid}_mean')
        report.append(f"\n### Highest Expression:")
        for _, row in high_expr.iterrows():
            report.append(f"- {row['miRNA']}: {row[f'{fluid}_mean']:.2f}")
    
    # Save report
    with open('results/forensic_signatures_report.md', 'w') as f:
        f.write('\n'.join(report))
    
    print("Forensic signature analysis complete!")
    print(f"Visualizations saved to: results/forensic_signatures_analysis.png")
    print(f"Report saved to: results/forensic_signatures_report.md")

def main():
    """Run the forensic signature analysis"""
    # Load data
    expr, metadata = load_gpr_data()
    
    # Calculate signatures
    signatures = calculate_fluid_signatures(expr, metadata)
    
    # Create visualizations
    Path('results').mkdir(exist_ok=True)
    fluid_markers = create_forensic_visualizations(expr, metadata, signatures)
    
    # Generate report
    generate_forensic_report(signatures, fluid_markers)
    
    # Print summary
    print("\nFluid-Specific Markers (Specificity > 1.5x):")
    for fluid, markers in fluid_markers.items():
        print(f"\n{fluid}: {len(markers)} markers")
        print(f"  Top 3: {', '.join(markers[:3])}")

if __name__ == "__main__":
    main()