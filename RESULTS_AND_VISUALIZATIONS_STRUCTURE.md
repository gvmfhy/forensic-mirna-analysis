# Results and Visualizations File Structure

## Overview
All analysis results are organized in the `results/` directory with clear subdirectories for each analysis type.

## Directory Structure

```
results/
├── continuous_expression_analysis/     # GPR continuous expression analysis
│   ├── expression_tiers.csv           # Multi-tier expression levels for each miRNA
│   ├── wilcoxon_results.csv          # Statistical test results with p-values and fold changes
│   ├── specificity_scores.csv        # Fluid-specific expression scores
│   ├── forensic_markers.csv          # Identified forensic marker candidates
│   ├── continuous_expression_report.md # Detailed analysis report
│   └── visualizations/
│       ├── expression_tiers_heatmap.png    # Heatmap of tiered expression levels
│       ├── expression_profiles_boxplot.png  # Box plots by body fluid
│       ├── specificity_volcano.png         # Volcano plots of fold change vs significance
│       └── forensic_markers_barplot.png    # Top markers visualization
│
├── cel_continuous_expression/         # CEL continuous expression analysis
│   ├── expression_matrix.csv         # Normalized expression values
│   ├── wilcoxon_results.csv         # Statistical test results
│   ├── specificity_scores.csv       # Fluid-specific scores
│   ├── forensic_candidates.csv      # Marker candidates meeting criteria
│   └── cel_continuous_report.md     # Analysis report
│
├── cel_forensic_analysis/            # CEL forensic-focused analysis (FDR)
│   ├── forensic_markers_fdr.csv     # Markers passing FDR correction
│   ├── forensic_analysis_report.md  # Report (0 markers passed strict FDR)
│   ├── forensic_markers_heatmap.png # Visualization (if markers found)
│   └── forensic_markers_by_fluid.png # Per-fluid visualization
│
└── cel_practical_forensic/           # CEL practical forensic analysis
    ├── forensic_candidates.csv       # 393 marker candidates
    ├── practical_forensic_report.md  # Detailed marker descriptions
    ├── top_markers_heatmap.png       # Top markers heatmap
    ├── effect_sizes.png              # Cohen's d vs fold change scatter
    └── expression_comparison.png     # Expression level comparisons

```

## Key Results Files

### 1. GPR Analysis (Pooled Samples)
- **Main Results**: `results/continuous_expression_analysis/forensic_markers.csv`
  - 36 total markers identified
  - Blood: 5 markers, Semen: 8 markers, Saliva: 11 markers, Vaginal: 12 markers
- **Key Visualization**: `expression_tiers_heatmap.png` showing multi-tier detection framework

### 2. CEL Analysis (Individual Samples)
- **Main Results**: `results/cel_practical_forensic/forensic_candidates.csv`
  - 393 total marker candidates
  - Blood: 155 markers, Semen: 167 markers, Saliva: 49 markers, Vaginal: 22 markers
- **Key Visualization**: `top_markers_heatmap.png` showing fluid-specific expression patterns

### 3. Statistical Reports
- **GPR Report**: `results/continuous_expression_analysis/continuous_expression_report.md`
  - Explains pooled sample limitations
  - Details multi-tier framework rationale
- **CEL Report**: `results/cel_practical_forensic/practical_forensic_report.md`
  - Lists top 5 markers per fluid with detailed statistics
  - Explains why strict FDR < 0.05 not achieved with n=5

## Visualization Highlights

### Best Visualizations for Presentations:
1. **`cel_practical_forensic/top_markers_heatmap.png`**
   - Shows clear fluid-specific patterns
   - Includes fold changes and significance levels
   - Professional appearance for forensic context

2. **`cel_practical_forensic/effect_sizes.png`**
   - Demonstrates large effect sizes (Cohen's d > 1.5)
   - Shows practical significance beyond p-values
   - Color-coded by body fluid

3. **`continuous_expression_analysis/expression_tiers_heatmap.png`**
   - Illustrates multi-tier detection concept
   - Shows continuous expression nature
   - Good for explaining methodology

## Data Formats

### CSV Files
All CSV files include headers with clear column names:
- `mirna`: miRNA identifier
- `fluid`: Target body fluid
- `pvalue`: Raw p-value
- `padj` or `fdr`: Adjusted p-value
- `log2fc`: Log2 fold change
- `cohens_d`: Effect size
- `fluid_mean`: Mean expression in target fluid
- `other_mean`: Mean expression in other fluids
- `detection_rate`: Fraction of samples with expression

### Markdown Reports
Reports follow consistent structure:
1. Executive Summary
2. Results by Body Fluid
3. Statistical Considerations
4. Forensic Implementation Recommendations
5. Cross-Platform Validation (where applicable)

## Usage Notes

1. **For Publication**: Use results from `cel_practical_forensic/` as they represent individual samples with proper statistics
2. **For Method Development**: Reference `continuous_expression_analysis/` for the multi-tier framework concept
3. **For Forensic Implementation**: Focus on markers with >100-fold enrichment and Cohen's d > 3.0
4. **For Validation**: Cross-reference markers appearing in both GPR and CEL analyses

## File Size Summary
- Total results size: ~15 MB
- Largest files: Expression matrices (~2-3 MB each)
- Visualizations: ~200-500 KB per PNG file
- Reports: ~10-50 KB per markdown file

## Related Documentation
- **Detailed Analysis Walkthrough**: See `DETAILED_ANALYSIS_WALKTHROUGH.md` for step-by-step explanation of the entire analysis process, including errors and corrections