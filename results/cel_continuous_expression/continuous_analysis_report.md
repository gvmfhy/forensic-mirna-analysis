# Continuous Expression Analysis Report (GPR Data)

## Dataset Overview
- Total samples: 20
- Total miRNAs: 22676
- Expression range: -0.15 to 15.02 (log2)

## Sample Distribution
- blood: 5 samples (pooled from 12 donors)
- semen: 5 samples (pooled from 12 donors)
- vaginal: 5 samples (pooled from 12 donors)
- saliva: 5 samples (pooled from 12 donors)

## Expression Tier Distribution
- High confidence (>-2.0): 453520 (100.0%)
- Moderate confidence: 0 (0.0%)
- Low confidence: 0 (0.0%)
- Below detection: 0 (0.0%)

## Statistical Analysis Results
- Total comparisons: 90704
- Significant (p<0.05 before correction): 18999
- Significant (Bonferroni corrected): 0

## Forensic Marker Candidates

**WARNING**: No markers met all forensic criteria.
This is expected with n=2 pooled samples per fluid.

Relaxed criteria results (p<0.05, |log2FC|>1):

### Blood
- hsa-miR-500a-star_st: log2FC=5.86, p=0.001
- mmu-miR-151-3p_st: log2FC=2.45, p=0.001
- mmu-miR-134_st: log2FC=2.85, p=0.001

### Saliva
- hvt-miR-H18-5p_st: log2FC=-4.74, p=0.001
- lca-miR-23a_st: log2FC=5.59, p=0.001
- lgi-miR-375_st: log2FC=3.70, p=0.001

### Semen
- dme-miR-974-5p_st: log2FC=3.38, p=0.001
- dme-miR-986-3p_st: log2FC=2.78, p=0.001
- dmo-miR-210_st: log2FC=-6.91, p=0.001

## Limitations
- Small sample size (n=2 arrays per fluid) limits statistical power
- Arrays represent pooled samples, reducing biological variability
- Results require validation in larger, individual sample cohorts
- CEL data analysis recommended for more robust findings
## CEL-Specific Notes
- Platform: Affymetrix miRNA 3.0 ST arrays
- Processing: oligo::rma() normalization
- Sample type: Individual donors (not pooled)
- Statistical power: Improved with n=5 per fluid
