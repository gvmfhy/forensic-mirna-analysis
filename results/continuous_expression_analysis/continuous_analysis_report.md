# Continuous Expression Analysis Report (GPR Data)

## Dataset Overview
- Total samples: 10
- Total miRNAs: 2057
- Expression range: -14.80 to 0.50 (log2)

## Sample Distribution
- saliva: 2 samples (pooled from 5 donors)
- blood: 2 samples (pooled from 5 donors)
- semen: 2 samples (pooled from 5 donors)
- vaginal: 2 samples (pooled from 5 donors)
- menstrual: 2 samples (pooled from 5 donors)

## Expression Tier Distribution
- High confidence (>-2.0): 2960 (14.4%)
- Moderate confidence: 8281 (40.3%)
- Low confidence: 7550 (36.7%)
- Below detection: 1779 (8.6%)

## Statistical Analysis Results
- Total comparisons: 10285
- Significant (p<0.05 before correction): 1122
- Significant (Bonferroni corrected): 0

## Forensic Marker Candidates

**WARNING**: No markers met all forensic criteria.
This is expected with n=2 pooled samples per fluid.

Relaxed criteria results (p<0.05, |log2FC|>1):

### Saliva
- ebv-miR-BART10-3p: log2FC=1.46, p=0.037
- hsa-miR-3115: log2FC=-1.19, p=0.037
- hsa-miR-3130-3p: log2FC=-1.36, p=0.037

### Vaginal
- hsa-miR-4764-3p: log2FC=1.04, p=0.037
- hsa-miR-4767: log2FC=-1.16, p=0.037
- hsa-miR-4779: log2FC=-2.38, p=0.037

### Menstrual
- hsa-miR-4662a-5p: log2FC=1.17, p=0.037
- hsa-miR-4678: log2FC=2.70, p=0.037
- hsa-miR-4684-5p: log2FC=2.44, p=0.037

## Limitations
- Small sample size (n=2 arrays per fluid) limits statistical power
- Arrays represent pooled samples, reducing biological variability
- Results require validation in larger, individual sample cohorts
- CEL data analysis recommended for more robust findings