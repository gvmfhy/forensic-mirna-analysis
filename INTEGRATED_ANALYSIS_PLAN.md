# Integrated Analysis Plan for Forensic miRNA

## Key Insights from Multi-Agent Analysis

### 1. Data Structure Corrections
- **GPR**: 10 arrays representing 26 pooled donors (not suitable for individual-level statistics)
- **CEL**: 20 arrays representing 20 individual donors (suitable for proper statistics)
- Focus on CEL for discovery, use GPR for validation only

### 2. Bioconductor Processing (CEL)
- Use `oligo::rma()` NOT `affy::rma()`
- Affymetrix miRNA ST arrays require oligo package
- Include QC plots: PCA, boxplots, correlation heatmaps
- Filter for human miRNAs (hsa-miR- prefix)

### 3. Statistical Approach
- **Primary**: Wilcoxon rank-sum tests (non-parametric, n=5 suitable)
- **Correction**: Bonferroni for forensic stringency (not FDR)
- **Effect size**: Require |log2FC| > 2 AND Cohen's d > 0.8
- **No machine learning** with these sample sizes

### 4. Continuous Expression Framework
Instead of binary presence/absence:
- **High confidence**: > -2.0 log2
- **Moderate confidence**: -6.0 to -2.0 log2
- **Low confidence**: -10.0 to -6.0 log2
- **Below detection**: < -10.0 log2

### 5. Forensic Marker Criteria
A marker must show:
- >5-fold higher in target fluid (log2FC > 2.32)
- Detected in >80% of target fluid samples
- AUC > 0.99 for specificity
- Biological plausibility

## Immediate Action Plan

### Phase 1: CEL Processing (Priority)
```R
# Use the oligo-based script from agent
library(oligo)
celFiles <- list.celfiles("data/raw/GSE49630_CEL", full.names=TRUE)
rawData <- read.celfiles(celFiles)
eset <- rma(rawData)
```

### Phase 2: Quality Assessment
1. Run PCA to check sample clustering by fluid
2. Calculate detection rates per miRNA
3. Assess expression distributions
4. Identify and remove outliers if necessary

### Phase 3: Statistical Analysis
```
For each fluid vs all others:
1. Wilcoxon test with Bonferroni correction
2. Calculate fold changes and Cohen's d
3. Apply forensic filters (FC, effect size, detection rate)
4. Calculate specificity scores and AUC
```

### Phase 4: Continuous Expression Analysis
1. Implement multi-tier confidence system
2. Calculate likelihood ratios (not binary calls)
3. Create gradient visualizations
4. Account for sample quality/degradation

### Phase 5: Validation
1. Check if CEL markers show same direction in GPR
2. Compare effect sizes between platforms
3. Literature validation of top candidates

## Why This Approach?

1. **Scientifically Sound**: Respects continuous biology, proper statistics for small n
2. **Forensically Appropriate**: High specificity, conservative corrections, court-defensible
3. **Practical**: Works with actual sample sizes, accounts for pooling
4. **Transparent**: Clear about limitations, probabilistic reporting

## Expected Outcomes

With n=5 (CEL), we can likely identify:
- 2-5 high-confidence markers per fluid
- Effect sizes large enough for forensic use
- Probabilistic models for mixed samples
- Clear quality thresholds for degraded samples

## Key Differences from Previous Approach

1. No premature machine learning
2. Treat expression as continuous, not binary
3. CEL for discovery, GPR for validation only
4. Conservative statistics appropriate for court
5. Account for pooled nature of GPR data