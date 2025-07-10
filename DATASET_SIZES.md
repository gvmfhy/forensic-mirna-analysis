# Dataset Sizes and Analysis Potential

## What We Have

### GSE153135 (GenePix GPR format) - ANALYZED
- **Saliva**: 2 samples
- **Peripheral blood**: 2 samples  
- **Semen**: 2 samples
- **Vaginal secretion**: 2 samples
- **Menstrual blood**: 2 samples
- **Total**: 10 samples (ALL available samples analyzed)

### GSE49630 (Affymetrix CEL format) - NOT YET ANALYZED
- **Blood**: 5 samples
- **Saliva**: 5 samples
- **Semen**: 5 samples  
- **Vaginal**: 5 samples
- **Menstrual blood**: 0 samples (not included in this dataset)
- **Total**: 20 samples

## Statistical Implications

### Current Analysis (n=2 per fluid)
- Too small for machine learning
- No statistical significance possible
- Results are essentially anecdotal
- Cross-validation impossible without overfitting

### Potential with CEL Data (n=5 per fluid)
- Still small but 2.5x improvement
- Could use 3-fold cross-validation
- Better chance of finding patterns
- Still below ideal forensic standards (n=20-50)

### Combined Analysis Potential
- **Blood**: 7 samples total (2 GPR + 5 CEL)
- **Saliva**: 7 samples total
- **Semen**: 7 samples total
- **Vaginal**: 7 samples total
- **Menstrual**: 2 samples (GPR only)

## Key Insight
The analysis pipeline we developed for GPR data can be applied to CEL data once properly processed through R/Bioconductor. The CEL dataset offers:
1. More samples per fluid type
2. No menstrual blood (limitation)
3. Different platform requiring cross-platform normalization

## Recommendation
Process CEL data and either:
1. Analyze CEL separately (n=5 per fluid)
2. Attempt cross-platform integration (n=7 per fluid for 4 fluids)
3. Use one dataset for training, other for validation