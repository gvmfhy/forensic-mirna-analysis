# Forensic miRNA Analysis Plan: Cross-Platform Integration

## Essential Elements of Information (EEIs)

1. **Minimum Viable Panel**: Identify smallest miRNA set achieving >99% specificity
2. **Mixture Resolution**: Detect minor contributors at 1:10 to 1:100 ratios
3. **Platform Independence**: Find markers stable across CEL/GPR platforms
4. **Degradation Robustness**: Select miRNAs surviving real-world conditions
5. **Implementation Feasibility**: Ensure deployability in forensic labs

## Data Overview

### GSE49630 (Affymetrix CEL)
- **Platform**: miRNA-3_0 array
- **Samples**: 20 (5 each: blood, saliva, semen, vaginal)
- **Format**: Binary CEL files with annotation CSV
- **Missing**: Menstrual blood

### GSE153135 (Two-color GPR)
- **Platform**: miRCURY LNA array (Exiqon)
- **Samples**: 10 (pools from ~13 donors)
- **Format**: Text-based GPR with dual-channel intensities
- **Includes**: Menstrual blood (critical addition)

## Selected Approach: Integrated Multi-Platform Ensemble

### Phase 1: Data Preprocessing
1. **CEL Processing**
   - Use `affy` package for reading
   - RMA normalization
   - Log2 transformation

2. **GPR Processing**
   - Extract Cy5/Cy3 intensities
   - Loess normalization within arrays
   - Quantile normalization between arrays

3. **Annotation Mapping**
   - Map probes to miRBase v21
   - Retain only shared miRNAs
   - Handle multiple probes per miRNA

### Phase 2: Cross-Platform Integration
1. **Initial QC**
   - PCA by platform and fluid type
   - Identify platform-specific batch effects

2. **Batch Correction Strategy**
   - Combat with platform as batch
   - Preserve biological variation (fluid types)
   - Validate with held-out samples

3. **Feature Selection**
   - limma differential expression
   - Information gain ranking
   - Boruta wrapper selection
   - Consensus: miRNAs in ≥2/3 methods

### Phase 3: Ensemble Classifier
1. **Base Models**
   - Random Forest (captures interactions)
   - XGBoost (handles imbalanced classes)
   - Elastic Net (interpretable coefficients)

2. **Ensemble Strategy**
   - Soft voting with optimized weights
   - Cross-validation within platforms
   - Leave-one-platform-out validation

3. **Mixture Handling**
   - Generate synthetic mixtures (1:1 to 1:100)
   - Train mixture deconvolution model
   - Validate detection limits

### Phase 4: Forensic Validation
1. **Marker Panel**
   - Select top 15 consensus miRNAs
   - Include miR-191-5p as housekeeping
   - Document platform-specific biases

2. **Performance Metrics**
   - Specificity >99% (court requirement)
   - Sensitivity by fluid type
   - Mixture detection thresholds
   - ROC curves with confidence intervals

3. **Implementation Guide**
   - qPCR primer design priorities
   - Sample QC thresholds
   - Decision tree for interpretation

## Technical Implementation

### Environment Setup
```bash
# R packages
install.packages(c("tidyverse", "caret", "xgboost", "glmnet"))
BiocManager::install(c("affy", "limma", "sva", "GEOquery"))

# Python packages
pip install lightgbm scikit-learn pandas numpy matplotlib
```

### File Organization
```
analysis/
├── 01_preprocessing.R      # Read CEL/GPR, normalize
├── 02_integration.R        # Batch correction, QC
├── 03_feature_selection.R  # Multi-method selection
├── 04_ensemble_model.R     # Train classifiers
├── 05_mixture_analysis.R   # Synthetic mixtures
└── 06_report_generation.R  # Final results

scripts/
├── utils/
│   ├── cel_reader.R
│   ├── gpr_parser.R
│   └── plotting_functions.R
└── models/
    ├── ensemble_trainer.py
    └── mixture_deconvolution.py
```

## Risk Mitigation

1. **Platform Bias**: Use leave-platform-out validation
2. **Overfitting**: Nested cross-validation, simple models
3. **Class Imbalance**: SMOTE for minority classes
4. **Reproducibility**: Set seeds, document versions

## Success Criteria

- [ ] Cross-platform correlation >0.8 for top markers
- [ ] Individual fluid specificity >99%
- [ ] Mixture detection at 1:10 ratio >90% accuracy
- [ ] Complete analysis in 3 days
- [ ] Generate publication-ready figures

## Next Steps

1. Parse both data formats into common structure
2. Perform initial QC and visualization
3. Implement preprocessing pipeline
4. Begin feature selection process