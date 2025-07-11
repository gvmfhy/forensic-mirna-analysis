# Practical Forensic miRNA Analysis - CEL Data

## Executive Summary
With n=5 individual samples per body fluid, we can identify
scientifically meaningful markers for forensic identification.

## Blood
Found 155 strong candidates (p<0.01, |FC|>4, |d|>1.5)

### Top 5 Markers:

1. **hsa-miR-126**
   - Fold change: 476.5x (higher in blood)
   - Effect size (Cohen's d): 15.83
   - p-value: 0.0011
   - FDR q-value: 0.071
   - Expression: 10.2 vs 1.3

2. **hsa-miR-3200-3p**
   - Fold change: 57.4x (higher in blood)
   - Effect size (Cohen's d): 23.33
   - p-value: 0.0011
   - FDR q-value: 0.071
   - Expression: 6.6 vs 0.7

3. **hsa-miR-363**
   - Fold change: 196.7x (higher in blood)
   - Effect size (Cohen's d): 13.05
   - p-value: 0.0011
   - FDR q-value: 0.071
   - Expression: 8.7 vs 1.1

4. **hsa-miR-486-5p**
   - Fold change: 3243.2x (higher in blood)
   - Effect size (Cohen's d): 8.36
   - p-value: 0.0011
   - FDR q-value: 0.071
   - Expression: 13.8 vs 2.1

5. **hsa-miR-200c**
   - Fold change: 0.0x (lower in blood)
   - Effect size (Cohen's d): -13.64
   - p-value: 0.0011
   - FDR q-value: 0.071
   - Expression: 6.0 vs 12.4

## Saliva
Found 49 strong candidates (p<0.01, |FC|>4, |d|>1.5)

### Top 5 Markers:

1. **hsa-miR-223**
   - Fold change: 119.4x (higher in saliva)
   - Effect size (Cohen's d): 4.81
   - p-value: 0.0011
   - FDR q-value: 0.071
   - Expression: 9.8 vs 2.9

2. **hsa-miR-3651**
   - Fold change: 57.9x (higher in saliva)
   - Effect size (Cohen's d): 3.79
   - p-value: 0.0011
   - FDR q-value: 0.071
   - Expression: 8.8 vs 2.9

3. **hsa-miR-145**
   - Fold change: 104.6x (higher in saliva)
   - Effect size (Cohen's d): 3.49
   - p-value: 0.0019
   - FDR q-value: 0.076
   - Expression: 9.6 vs 2.9

4. **hsa-miR-939**
   - Fold change: 0.1x (lower in saliva)
   - Effect size (Cohen's d): -3.81
   - p-value: 0.0011
   - FDR q-value: 0.071
   - Expression: 4.9 vs 8.5

5. **hsa-miR-125a-5p**
   - Fold change: 20.3x (higher in saliva)
   - Effect size (Cohen's d): 3.58
   - p-value: 0.0026
   - FDR q-value: 0.084
   - Expression: 6.2 vs 1.8

## Semen
Found 167 strong candidates (p<0.01, |FC|>4, |d|>1.5)

### Top 5 Markers:

1. **hsa-miR-222**
   - Fold change: 0.0x (lower in semen)
   - Effect size (Cohen's d): -13.68
   - p-value: 0.0011
   - FDR q-value: 0.071
   - Expression: 0.9 vs 9.9

2. **hsa-miR-221**
   - Fold change: 0.0x (lower in semen)
   - Effect size (Cohen's d): -7.61
   - p-value: 0.0011
   - FDR q-value: 0.071
   - Expression: 1.1 vs 9.9

3. **hsa-miR-16**
   - Fold change: 0.0x (lower in semen)
   - Effect size (Cohen's d): -6.14
   - p-value: 0.0011
   - FDR q-value: 0.071
   - Expression: 2.6 vs 11.5

4. **hsa-miR-2392**
   - Fold change: 79.1x (higher in semen)
   - Effect size (Cohen's d): 6.45
   - p-value: 0.0011
   - FDR q-value: 0.071
   - Expression: 10.5 vs 4.1

5. **hsa-miR-483-5p**
   - Fold change: 137.9x (higher in semen)
   - Effect size (Cohen's d): 5.19
   - p-value: 0.0011
   - FDR q-value: 0.071
   - Expression: 11.0 vs 3.8

## Vaginal
Found 22 strong candidates (p<0.01, |FC|>4, |d|>1.5)

### Top 5 Markers:

1. **hsa-miR-20b**
   - Fold change: 0.0x (lower in vaginal)
   - Effect size (Cohen's d): -2.41
   - p-value: 0.0014
   - FDR q-value: 0.072
   - Expression: 1.3 vs 7.1

2. **hsa-miR-151-5p**
   - Fold change: 0.0x (lower in vaginal)
   - Effect size (Cohen's d): -2.52
   - p-value: 0.0026
   - FDR q-value: 0.084
   - Expression: 1.5 vs 7.5

3. **hsa-miR-138-1-star**
   - Fold change: 7.0x (higher in vaginal)
   - Effect size (Cohen's d): 3.79
   - p-value: 0.0011
   - FDR q-value: 0.071
   - Expression: 3.9 vs 1.1

4. **hsa-miR-1260b**
   - Fold change: 16.7x (higher in vaginal)
   - Effect size (Cohen's d): 2.42
   - p-value: 0.0011
   - FDR q-value: 0.071
   - Expression: 6.9 vs 2.9

5. **hsa-miR-574-3p**
   - Fold change: 10.5x (higher in vaginal)
   - Effect size (Cohen's d): 2.82
   - p-value: 0.0014
   - FDR q-value: 0.072
   - Expression: 10.6 vs 7.2

## Cross-Platform Validation

## Statistical Considerations

### Multiple Testing
- Total human miRNA tests: 13,564
- Minimum p-value: 0.001063
- Tests with p < 0.01: 906
- Minimum FDR q-value: 0.071

### Why Strict FDR < 0.05 Not Achieved
- With rank-based tests and n=5, minimum possible p-value is ~0.001
- Testing 13,564 miRNAs requires p < 0.0000037 for FDR < 0.05
- This is mathematically impossible with n=5

### Practical Approach
- Focus on effect size (fold change) and consistency
- Validate findings across platforms
- Use targeted assays for final forensic application

## Forensic Implementation

### Recommended Marker Panel

**Blood:**
- hsa-miR-486-5p: 3243-fold enriched
- hsa-miR-126: 476-fold enriched

**Saliva:**
- hsa-miR-205: 275-fold enriched
- hsa-miR-203: 209-fold enriched

**Semen:**
- hsa-miR-891a: 175-fold enriched
- hsa-miR-483-5p: 138-fold enriched

**Vaginal:**
- hp_hsa-mir-1299: 22-fold enriched
- hsa-miR-1915-star: 17-fold enriched

### Validation Strategy
1. Confirm expression patterns with qRT-PCR
2. Test on independent sample cohort
3. Assess stability in degraded samples
4. Develop multiplex assay for simultaneous detection
