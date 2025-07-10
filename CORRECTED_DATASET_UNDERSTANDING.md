# Corrected Dataset Understanding

## Critical Update: Pooled vs Individual Samples

### GSE153135 (GPR) - POOLED SAMPLES
- **10 arrays** = what we download and analyze
- **26 individual donors** represented in those pools
- Each array is a pool of 2-3 donors:
  - Array 1 per fluid: 3 donors pooled
  - Array 2 per fluid: 2-3 donors pooled
- Breakdown by fluid:
  - Peripheral blood: 5 donors → 2 arrays
  - Saliva: 5 donors → 2 arrays
  - Semen: 5 donors → 2 arrays
  - Vaginal: 5 donors → 2 arrays
  - Menstrual blood: 6 donors → 2 arrays

### GSE49630 (CEL) - INDIVIDUAL SAMPLES
- **20 arrays** = 20 individual donors (1:1 mapping)
- No pooling
- Breakdown by fluid:
  - Blood: 5 individual arrays
  - Saliva: 5 individual arrays
  - Semen: 5 individual arrays
  - Vaginal: 5 individual arrays

## Total: 46 Individual Donors Represented

## Critical Implications

### For GPR Analysis:
1. **Reduced biological variability** - pooling averages out individual differences
2. **Cannot do donor-level statistics** - we have array-level data only
3. **Cross-validation must be by array** (n=2 per fluid), not by donor
4. **May see stronger signals** - pooling can reduce noise
5. **May miss rare variants** - pooling dilutes individual-specific markers

### For Statistical Analysis:
- GPR: n=2 observations per fluid (arrays), representing 5-6 donors each
- CEL: n=5 observations per fluid (individual donors)
- **Cannot treat GPR arrays as if they were individual samples**

### For Method Reporting:
- Must state: "26 individual samples profiled on 10 pooled microarrays"
- Not: "10 samples analyzed"

## Analysis Implications

### Why Binary Analysis May Be Limiting:
1. **Biology is continuous** - miRNAs act as dimmer switches, not on/off
2. **Forensic degradation** - samples may have reduced but detectable levels
3. **Pooling effects** - GPR data already averaged, binary cutoff loses information
4. **Platform differences** - different dynamic ranges need different thresholds

### Better Approaches:
1. **Continuous expression analysis** with robust statistics
2. **Rank-based methods** for cross-platform comparison
3. **Effect size focus** rather than just p-values
4. **Detection confidence scores** rather than binary calls