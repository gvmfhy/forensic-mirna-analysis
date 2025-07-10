# oligo vs affy Package Comparison for Affymetrix miRNA Arrays

## Executive Summary

**For Affymetrix miRNA 3.0/4.0 arrays, use the `oligo` package, NOT the `affy` package.**

## Detailed Comparison

### 1. **Array Support**

| Feature | affy | oligo |
|---------|------|-------|
| 3' IVT arrays (e.g., HG-U133) | ✓ | ✓ |
| Gene ST arrays | ✗ | ✓ |
| Exon ST arrays | ✗ | ✓ |
| miRNA arrays | ✗ | ✓ |
| SNP arrays | ✗ | ✓ |
| Tiling arrays | ✗ | ✓ |

### 2. **CEL File Format Support**

- **affy**: Only supports text-based CEL files (GCOS/MAS5 format)
- **oligo**: Supports both text-based and binary CEL files (Calvin/AGCC format)

The GSE49630 dataset uses **Calvin format** CEL files, which means:
```r
# This will FAIL with affy:
library(affy)
data <- ReadAffy(filenames = celFiles)  # ERROR!

# This will WORK with oligo:
library(oligo)
data <- read.celfiles(celFiles)  # Success!
```

### 3. **Platform Design (PD) Packages**

**affy** uses CDF (Chip Definition File) packages:
```r
# Example for older arrays
library(hgu133plus2cdf)
```

**oligo** uses Platform Design (PD) packages:
```r
# For miRNA 3.0 arrays
library(pd.mirna.3.0)
# For miRNA 4.0 arrays
library(pd.mirna.4.0)
```

### 4. **Memory Efficiency**

oligo is more memory-efficient for large datasets:

```r
# Memory usage comparison (approximate)
# For 20 miRNA arrays:
# affy: ~2-3 GB RAM (if it could read the files)
# oligo: ~1-2 GB RAM
```

### 5. **Normalization Methods**

Both packages support RMA, but with different implementations:

```r
# affy (won't work for miRNA arrays)
eset <- rma(affybatch)

# oligo (recommended)
eset <- rma(oligobatch)
# More options available:
eset <- rma(oligobatch, 
            background = TRUE,
            normalize = TRUE, 
            target = "probeset")  # or "core", "extended", "full"
```

### 6. **Probe Set Definitions**

For miRNA arrays, oligo handles multiple probe set types better:

```r
# oligo can handle different annotation levels
getNetAffx(eset, "probeset")  # Default
getNetAffx(eset, "transcript")  # Transcript level
```

### 7. **Code Example: Why affy Fails**

```r
# Attempting to use affy with miRNA 3.0 arrays
library(affy)

# This will fail at multiple points:
tryCatch({
  # 1. CEL file reading fails
  celFiles <- list.files("data/raw/GSE49630_CEL", 
                        pattern = "*.CEL.gz$", 
                        full.names = TRUE)
  data <- ReadAffy(filenames = celFiles)
}, error = function(e) {
  print("Error: Cannot read Calvin format CEL files")
})

# Even if it could read the files:
# 2. No CDF package exists for miRNA arrays
# 3. Probe set definitions incompatible
# 4. RMA implementation doesn't support ST array structure
```

### 8. **Platform-Specific Features**

oligo provides features specifically for ST arrays:

```r
# Background correction options specific to ST arrays
eset <- rma(data, target = "core")  # Core probe sets only

# Access to different probe types
pmSeq <- pmSequence(data)  # Get PM probe sequences

# Quality assessment specific to ST arrays
bbox <- boxplot(data, target = "core")
```

### 9. **Integration with Bioconductor**

Both integrate well, but oligo has better support for modern workflows:

```r
# oligo works seamlessly with newer Bioconductor packages
library(oligo)
library(arrayQualityMetrics)
library(limma)

# Standard Bioconductor workflow
eset <- rma(data)
arrayQualityMetrics(eset)  # Works perfectly
```

### 10. **When to Use Each Package**

**Use affy when:**
- Working with older 3' IVT arrays (HG-U133, HG-U95, MG-U74, etc.)
- CEL files are in text format (GCOS)
- You need MAS5 preprocessing
- Working with legacy code/pipelines

**Use oligo when:**
- Working with any ST array (Gene ST, Exon ST, miRNA ST)
- Working with SNP arrays
- CEL files are in binary format (Calvin/AGCC)
- Need memory-efficient processing
- Working with modern Affymetrix arrays (post-2007)

## Conclusion

For the GSE49630 dataset (Affymetrix miRNA 3.0 ST arrays), **oligo is the only viable option**. The affy package simply cannot process these arrays due to:

1. Incompatible CEL file format (Calvin vs GCOS)
2. Lack of CDF support for miRNA arrays
3. Different probe set structure
4. Missing ST array-specific features

Always use oligo for modern Affymetrix arrays, including all miRNA arrays.