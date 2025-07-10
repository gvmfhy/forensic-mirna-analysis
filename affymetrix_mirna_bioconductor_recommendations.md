# Bioconductor Analysis Recommendations for Affymetrix miRNA Arrays (GSE49630)

## Dataset Overview
- **Platform**: Affymetrix miRNA-3.0 ST Array (GPL16384)
- **Array Type**: miRNA-3_1 (based on annotation file header)
- **Samples**: 20 CEL files (5 replicates each of Blood, Semen, Vaginal, Saliva)
- **miRBase Version**: v17 (as per annotation file)

## Key Recommendations

### 1. **Recommended Bioconductor Package: oligo**

The **oligo** package is the preferred choice for Affymetrix miRNA arrays for several reasons:

```r
# Installation
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("oligo")
BiocManager::install("pd.mirna.3.0")  # Platform design package
```

**Why oligo over affy?**
- The `affy` package is designed for older 3' IVT arrays (e.g., HG-U133)
- `oligo` supports newer Affymetrix arrays including ST arrays and miRNA arrays
- `oligo` handles the Calvin CEL file format (AGCC) used by these arrays
- Better memory efficiency for large datasets

### 2. **Processing Workflow**

```r
library(oligo)
library(pd.mirna.3.0)

# Read CEL files
celFiles <- list.files("data/raw/GSE49630_CEL", pattern = "*.CEL.gz$", full.names = TRUE)
rawData <- read.celfiles(celFiles)

# Background correction, normalization, and summarization
eset <- rma(rawData)

# Alternative: if you need more control
eset <- rma(rawData, 
            background = TRUE,
            normalize = TRUE,
            target = "probeset")
```

### 3. **Quality Control Steps**

Essential QC steps for miRNA arrays:

```r
# 1. Array quality metrics
library(arrayQualityMetrics)
arrayQualityMetrics(rawData, outdir = "QC_report")

# 2. MA plots
MAplot(rawData, pairs = TRUE)

# 3. Density plots
hist(rawData, main = "Density plots of raw intensities")

# 4. Box plots
boxplot(rawData, main = "Boxplot of raw intensities")

# 5. RNA degradation (less relevant for miRNA)
# miRNAs are short and stable, so traditional RNA degradation plots are not applicable

# 6. Principal Component Analysis
plotPCA(eset, groups = pData(eset)$tissue_type)

# 7. Hierarchical clustering
library(gplots)
heatmap.2(cor(exprs(eset)), trace = "none")
```

### 4. **Probe-to-miRNA Mapping**

The annotation file contains multiple probe types:

```r
# Load annotation
library(data.table)
annot <- fread("GPL16384_miRNA-3_1-st-v1.annotations.20140513.csv.gz", skip = 4)

# Filter for different RNA types
mirna_probes <- annot[`Sequence Type` == "miRNA"]
snoRNA_probes <- annot[`Sequence Type` %in% c("CDBox", "HAcaBox", "scaRna")]
pre_mirna_probes <- annot[`Sequence Type` == "stem-loop"]

# Map probes to miRNAs
featureData(eset) <- new("AnnotatedDataFrame", 
                         data = annot[match(featureNames(eset), annot$`Probe Set ID`),])

# Filter for human miRNAs only
human_mirnas <- grep("^hsa-", featureNames(eset), value = TRUE)
eset_human <- eset[human_mirnas,]
```

### 5. **Special Considerations for miRNA Arrays**

#### a. **Multiple probe types**
- The array contains probes for mature miRNAs, pre-miRNAs (stem-loops), and other small RNAs
- Decide whether to analyze all or filter for specific types

#### b. **Cross-hybridization**
- miRNAs from the same family can cross-hybridize due to sequence similarity
- Consider using probe sets with "_st" suffix which are typically more specific

#### c. **Low expression filtering**
```r
# Remove lowly expressed miRNAs
library(genefilter)
eset_filtered <- nsFilter(eset_human, 
                         require.entrez = FALSE,
                         remove.dupEntrez = FALSE,
                         var.filter = FALSE,
                         feature.exclude = "^AFFX")$eset
```

#### d. **Normalization considerations**
- RMA is generally appropriate for miRNA arrays
- Alternative: quantile normalization with background correction
- Some studies use VSN (variance stabilizing normalization) for miRNA data

```r
# Alternative normalization
library(vsn)
eset_vsn <- vsnrma(rawData)
```

#### e. **Batch effects**
- Check for batch effects between tissue types
- Use ComBat or limma's removeBatchEffect if needed

```r
library(sva)
batch <- pData(eset)$scan_date  # if available
modcombat <- model.matrix(~tissue_type, data = pData(eset))
combat_eset <- ComBat(dat = exprs(eset), batch = batch, mod = modcombat)
```

### 6. **Differential Expression Analysis**

```r
library(limma)

# Create design matrix
tissue_type <- factor(pData(eset)$tissue_type)
design <- model.matrix(~0 + tissue_type)
colnames(design) <- levels(tissue_type)

# Fit linear model
fit <- lmFit(eset, design)

# Define contrasts
contrast.matrix <- makeContrasts(
  Blood_vs_Saliva = Blood - Saliva,
  Blood_vs_Semen = Blood - Semen,
  Blood_vs_Vaginal = Blood - Vaginal,
  levels = design
)

# Apply contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get results
results <- topTable(fit2, coef = "Blood_vs_Saliva", adjust = "BH", number = Inf)
```

### 7. **Data Export and Visualization**

```r
# Export normalized expression matrix
write.csv(exprs(eset), "normalized_expression.csv")

# Volcano plot
volcanoplot(fit2, coef = "Blood_vs_Saliva", highlight = 10)

# Heatmap of top differentially expressed miRNAs
top_mirnas <- rownames(results)[1:50]
heatmap.2(exprs(eset)[top_mirnas,], 
          scale = "row",
          trace = "none",
          col = bluered(100))
```

## Summary of Key Differences from mRNA Arrays

1. **Probe design**: miRNA probes are shorter and may have different hybridization kinetics
2. **Dynamic range**: miRNA expression typically has a wider dynamic range
3. **Annotation complexity**: Multiple probe types (mature, precursor, star sequences)
4. **Background correction**: May need adjustment due to shorter probe sequences
5. **Quality metrics**: Traditional RNA degradation plots not applicable
6. **Normalization**: Consider miRNA-specific normalization methods
7. **Filtering**: More stringent filtering may be needed due to cross-hybridization

## Recommended Complete Workflow

```r
# 1. Load libraries
library(oligo)
library(pd.mirna.3.0)
library(limma)
library(data.table)

# 2. Read data
celFiles <- list.files("data/raw/GSE49630_CEL", pattern = "*.CEL.gz$", full.names = TRUE)
rawData <- read.celfiles(celFiles)

# 3. Quality control
hist(rawData)
boxplot(rawData)

# 4. Normalization
eset <- rma(rawData)

# 5. Annotation
annot <- fread("GPL16384_miRNA-3_1-st-v1.annotations.20140513.csv.gz", skip = 4)
featureData(eset) <- new("AnnotatedDataFrame", 
                         data = annot[match(featureNames(eset), annot$`Probe Set ID`),])

# 6. Filter for human miRNAs
human_mirnas <- grep("^hsa-mir|^hsa-miR", featureNames(eset), value = TRUE)
eset_human <- eset[human_mirnas,]

# 7. Differential expression analysis
# ... (as shown above)
```

This approach ensures robust and reproducible analysis of Affymetrix miRNA array data.