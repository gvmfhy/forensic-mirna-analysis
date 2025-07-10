#!/usr/bin/env Rscript

# Affymetrix miRNA Array Processing Script for GSE49630
# Platform: Affymetrix miRNA-3.0 ST Array (GPL16384)

# Load required libraries
suppressPackageStartupMessages({
  library(oligo)
  library(pd.mirna.3.0)
  library(limma)
  library(data.table)
  library(ggplot2)
  library(gplots)
})

# Set working directory
setwd("/Users/austinmorrissey/geo_downloads")

# 1. Read CEL files
cat("Reading CEL files...\n")
celFiles <- list.files("data/raw/GSE49630_CEL", 
                      pattern = "*.CEL.gz$", 
                      full.names = TRUE)

# Extract sample information from filenames
sample_info <- data.frame(
  FileName = basename(celFiles),
  SampleID = gsub("_.*", "", basename(celFiles)),
  Tissue = gsub(".*_(.*)_[0-9]+\\.CEL\\.gz", "\\1", basename(celFiles)),
  Replicate = gsub(".*_([0-9]+)\\.CEL\\.gz", "\\1", basename(celFiles))
)

# Read CEL files
rawData <- read.celfiles(celFiles)

# Add phenotype data
pData(rawData) <- sample_info

# 2. Quality Control - Raw Data
cat("Generating QC plots for raw data...\n")

# Create QC directory
dir.create("QC_plots", showWarnings = FALSE)

# Density plots
pdf("QC_plots/density_plots_raw.pdf", width = 10, height = 8)
hist(rawData, main = "Density plots of raw intensities")
dev.off()

# Box plots
pdf("QC_plots/boxplot_raw.pdf", width = 10, height = 8)
boxplot(rawData, main = "Boxplot of raw intensities", las = 2)
dev.off()

# MA plots (first 4 samples)
pdf("QC_plots/MA_plots.pdf", width = 10, height = 10)
MAplot(rawData[, 1:4], pairs = TRUE)
dev.off()

# 3. Background correction and normalization using RMA
cat("Performing RMA normalization...\n")
eset <- rma(rawData)

# 4. Load and apply annotations
cat("Loading annotations...\n")
annot <- fread("data/raw/GSE49630_CEL/GPL16384_miRNA-3_1-st-v1.annotations.20140513.csv.gz", 
               skip = 4)

# Match annotations to features
matched_annot <- annot[match(featureNames(eset), annot$`Probe Set ID`), ]
featureData(eset) <- new("AnnotatedDataFrame", data = as.data.frame(matched_annot))

# 5. Filter for different probe types
cat("Filtering probe types...\n")

# Get different probe types
all_probes <- featureNames(eset)
human_mirna_mature <- grep("^hsa-miR-", all_probes, value = TRUE)
human_mirna_precursor <- grep("^hp_hsa-mir-", all_probes, value = TRUE)
control_probes <- grep("^AFFX", all_probes, value = TRUE)

cat(sprintf("Total probes: %d\n", length(all_probes)))
cat(sprintf("Human mature miRNAs: %d\n", length(human_mirna_mature)))
cat(sprintf("Human precursor miRNAs: %d\n", length(human_mirna_precursor)))
cat(sprintf("Control probes: %d\n", length(control_probes)))

# Create subset with human miRNAs only (mature + precursor)
human_mirnas <- c(human_mirna_mature, human_mirna_precursor)
eset_human <- eset[human_mirnas, ]

# 6. Post-normalization QC
cat("Generating post-normalization QC plots...\n")

# Box plot of normalized data
pdf("QC_plots/boxplot_normalized.pdf", width = 10, height = 8)
boxplot(exprs(eset_human), main = "Boxplot of normalized intensities", 
        las = 2, names = pData(eset_human)$SampleID)
dev.off()

# PCA plot
pca_data <- prcomp(t(exprs(eset_human)), scale = TRUE)
pca_df <- data.frame(
  PC1 = pca_data$x[, 1],
  PC2 = pca_data$x[, 2],
  Tissue = pData(eset_human)$Tissue
)

pdf("QC_plots/PCA_plot.pdf", width = 8, height = 6)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 3) +
  theme_bw() +
  ggtitle("PCA of normalized miRNA expression") +
  xlab(paste0("PC1 (", round(summary(pca_data)$importance[2, 1] * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca_data)$importance[2, 2] * 100, 1), "%)"))
dev.off()

# Hierarchical clustering
pdf("QC_plots/sample_clustering.pdf", width = 10, height = 8)
sample_cor <- cor(exprs(eset_human))
heatmap.2(sample_cor, 
          trace = "none", 
          col = colorRampPalette(c("blue", "white", "red"))(100),
          main = "Sample correlation heatmap",
          margins = c(10, 10))
dev.off()

# 7. Filter lowly expressed miRNAs
cat("Filtering lowly expressed miRNAs...\n")

# Calculate median expression for each miRNA
median_expr <- apply(exprs(eset_human), 1, median)

# Keep miRNAs with median expression > 4 (log2 scale)
keep <- median_expr > 4
eset_filtered <- eset_human[keep, ]

cat(sprintf("Retained %d of %d miRNAs after filtering\n", 
            sum(keep), length(keep)))

# 8. Differential expression analysis
cat("Performing differential expression analysis...\n")

# Create design matrix
tissue_type <- factor(pData(eset_filtered)$Tissue)
design <- model.matrix(~0 + tissue_type)
colnames(design) <- levels(tissue_type)

# Fit linear model
fit <- lmFit(eset_filtered, design)

# Define contrasts (all pairwise comparisons)
contrast.matrix <- makeContrasts(
  Blood_vs_Saliva = Blood - Saliva,
  Blood_vs_Semen = Blood - Semen,
  Blood_vs_Vaginal = Blood - Vaginal,
  Saliva_vs_Semen = Saliva - Semen,
  Saliva_vs_Vaginal = Saliva - Vaginal,
  Semen_vs_Vaginal = Semen - Vaginal,
  levels = design
)

# Apply contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get results for each comparison
comparisons <- colnames(contrast.matrix)
results_list <- list()

for (comp in comparisons) {
  results <- topTable(fit2, coef = comp, adjust.method = "BH", 
                     number = Inf, sort.by = "P")
  results$miRNA <- rownames(results)
  results$Comparison <- comp
  results_list[[comp]] <- results
  
  # Save individual comparison results
  write.csv(results, 
           file = paste0("results/DE_", comp, ".csv"),
           row.names = FALSE)
}

# Combine all results
all_results <- do.call(rbind, results_list)

# 9. Create visualizations for top comparison
cat("Creating visualization plots...\n")

# Volcano plot for Blood vs Saliva
results_blood_saliva <- results_list[["Blood_vs_Saliva"]]

pdf("QC_plots/volcano_Blood_vs_Saliva.pdf", width = 8, height = 6)
plot(results_blood_saliva$logFC, -log10(results_blood_saliva$P.Value),
     xlab = "log2 Fold Change", 
     ylab = "-log10(p-value)",
     main = "Volcano plot: Blood vs Saliva",
     pch = 20, col = "gray")

# Highlight significant miRNAs
sig <- results_blood_saliva$adj.P.Val < 0.05
points(results_blood_saliva$logFC[sig], 
       -log10(results_blood_saliva$P.Value)[sig],
       col = "red", pch = 20)

# Add labels for top 10 miRNAs
top10 <- head(results_blood_saliva[sig, ], 10)
text(top10$logFC, -log10(top10$P.Value), 
     labels = top10$miRNA, cex = 0.5, pos = 4)
dev.off()

# Heatmap of top 50 differentially expressed miRNAs
top_mirnas <- unique(unlist(lapply(results_list, function(x) {
  sig <- x$adj.P.Val < 0.05
  if (sum(sig) > 0) {
    head(x$miRNA[sig], 10)
  }
})))

if (length(top_mirnas) > 0) {
  pdf("QC_plots/heatmap_top_DE_miRNAs.pdf", width = 10, height = 12)
  
  # Prepare data for heatmap
  heatmap_data <- exprs(eset_filtered)[top_mirnas, ]
  
  # Scale by row (z-score)
  heatmap_data_scaled <- t(scale(t(heatmap_data)))
  
  # Create column annotations
  col_annot <- data.frame(Tissue = pData(eset_filtered)$Tissue)
  rownames(col_annot) <- colnames(heatmap_data)
  
  # Define colors
  tissue_colors <- c(Blood = "red", Saliva = "blue", 
                    Semen = "green", Vaginal = "purple")
  
  heatmap.2(heatmap_data_scaled,
            trace = "none",
            scale = "none",
            col = colorRampPalette(c("blue", "white", "red"))(100),
            ColSideColors = tissue_colors[col_annot$Tissue],
            margins = c(10, 10),
            main = "Top differentially expressed miRNAs",
            key.title = "Z-score")
  
  # Add legend
  legend("topright", legend = names(tissue_colors), 
         fill = tissue_colors, title = "Tissue", bty = "n")
  
  dev.off()
}

# 10. Export normalized expression matrix
cat("Exporting results...\n")
dir.create("results", showWarnings = FALSE)

# Export normalized expression values
write.csv(exprs(eset_filtered), 
         file = "results/normalized_expression_filtered.csv")

# Export sample information
write.csv(pData(eset_filtered), 
         file = "results/sample_info.csv",
         row.names = FALSE)

# Summary statistics
summary_stats <- data.frame(
  Comparison = names(results_list),
  Total_Tested = sapply(results_list, nrow),
  Significant_p05 = sapply(results_list, function(x) sum(x$P.Value < 0.05)),
  Significant_FDR05 = sapply(results_list, function(x) sum(x$adj.P.Val < 0.05)),
  Significant_FDR01 = sapply(results_list, function(x) sum(x$adj.P.Val < 0.01))
)

write.csv(summary_stats, 
         file = "results/DE_summary_statistics.csv",
         row.names = FALSE)

cat("\nAnalysis complete!\n")
cat("Results saved in 'results' directory\n")
cat("QC plots saved in 'QC_plots' directory\n")

# Print summary
print(summary_stats)

# Session info
sink("results/session_info.txt")
sessionInfo()
sink()