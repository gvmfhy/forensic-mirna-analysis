
# Simplified CEL processing script
# Run this in R after installing: BiocManager::install(c("affy", "oligo"))

library(oligo)

# File paths
cel_files <- c("data/raw/GSE49630_CEL/GSM1202839_Vaginal_02.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202836_Semen_04.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202833_Semen_01.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202841_Vaginal_04.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202831_Blood_04.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202846_Saliva_04.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202845_Saliva_03.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202835_Semen_03.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202829_Blood_02.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202832_Blood_05.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202838_Vaginal_01.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202844_Saliva_02.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202834_Semen_02.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202847_Saliva_05.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202828_Blood_01.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202843_Saliva_01.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202830_Blood_03.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202842_Vaginal_05.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202837_Semen_05.CEL.gz",
  "data/raw/GSE49630_CEL/GSM1202840_Vaginal_03.CEL.gz")
sample_names <- c("GSM1202839_Vaginal_02",
  "GSM1202836_Semen_04",
  "GSM1202833_Semen_01",
  "GSM1202841_Vaginal_04",
  "GSM1202831_Blood_04",
  "GSM1202846_Saliva_04",
  "GSM1202845_Saliva_03",
  "GSM1202835_Semen_03",
  "GSM1202829_Blood_02",
  "GSM1202832_Blood_05",
  "GSM1202838_Vaginal_01",
  "GSM1202844_Saliva_02",
  "GSM1202834_Semen_02",
  "GSM1202847_Saliva_05",
  "GSM1202828_Blood_01",
  "GSM1202843_Saliva_01",
  "GSM1202830_Blood_03",
  "GSM1202842_Vaginal_05",
  "GSM1202837_Semen_05",
  "GSM1202840_Vaginal_03")

# Read CEL files
cat("Reading CEL files...\n")
raw_data <- read.celfiles(cel_files)

# RMA normalization
cat("Performing RMA normalization...\n")
eset <- rma(raw_data)

# Extract expression matrix
expr_matrix <- exprs(eset)
colnames(expr_matrix) <- sample_names

# Save as CSV
write.csv(expr_matrix, "data/processed/cel/cel_raw_expression.csv")
cat("Saved expression matrix to data/processed/cel/cel_raw_expression.csv\n")

# Print summary
cat(sprintf("Processed %d probes x %d samples\n", nrow(expr_matrix), ncol(expr_matrix)))
