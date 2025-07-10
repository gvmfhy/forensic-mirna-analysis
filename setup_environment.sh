#!/bin/bash

# Environment Setup Script for Forensic miRNA Analysis

echo "=== Setting up project environments ==="

# Python virtual environment
echo "1. Creating Python virtual environment..."
python3 -m venv venv
source venv/bin/activate

echo "2. Installing Python dependencies..."
pip install --upgrade pip
pip install numpy pandas scikit-learn matplotlib seaborn
pip install lightgbm xgboost shap
pip install jupyterlab scipy statsmodels

# Save Python requirements
pip freeze > requirements.txt

echo "3. Python environment ready!"

# R environment setup
echo "4. Creating R environment file..."
cat > setup_r_env.R << 'EOF'
# R Environment Setup

# Check if renv is installed
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Initialize renv for this project
renv::init()

# Install required packages
packages <- c(
  # Data manipulation
  "tidyverse",
  "data.table",
  
  # Bioconductor packages
  "BiocManager"
)

# Install CRAN packages
install.packages(packages)

# Install Bioconductor packages
BiocManager::install(c(
  "affy",
  "oligo", 
  "limma",
  "sva",
  "GEOquery",
  "Biobase"
))

# Machine learning packages
install.packages(c(
  "caret",
  "randomForest",
  "glmnet",
  "e1071"
))

# Save the environment
renv::snapshot()

cat("R environment setup complete!\n")
EOF

echo "5. Project structure verification..."
cat > project_info.txt << 'EOF'
Project: Forensic miRNA Cross-Platform Analysis
Python: venv (virtual environment)
R: renv (reproducible environment)

Key directories:
- data/: Raw and processed data
- scripts/: Analysis scripts
- results/: Output files
- reports: Final reports

To activate environments:
Python: source venv/bin/activate
R: renv::restore() in R session
EOF

echo "=== Setup complete! ==="
echo ""
echo "To activate Python environment: source venv/bin/activate"
echo "To setup R environment: Rscript setup_r_env.R"