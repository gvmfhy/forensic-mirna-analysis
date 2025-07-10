# Manual CEL Processing Instructions

These Affymetrix miRNA ST arrays require the 'oligo' package.
Since automatic installation is failing, here are manual steps:

## Option 1: Use a different R environment
1. Install R from CRAN (not homebrew)
2. Install BiocManager: install.packages('BiocManager')
3. Install oligo: BiocManager::install('oligo')
4. Run: Rscript scripts/preprocessing/cel_processor_minimal.R

## Option 2: Use Docker
```bash
docker run -it -v $(pwd):/data bioconductor/bioconductor_docker:latest R
# Inside container:
BiocManager::install('oligo')
source('/data/scripts/preprocessing/cel_processor_minimal.R')
```

## Option 3: Process on another machine
Transfer the CEL files to a machine with working Bioconductor

## What the data contains:
-  20 CEL files ( 5 blood, 5 saliva, 5 semen, 5 vaginal)
- Platform: Affymetrix miRNA 3.0/4.0 ST arrays
- These are individual samples (not pooled like GPR)

