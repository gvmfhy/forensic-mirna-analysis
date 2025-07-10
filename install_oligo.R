#!/usr/bin/env Rscript

# Simple script to install just oligo package
cat("Installing oligo package...\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

BiocManager::install("oligo", ask = FALSE, update = FALSE)

# Test if installed
if ("oligo" %in% installed.packages()[,"Package"]) {
    cat("\nSUCCESS: oligo is installed!\n")
    library(oligo)
    cat("oligo version:", packageVersion("oligo"), "\n")
} else {
    cat("\nERROR: oligo installation failed\n")
    quit(status = 1)
}