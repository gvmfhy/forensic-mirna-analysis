#!/bin/bash
# Script to clean up redundant preprocessing files

echo "Creating archive directory..."
mkdir -p scripts/preprocessing/archive

echo "Moving abandoned files to archive..."
mv scripts/preprocessing/cel_processor.R scripts/preprocessing/archive/ 2>/dev/null || true
mv scripts/preprocessing/cel_processor_simple.py scripts/preprocessing/archive/ 2>/dev/null || true
mv scripts/preprocessing/cel_processor_rpy2.py scripts/preprocessing/archive/ 2>/dev/null || true
mv scripts/preprocessing/cel_processor_final.py scripts/preprocessing/archive/ 2>/dev/null || true
mv minimal_cel_process.R scripts/preprocessing/archive/ 2>/dev/null || true
mv process_affymetrix_mirna.R scripts/preprocessing/archive/ 2>/dev/null || true

echo "Cleanup complete! Abandoned files moved to scripts/preprocessing/archive/"
