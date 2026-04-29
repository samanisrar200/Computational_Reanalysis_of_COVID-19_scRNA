#!/bin/bash
# run_pipeline.sh
# Master script to run the full scRNA-seq analysis pipeline

set -e  # Exit on error

echo "=========================================="
echo "COVID-19 scRNA-seq Reanalysis Pipeline"
echo "=========================================="

# Create directories
mkdir -p data/raw data/processed results

# 1: Create Seurat object
echo "[1/7] Creating Seurat object..."
Rscript scripts/01_create_seurat_object.R

# 2: QC filtering
echo "[2/7] Quality control..."
Rscript scripts/02_qc_filtering.R

# 3: Normalization
echo "[3/7] Normalization and feature selection..."
Rscript scripts/03_normalization.R

# 4: Dimensionality reduction and clustering
echo "[4/7] Dimensionality reduction and clustering..."
Rscript scripts/04_dimensionality_reduction.R

# 5: Marker identification
echo "[5/7] Identifying cluster markers..."
Rscript scripts/05_marker_identification.R

# 6: Functional enrichment
echo "[6/7] GO enrichment analysis..."
Rscript scripts/06_functional_enrichment.R


echo "=========================================="
echo "Pipeline complete! Results in results/"
echo "=========================================="
