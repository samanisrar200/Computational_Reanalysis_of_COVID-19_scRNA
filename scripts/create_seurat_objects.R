#!/usr/bin/env Rscript
# 01_create_seurat_object.R
# Create Seurat object from GEO GSM4557327 raw counts

library(Seurat)
library(Matrix)

# ---- Configuration ----
DATA_DIR <- "data/"
RDS_FILE <- file.path(DATA_DIR, "GSM4557327.rds") 
OUT_DIR <- "data/processed/"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Load Raw Data ----
message("Loading raw count matrix...")
raw_data <- readRDS(RDS_FILE)

# Extract exon counts (most common for scRNA-seq analysis)
# Adjust based on actual structure of your RDS file
if (is.list(raw_data)) {
  counts <-  raw_data[["exon"]]
} else {
  counts <- raw_data 
}

message("Matrix dimensions: ", nrow(counts), " genes x ", ncol(counts), " cells")

# ---- Create Seurat Object ----
# Initialize with raw counts (Seurat expects genes as rows, cells as columns)
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "COVID19_scRNA",
  min.cells = 3,      # Keep genes expressed in at least 3 cells
  min.features = 200  # Keep cells with at least 200 genes detected
)

# Add metadata
seurat_obj$sample <- "GSM4557327"
seurat_obj$condition <- "COVID19"

message("Initial Seurat object:")
print(seurat_obj)

# ---- Save ----
saveRDS(seurat_obj, file.path(OUT_DIR, "seurat_object_raw.rds"))
message("Saved: ", file.path(OUT_DIR, "seurat_object_raw.rds"))
