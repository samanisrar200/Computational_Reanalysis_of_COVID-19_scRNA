#!/usr/bin/env Rscript
# 03_normalization.R
# Log normalization and feature selection

library(Seurat)

# ---- Load Data ----
seurat_obj <- readRDS("data/processed/seurat_object_qc.rds")

# ---- Normalization ----
message("Performing log normalization...")
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000  # Per 10,000 (CPM-like)
)

# ---- Feature Selection ----
message("Finding variable features...")
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",  # Variance stabilizing transformation
  nfeatures = 2000           # Top 2000 variable genes
)

# Visualize top variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

ggsave("results/variable_features.pdf", plot2, width = 10, height = 6)
ggsave("results/variable_features.png", plot2, width = 10, height = 6, dpi = 300)

message("Top 10 variable features: ", paste(top10, collapse = ", "))

# ---- Scaling ----
message("Scaling data...")
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

# ---- Save ----
saveRDS(seurat_obj, "data/processed/seurat_object_normalized.rds")
message("Normalization complete.")
