#!/usr/bin/env Rscript
# 04_dimensionality_reduction.R
# PCA, UMAP, and clustering

library(Seurat)
library(ggplot2)

# ---- Load Data ----
seurat_obj <- readRDS("data/processed/seurat_object_normalized.rds")

# ---- PCA ----
message("Running PCA...")
seurat_obj <- RunPCA(
  seurat_obj,
  features = VariableFeatures(object = seurat_obj),
  npcs = 50
)

# Visualize PCA results
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)

# Elbow plot to determine significant PCs
elbow_plot <- ElbowPlot(seurat_obj, ndims = 50)
ggsave("results/pca_elbow_plot.pdf", elbow_plot, width = 8, height = 6)
ggsave("results/pca_elbow_plot.png", elbow_plot, width = 8, height = 6, dpi = 300)

# PCA heatmap (first 15 PCs)
pca_heatmap <- DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)
ggsave("results/pca_heatmap.pdf", pca_heatmap, width = 12, height = 10)


# Based on elbow plot, typically use 20-30 PCs for COVID data
NDIMS <- 25  # Adjust based on elbow plot

# ---- Clustering ----
message("Finding neighbors and clustering...")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:NDIMS)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)  # Adjust resolution as needed

message("Identified clusters: ", length(unique(Idents(seurat_obj))))

# ---- UMAP ----
message("Running UMAP...")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:NDIMS, n.neighbors = 30, min.dist = 0.3)

# Visualize clusters
umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5)
ggsave("results/umap_clusters.pdf", umap_plot, width = 8, height = 7)
ggsave("results/umap_clusters.png", umap_plot, width = 8, height = 7, dpi = 300)

# UMAP without labels
umap_no_label <- DimPlot(seurat_obj, reduction = "umap", label = FALSE, pt.size = 0.5)
ggsave("results/umap_clusters_no_label.pdf", umap_no_label, width = 8, height = 7)

# ---- Save ----
saveRDS(seurat_obj, "data/processed/seurat_object_clustered.rds")
message("Dimensionality reduction and clustering complete.")
