#!/usr/bin/env Rscript
# 05_marker_identification.R
# Identify cluster-specific marker genes

library(Seurat)
library(ggplot2)
library(dplyr)

# ---- Load Data ----
seurat_obj <- readRDS("data/processed/seurat_object_clustered.rds")

# ---- Find All Markers ----
message("Finding cluster markers with FindAllMarkers...")
cluster_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,        # Only positive markers
  min.pct = 0.25,         # Expressed in at least 25% of cells in cluster
  logfc.threshold = 0.25  # Minimum log2 fold change
)

# Save marker table
write.csv(cluster_markers, "results/cluster_markers_all.csv", row.names = FALSE)

# ---- Top Markers per Cluster ----
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(top_markers, "results/top10_markers_per_cluster.csv", row.names = FALSE)

# ---- Highlight Known Markers ----
# Genes mentioned in your README: S100A8, HBB, IFI27
known_markers <- c("S100A8", "S100A9", "HBB", "IFI27", "IL6", "TNF", "IFNG")

# Check which markers are present in the dataset
available_markers <- known_markers[known_markers %in% rownames(seurat_obj)]
message("Plotting known markers: ", paste(available_markers, collapse = ", "))

if (length(available_markers) > 0) {
  feature_plot <- FeaturePlot(
    seurat_obj,
    features = available_markers,
    ncol = 3,
    pt.size = 0.5,
    order = TRUE
  )
  ggsave("results/known_markers_featureplot.pdf", feature_plot, width = 12, height = 10)
  ggsave("results/known_markers_featureplot.png", feature_plot, width = 12, height = 10, dpi = 300)
}

# ---- Heatmap of Top Markers ----
top10_for_heatmap <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

heatmap_plot <- DoHeatmap(
  seurat_obj,
  features = top10_for_heatmap,
  group.by = "ident",
  group.bar = TRUE,
  size = 3
) + NoLegend()

ggsave("results/cluster_marker_heatmap.pdf", heatmap_plot, width = 14, height = 12)
ggsave("results/cluster_marker_heatmap.png", heatmap_plot, width = 14, height = 12, dpi = 300)

# ---- Violin Plots for Key Genes ----
vln_plot <- VlnPlot(
  seurat_obj,
  features = available_markers,
  ncol = 3,
  pt.size = 0
)
ggsave("results/marker_violin_plots.pdf", vln_plot, width = 12, height = 10)

# ---- Save ----
saveRDS(cluster_markers, "data/processed/cluster_markers.rds")
message("Marker identification complete. Results saved to results/")
