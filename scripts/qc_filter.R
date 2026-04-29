#!/usr/bin/env Rscript
# 02_qc_filtering.R
# Quality control: filter low-quality cells and calculate mitochondrial percentage

library(Seurat)
library(ggplot2)

# ---- Load Data ----
seurat_obj <- readRDS("data/processed/seurat_object_raw.rds")

# ---- Calculate QC Metrics ----
# Identify mitochondrial genes (human: MT- prefix)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# View QC metrics
head(seurat_obj@meta.data, 5)

# ---- Visualize QC Distributions ----
qc_plots <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)

ggsave("results/qc_violin_plots.pdf", qc_plots, width = 12, height = 5)
ggsave("results/qc_violin_plots.png", qc_plots, width = 12, height = 5, dpi = 300)

# Scatter plots for feature relationships
scatter_plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
scatter_plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
combined_scatter <- scatter_plot1 + scatter_plot2

ggsave("results/qc_scatter_plots.pdf", combined_scatter, width = 10, height = 5)
ggsave("results/qc_scatter_plots.png", combined_scatter, width = 10, height = 5, dpi = 300)

# ---- Filter Cells ----
# Adjust thresholds based on your data distribution
# Typical thresholds for 10x Genomics data:
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 200 & 
           nFeature_RNA < 7500 &    # Remove doublets/high outliers
           nCount_RNA > 500 & 
           nCount_RNA < 50000 &
           percent.mt < 15            # Remove stressed/dead cells
)

message("Cells after QC filtering: ", ncol(seurat_obj))

# ---- Save ----
saveRDS(seurat_obj, "data/processed/seurat_object_qc.rds")
message("QC complete. Saved filtered object.")
