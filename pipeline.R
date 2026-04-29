#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

set.seed(123)

# -----------------------------
# CONFIG
# -----------------------------
RDS_FILE <- "data/GSM4557327_555_1_cell.counts.matrices.rds"
OUT_DIR <- "results"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# LOAD DATA (FIXED)
# -----------------------------
message("Loading RDS file...")

if (!file.exists(RDS_FILE)) {
  stop("File not found: ", RDS_FILE)
}

raw_data <- readRDS(RDS_FILE)

message("Inspecting structure...")
print(str(raw_data, max.level = 1))

# -----------------------------
# EXTRACT COUNTS (ROBUST)
# -----------------------------
if (is.list(raw_data)) {
  
  slots <- names(raw_data)
  message("Detected slots: ", paste(slots, collapse = ", "))
  
  if ("spliced" %in% slots) {
    counts <- raw_data$spliced
    message("Using: spliced matrix")
  } else if ("exon" %in% slots) {
    counts <- raw_data$exon
    message("Using: exon matrix")
  } else if ("counts" %in% slots) {
    counts <- raw_data$counts
    message("Using: counts matrix")
  } else {
    stop("Unknown structure. Inspect raw_data manually.")
  }
  
} else {
  counts <- raw_data
  message("Using raw matrix directly")
}

counts <- as(counts, "dgCMatrix")

# Check names
if (is.null(rownames(counts)) || is.null(colnames(counts))) {
  stop("Counts matrix must have gene names and cell barcodes")
}

message("Matrix dimensions: ", nrow(counts), " genes x ", ncol(counts), " cells")

# -----------------------------
# CREATE SEURAT OBJECT
# -----------------------------
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "COVID19_scRNA",
  min.cells = 3,
  min.features = 200
)

seurat_obj$sample <- "GSM4557327"
seurat_obj$condition <- "COVID19"

# -----------------------------
# QC METRICS
# -----------------------------
message("Calculating QC metrics...")

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

pdf(file.path(OUT_DIR, "qc_violin.pdf"))
print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")))
dev.off()

# -----------------------------
# FILTER CELLS
# -----------------------------
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 7500 &
    nCount_RNA > 500 &
    nCount_RNA < 50000 &
    percent.mt < 15
)

message("Cells after QC: ", ncol(seurat_obj))

# -----------------------------
# NORMALIZATION
# -----------------------------
message("Normalizing...")

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)

# -----------------------------
# PCA
# -----------------------------
message("Running PCA...")

seurat_obj <- RunPCA(seurat_obj, npcs = 50)

pdf(file.path(OUT_DIR, "pca_elbow.pdf"))
print(ElbowPlot(seurat_obj))
dev.off()

NDIMS <- 25

# -----------------------------
# CLUSTERING
# -----------------------------
message("Clustering...")

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:NDIMS)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)

# -----------------------------
# UMAP
# -----------------------------
message("Running UMAP...")

seurat_obj <- RunUMAP(seurat_obj, dims = 1:NDIMS, seed.use = 123)

pdf(file.path(OUT_DIR, "umap.pdf"))
print(DimPlot(seurat_obj, label = TRUE))
dev.off()

# -----------------------------
# MARKER GENES
# -----------------------------
message("Finding markers...")

cluster_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(cluster_markers, file.path(OUT_DIR, "cluster_markers.csv"), row.names = FALSE)

top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(top_markers, file.path(OUT_DIR, "top_markers.csv"), row.names = FALSE)

# -----------------------------
# FEATURE PLOTS
# -----------------------------
known_genes <- c("S100A8", "S100A9", "HBB", "IFI27", "IL6", "TNF", "IFNG")
available <- known_genes[known_genes %in% rownames(seurat_obj)]

if (length(available) > 0) {
  pdf(file.path(OUT_DIR, "known_markers.pdf"))
  print(FeaturePlot(seurat_obj, features = available))
  dev.off()
}

# -----------------------------
# GO ENRICHMENT
# -----------------------------
message("Running GO enrichment...")

background <- rownames(seurat_obj)
bg_entrez <- bitr(background, "SYMBOL", "ENTREZID", org.Hs.eg.db)

clusters <- unique(cluster_markers$cluster)

for (cl in clusters) {
  
  genes <- cluster_markers %>%
    filter(cluster == cl) %>%
    slice_max(order_by = avg_log2FC, n = 100) %>%
    pull(gene)
  
  gene_entrez <- bitr(genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)
  
  if (nrow(gene_entrez) < 5) next
  
  ego <- enrichGO(
    gene = gene_entrez$ENTREZID,
    universe = bg_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    
    write.csv(as.data.frame(ego),
              file.path(OUT_DIR, paste0("go_cluster_", cl, ".csv")),
              row.names = FALSE)
    
    pdf(file.path(OUT_DIR, paste0("go_dotplot_", cl, ".pdf")))
    print(dotplot(ego, showCategory = 15))
    dev.off()
  }
}

# -----------------------------
# SAVE FINAL OBJECT
# -----------------------------
saveRDS(seurat_obj, file.path(OUT_DIR, "final_seurat.rds"))

message("Pipeline completed successfully ")

