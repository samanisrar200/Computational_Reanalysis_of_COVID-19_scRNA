# Single-cell RNA-seq analysis of a COVID-19 blood sample (GSM4557327)
# Goal: find out which types of cells are in the sample

library(Seurat)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

set.seed(123)

data_file <- data_file <- "C:/Users/saman/Desktop/independant projects/scRNA-COVID19/data/GSM4557327_555_1_cell.counts.matrices.rds"
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)


# STEP 1: Load the data
# Each number = how many times that gene was detected in that cell.
raw_data <- readRDS(data_file)

counts <- raw_data$spliced
print(dim(counts))   # genes x cells


# STEP 2: Make a Seurat object
#min.cells = 3 : ignore genes found in fewer than 3 cells
# min.features = 200 : ignore cells with fewer than 200 genes
counts <- raw_data$exon
print(class(counts))
print(dim(counts))
seurat_obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
names(raw_data)        
class(counts)          
dim(counts)           
head(rownames(counts))
head(colnames(counts))

# STEP 3: Quality control
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

pdf(file.path(out_dir, "qc_violin.pdf"), width = 10, height = 6)
print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")))
dev.off()

# Remove bad cells:
seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA > 200 & nFeature_RNA < 7500 &
                       nCount_RNA > 500 & nCount_RNA < 50000 &
                       percent.mt < 15)
print(ncol(seurat_obj))   # cells left after QC


# STEP 4: Normalize
seurat_obj <- NormalizeData(seurat_obj)

# Keep the 2000 genes that differ most between cells
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)

# Scale genes so highly expressed genes don't dominate
seurat_obj <- ScaleData(seurat_obj)


# STEP 5: PCA
# Shrinks 2000 genes down to principal components
seurat_obj <- RunPCA(seurat_obj, npcs = 50)

# The elbow plot helps choose how many PCs to keep
pdf(file.path(out_dir, "pca_elbow.pdf"), width = 8, height = 6)
print(ElbowPlot(seurat_obj, ndims = 50))
dev.off()

num_pcs <- 25   


# STEP 6: Clustering and UMAP
# Group cells that look alike
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:num_pcs)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)

# UMAP = 2D picture of the cells, similar cells are close together
seurat_obj <- RunUMAP(seurat_obj, dims = 1:num_pcs)

pdf(file.path(out_dir, "umap_clusters.pdf"), width = 10, height = 8)
print(DimPlot(seurat_obj, label = TRUE))
dev.off()


# STEP 7: Marker genes
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file.path(out_dir, "cluster_markers.csv"), row.names = FALSE)

# Top 10 per cluster
top_markers <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10)
write.csv(top_markers, file.path(out_dir, "top_markers.csv"), row.names = FALSE)


# STEP 8: Check known cell-type markers
known_genes <- c("CD14", "LYZ",       # monocytes
                 "CD3D", "IL7R",      # T cells
                 "MS4A1", "CD79A",    # B cells
                 "GNLY", "NKG7",      # NK cells
                 "IFI27", "IFITM3")   # interferon response

pdf(file.path(out_dir, "known_markers.pdf"), width = 12, height = 10)
print(FeaturePlot(seurat_obj, features = known_genes, ncol = 3))
dev.off()


# STEP 9: Cluster sizes
cluster_sizes <- as.data.frame(table(Idents(seurat_obj)))
colnames(cluster_sizes) <- c("Cluster", "Cells")
cluster_sizes$Percent <- round(100 * cluster_sizes$Cells / sum(cluster_sizes$Cells), 1)
write.csv(cluster_sizes, file.path(out_dir, "cluster_sizes.csv"), row.names = FALSE)

pdf(file.path(out_dir, "cluster_sizes.pdf"), width = 8, height = 5)
print(ggplot(cluster_sizes, aes(x = Cluster, y = Cells)) +
        geom_col() +
        labs(title = "Number of cells in each cluster"))
dev.off()


# STEP 10: GO enrichment
# First convert gene names to ENTREZ IDs (what the GO database uses)
all_genes <- bitr(rownames(seurat_obj), fromType = "SYMBOL",
                  toType = "ENTREZID", OrgDb = org.Hs.eg.db)

for (cl in unique(markers$cluster)) {
  
  # top 100 marker genes of this cluster
  cluster_genes <- markers %>%
    filter(cluster == cl) %>%
    slice_max(avg_log2FC, n = 100) %>%
    pull(gene)
  
  gene_ids <- bitr(cluster_genes, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  go_result <- enrichGO(gene = gene_ids$ENTREZID,
                        universe = all_genes$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        ont = "BP",              # biological process
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05)
  
  # only save if something significant was found
  if (!is.null(go_result) && nrow(as.data.frame(go_result)) > 0) {
    write.csv(as.data.frame(go_result),
              file.path(out_dir, paste0("go_cluster_", cl, ".csv")), row.names = FALSE)
    
    pdf(file.path(out_dir, paste0("go_dotplot_", cl, ".pdf")), width = 10, height = 8)
    print(dotplot(go_result, showCategory = 15))
    dev.off()
  }
}


pdf(file.path(out_dir, "final_known_markers_dotplot.pdf"), width = 10, height = 6)
print(DotPlot(seurat_obj, features = known_genes) + RotatedAxis())
dev.off()


top3 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 3)
pdf(file.path(out_dir, "final_top_markers_dotplot.pdf"), width = 14, height = 6)
print(DotPlot(seurat_obj, features = unique(top3$gene)) + RotatedAxis())
dev.off()


summary_table <- cluster_sizes
summary_table$Top_genes <- NA
summary_table$Top_GO_term <- NA

for (i in 1:nrow(summary_table)) {
  cl <- as.character(summary_table$Cluster[i])
  
  top5 <- markers %>%
    filter(cluster == cl) %>%
    slice_max(avg_log2FC, n = 5) %>%
    pull(gene)
  summary_table$Top_genes[i] <- paste(top5, collapse = ", ")
  

  go_file <- file.path(out_dir, paste0("go_cluster_", cl, ".csv"))
  if (file.exists(go_file)) {
    go_table <- read.csv(go_file)
    summary_table$Top_GO_term[i] <- go_table$Description[1]
  } else {
    summary_table$Top_GO_term[i] <- "none significant"
  }
}

write.csv(summary_table, file.path(out_dir, "final_summary_table.csv"), row.names = FALSE)
print(summary_table)

print(paste("Cells after QC:", ncol(seurat_obj)))
print(paste("Number of clusters:", nrow(cluster_sizes)))
print(paste("Largest cluster:", cluster_sizes$Cluster[which.max(cluster_sizes$Cells)],
            "with", max(cluster_sizes$Percent), "% of cells"))


# STEP 11: Save the final object
saveRDS(seurat_obj, file.path(out_dir, "final_seurat.rds"))
print("Done")

