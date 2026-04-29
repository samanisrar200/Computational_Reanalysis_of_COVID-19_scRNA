#!/usr/bin/env Rscript
# 06_functional_enrichment.R
# GO Biological Process enrichment of marker genes

library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(enrichplot)
library(ggplot2)
library(dplyr)


cluster_markers <- readRDS("data/processed/cluster_markers.rds")

# Convert gene symbols to Entrez IDs for clusterProfiler
clusters <- unique(cluster_markers$cluster)


seurat_obj <- readRDS("data/processed/seurat_object_clustered.rds")
background_genes <- rownames(seurat_obj)


bg_entrez <- bitr(
  background_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)


for (cl in clusters) {
  message("Processing cluster ", cl, "...")
  
  # Get top markers for this cluster (e.g., top 100 by logFC)
  cluster_genes <- cluster_markers %>%
    filter(cluster == cl) %>%
    top_n(n = 100, wt = avg_log2FC) %>%
    pull(gene)
  
  if (length(cluster_genes) < 10) {
    message("  Too few markers, skipping...")
    next
  }
  

  gene_entrez <- bitr(
    cluster_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  if (nrow(gene_entrez) < 5) {
    message("  Too few mapped genes, skipping...")
    next
  }
  

  ego <- enrichGO(
    gene = gene_entrez$ENTREZID,
    universe = bg_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",              
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE          
  )
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    # Save results
    write.csv(
      as.data.frame(ego),
      paste0("results/go_bp_cluster_", cl, ".csv"),
      row.names = FALSE
    )
    
    # Dotplot
    dotplot <- dotplot(ego, showCategory = 15) +
      ggtitle(paste("GO BP Enrichment - Cluster", cl))
    
    ggsave(
      paste0("results/go_dotplot_cluster_", cl, ".pdf"),
      dotplot,
      width = 10,
      height = 8
    )
    ggsave(
      paste0("results/go_dotplot_cluster_", cl, ".png"),
      dotplot,
      width = 10,
      height = 8,
      dpi = 300
    )
    
    # Barplot
    barplot <- barplot(ego, showCategory = 15) +
      ggtitle(paste("GO BP Enrichment - Cluster", cl))
    
    ggsave(
      paste0("results/go_barplot_cluster_", cl, ".pdf"),
      barplot,
      width = 10,
      height = 8
    )
  }
}

# ---- Combined GO Analysis (All Significant Markers) ----
message("Running combined GO analysis...")

all_sig_markers <- cluster_markers %>%
  filter(p_val_adj < 0.05) %>%
  pull(gene) %>%
  unique()

all_entrez <- bitr(
  all_sig_markers,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

ego_combined <- enrichGO(
  gene = all_entrez$ENTREZID,
  universe = bg_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

if (!is.null(ego_combined) && nrow(as.data.frame(ego_combined)) > 0) {
  write.csv(as.data.frame(ego_combined), "results/go_bp_all_clusters.csv", row.names = FALSE)
  
  # Combined dotplot
  combined_dot <- dotplot(ego_combined, showCategory = 20) +
    ggtitle("GO BP Enrichment - All Clusters")
  
  ggsave("results/go_dotplot_all_clusters.pdf", combined_dot, width = 10, height = 10)
  ggsave("results/go_dotplot_all_clusters.png", combined_dot, width = 10, height = 10, dpi = 300)
  
  # Network plot
  ego_pairwise <- pairwise_termsim(ego_combined)
  emap <- emapplot(ego_pairwise, showCategory = 30)
  ggsave("results/go_emap_all_clusters.pdf", emap, width = 12, height = 10)
}

message("Functional enrichment analysis complete.")
