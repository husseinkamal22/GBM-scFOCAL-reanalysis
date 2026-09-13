# ============================================
# GBM scFOCAL Reanalysis - Final Pipeline
# ============================================

library(Seurat)
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)

# --- Paths ---
base_dir <- here::here()
data_dir <- file.path(base_dir, "data")
output_dir <- file.path(base_dir, "output")
fig_dir <- file.path(base_dir, "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================
# 1. Load Data
# ============================================
cat("Loading data...\n")
counts <- fread(file.path(data_dir, "GSE229779_countsMatrix.tsv"))
counts_mat <- as.matrix(counts[, -1])
rownames(counts_mat) <- counts$V1
counts_sparse <- as(counts_mat, "dgCMatrix")
rm(counts_mat)

meta <- fread(file.path(data_dir, "GSE229779_cellMetadata.tsv"))
rownames(meta) <- meta$cellID

# ============================================
# 2. Create Seurat Object
# ============================================
cat("Creating Seurat object...\n")
sc <- CreateSeuratObject(
  counts = counts_sparse,
  meta.data = meta,
  project = "GBM_reanalysis",
  min.cells = 3,
  min.features = 200
)
rm(counts_sparse)

# ============================================
# 3. QC & Filtering
# ============================================
cat("QC & filtering...\n")
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

p_qc <- VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(fig_dir, "01_qc_violin.png"), p_qc, width = 12, height = 4)

cat("  Before:", ncol(sc), "cells\n")
sc <- subset(sc, nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 10)
cat("  After:", ncol(sc), "cells\n")

# ============================================
# 4. Normalization & Scaling
# ============================================
cat("Normalizing...\n")
sc <- NormalizeData(sc)
sc <- FindVariableFeatures(sc, nfeatures = 2000)
sc <- ScaleData(sc)

# ============================================
# 5. Dimensionality Reduction
# ============================================
cat("Running PCA...\n")
sc <- RunPCA(sc, npcs = 50)

p_elbow <- ElbowPlot(sc, ndims = 50)
ggsave(file.path(fig_dir, "02_elbow_plot.png"), p_elbow, width = 6, height = 4)

# ============================================
# 6. Clustering & UMAP
# ============================================
cat("Clustering...\n")
sc <- FindNeighbors(sc, dims = 1:12)
sc <- FindClusters(sc, resolution = 0.5)
sc <- RunUMAP(sc, dims = 1:12)

p_umap_clusters <- DimPlot(sc, group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("Clusters")
ggsave(file.path(fig_dir, "03_umap_clusters.png"), p_umap_clusters, width = 10, height = 8)

# ============================================
# 7. Cell State Annotation
# ============================================
cat("Annotating cell states...\n")
canonical_markers <- list(
  "AC-like" = c("GFAP", "AQP4", "S100B"),
  "NPC-like" = c("SOX2", "NES", "VIM"),
  "OPC-like" = c("OLIG2", "PDGFRA", "CSPG4"),
  "MES-like" = c("CHI3L1", "CD44", "ANGPTL4"),
  "Oligodendrocytes" = c("MBP", "PLP1", "MOG"),
  "Myeloid" = c("CSF1R", "CD68", "C1QB"),
  "T-Cells" = c("CD3D", "CD3E", "IL7R")
)

p_dot <- DotPlot(sc, features = unlist(canonical_markers), group.by = "seurat_clusters") +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "04_dotplot_markers.png"), p_dot, width = 14, height = 6)

sc$cell_state <- NA
sc$cell_state[sc$seurat_clusters %in% c(12)] <- "NPC-like"
sc$cell_state[sc$seurat_clusters %in% c(4, 8, 13, 16)] <- "AC-like"
sc$cell_state[sc$seurat_clusters %in% c(10, 11, 19)] <- "Myeloid"
sc$cell_state[sc$seurat_clusters %in% c(3, 6, 7, 14, 17, 18, 23)] <- "OPC-like"
sc$cell_state[sc$seurat_clusters %in% c(0, 1, 5, 9, 21)] <- "MES-like"
sc$cell_state[sc$seurat_clusters %in% c(22)] <- "Oligodendrocytes"
sc$cell_state[sc$seurat_clusters %in% c(20)] <- "T-Cells"
sc$cell_state[sc$seurat_clusters %in% c(15)] <- "Immune"
sc$cell_state[sc$seurat_clusters %in% c(2, 24)] <- "Unknown"

p_umap_states <- DimPlot(sc, group.by = "cell_state", label = TRUE, repel = TRUE) +
  ggtitle("Cell States")
ggsave(file.path(fig_dir, "05_umap_cell_states.png"), p_umap_states, width = 10, height = 8)

# ============================================
# 8. Find Markers (uncomment when ready — slow on large datasets)
# ============================================
# cat("Finding markers...\n")
# markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#
# top5 <- markers %>%
#   group_by(cluster) %>%
#   slice_max(avg_log2FC, n = 5)
#
# write.csv(top5, file.path(output_dir, "top5_markers.csv"), row.names = FALSE)
# write.csv(markers, file.path(output_dir, "all_markers.csv"), row.names = FALSE)
saveRDS(sc, file.path(output_dir, "gbm_seurat_object.rds"))

cat("Done!\n")
