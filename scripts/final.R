# ============================================
# GBM scFOCAL Reanalysis - Final Pipeline
# ============================================

library(Seurat)
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scDblFinder)
library(SingleR)
library(BiocSingular)

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
# 3. QC
# ============================================
cat("QC...\n")
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

p_qc <- VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(fig_dir, "01_qc_violin.png"), p_qc, width = 12, height = 4)

# ============================================
# 4. Doublet Detection (scDblFinder)
# ============================================
cat("Detecting doublets...\n")
sc <- scDblFinder(sc, dbr = 0.01)

p_dbl <- DimPlot(sc, group.by = "scDblFinder.class") +
  ggtitle("Doublet Detection")
ggsave(file.path(fig_dir, "02_doublets.png"), p_dbl, width = 8, height = 6)

cat("  Before:", ncol(sc), "cells\n")
sc <- subset(sc, scDblFinder.class == "singlet")
cat("  After:", ncol(sc), "cells\n")

# ============================================
# 5. Filter
# ============================================
cat("Filtering...\n")
sc <- subset(sc, nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 10)
cat("  After filtering:", ncol(sc), "cells\n")

# ============================================
# 6. Normalize & Scale
# ============================================
cat("Normalizing...\n")
sc <- NormalizeData(sc)
sc <- FindVariableFeatures(sc, nfeatures = 2000)
sc <- ScaleData(sc)

# ============================================
# 7. Dimensionality Reduction
# ============================================
cat("Running PCA...\n")
sc <- RunPCA(sc, npcs = 50)

p_elbow <- ElbowPlot(sc, ndims = 50)
ggsave(file.path(fig_dir, "03_elbow_plot.png"), p_elbow, width = 6, height = 4)

# ============================================
# 8. Clustering & UMAP
# ============================================
cat("Clustering...\n")
sc <- FindNeighbors(sc, dims = 1:12)
sc <- FindClusters(sc, resolution = 0.5)
sc <- RunUMAP(sc, dims = 1:12)

p_umap_clusters <- DimPlot(sc, group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("Clusters")
ggsave(file.path(fig_dir, "04_umap_clusters.png"), p_umap_clusters, width = 10, height = 8)

# ============================================
# 9. Automated Annotation (SingleR)
# ============================================
cat("Running SingleR...\n")
sce <- as.SingleCellExperiment(sc)
ref <- celldex::HumanPrimaryCellAtlasData()
singler_results <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main,
  assay.type.test = "logcounts"
)
sc$singler_labels <- singler_results$pruned.labels

p_singler <- DimPlot(sc, group.by = "singler_labels", label = TRUE, repel = TRUE, size = 0.5) +
  ggtitle("SingleR Annotation") +
  NoLegend()
ggsave(file.path(fig_dir, "05_umap_singler.png"), p_singler, width = 12, height = 8)

# ============================================
# 10. Manual Cell State Annotation
# ============================================
cat("Annotating cell states...\n")
canonical_markers <- list(
  "Neoplastic" = c("EGFR", "VIM", "SOX2"),
  "Myeloid" = c("CSF1R", "CD68", "C1QB"),
  "T-Cells" = c("CD3D", "CD3E", "IL7R"),
  "Oligodendrocytes" = c("MBP", "PLP1", "MOG")
)

p_dot <- DotPlot(sc, features = unlist(canonical_markers), group.by = "seurat_clusters") +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "06_dotplot_markers.png"), p_dot, width = 10, height = 6)

# Manual labeling based on DotPlot + SingleR
sc$cell_state <- "Unknown"
sc$cell_state[sc$seurat_clusters %in% c(0, 1, 5, 9, 12, 21)] <- "Neoplastic"
sc$cell_state[sc$seurat_clusters %in% c(10, 11, 19)] <- "Myeloid"
sc$cell_state[sc$seurat_clusters %in% c(20)] <- "T-Cells"
sc$cell_state[sc$seurat_clusters %in% c(22)] <- "Oligodendrocytes"
sc$cell_state[sc$seurat_clusters %in% c(2, 24)] <- "Unknown"

p_umap_states <- DimPlot(sc, group.by = "cell_state", label = TRUE, repel = TRUE) +
  ggtitle("Cell States")
ggsave(file.path(fig_dir, "07_umap_cell_states.png"), p_umap_states, width = 10, height = 8)

# ============================================
# 11. Save
# ============================================
cat("Saving...\n")
saveRDS(sc, file.path(output_dir, "gbm_seurat_object.rds"))

cat("Done!\n")
