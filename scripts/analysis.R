# ============================================
# GBM scFOCAL Reanalysis - scRNA-seq Pipeline
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

# --- 1. Load Data ---
cat("Loading counts matrix...\n")
counts <- fread(file.path(data_dir, "GSE229779_countsMatrix.tsv"))
counts_mat <- as.matrix(counts[, -1])
rownames(counts_mat) <- counts$V1
counts_sparse <- as(counts_mat, "dgCMatrix")
rm(counts_mat)

cat("Loading metadata...\n")
meta <- fread(file.path(data_dir, "GSE229779_cellMetadata.tsv"))
rownames(meta) <- meta$cellID

# --- 2. Create Seurat Object ---
cat("Creating Seurat object...\n")
my_sc_data <- CreateSeuratObject(
  counts = counts_sparse,
  meta.data = meta,
  project = "GBM_reanalysis",
  min.cells = 3,
  min.features = 200
)
rm(counts_sparse)

# --- 3. QC ---
cat("Calculating QC metrics...\n")
my_sc_data[["percent.mt"]] <- PercentageFeatureSet(my_sc_data, pattern = "^MT-")

p_qc <- VlnPlot(my_sc_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(fig_dir, "qc_violin.png"), p_qc, width = 12, height = 4)

# --- 4. Filter ---
cat("Filtering cells: ", ncol(my_sc_data), " -> ", append = FALSE)
my_sc_data <- subset(
  my_sc_data,
  nFeature_RNA > 500 &
  nFeature_RNA < 4500 &
  percent.mt < 10
)
cat(ncol(my_sc_data), " cells\n")

# --- 5. Normalize & Scale ---
cat("Normalizing...\n")
my_sc_data <- NormalizeData(my_sc_data)
my_sc_data <- FindVariableFeatures(my_sc_data, nfeatures = 2000)
my_sc_data <- ScaleData(my_sc_data)

# --- 6. PCA ---
cat("Running PCA...\n")
my_sc_data <- RunPCA(my_sc_data, npcs = 50)

p_elbow <- ElbowPlot(my_sc_data, ndims = 50)
ggsave(file.path(fig_dir, "elbow_plot.png"), p_elbow, width = 6, height = 4)

# --- 7. Clustering & UMAP ---
cat("Clustering...\n")
my_sc_data <- FindNeighbors(my_sc_data, dims = 1:12)
my_sc_data <- FindClusters(my_sc_data, resolution = 0.5)
my_sc_data <- RunUMAP(my_sc_data, dims = 1:12)

# --- 8. Cell State Annotation ---
canonical_markers <- list(
  "AC-like" = c("GFAP", "AQP4", "S100B"),
  "NPC-like" = c("SOX2", "NES", "VIM"),
  "OPC-like" = c("OLIG2", "PDGFRA", "MBP"),
  "MES-like" = c("CHI3L1", "CD44", "ANGPTL4"),
  "Oligodendrocytes" = c("MBP", "PLP1", "MOG"),
  "Myeloid" = c("CSF1R", "CD68", "C1QB"),
  "T-Cells" = c("CD3D", "CD3E", "IL7R")
)

p_dot <- DotPlot(my_sc_data,
  features = unlist(canonical_markers),
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "dotplot_markers.png"), p_dot, width = 14, height = 6)

my_sc_data$cell_state <- NA
my_sc_data$cell_state[my_sc_data$seurat_clusters %in% c(12)] <- "NPC-like"
my_sc_data$cell_state[my_sc_data$seurat_clusters %in% c(4, 8, 13, 16)] <- "AC-like"
my_sc_data$cell_state[my_sc_data$seurat_clusters %in% c(10, 11, 19)] <- "Myeloid"
my_sc_data$cell_state[my_sc_data$seurat_clusters %in% c(3, 6, 7, 14, 17, 18, 23)] <- "OPC-like"
my_sc_data$cell_state[my_sc_data$seurat_clusters %in% c(0, 1, 5, 9, 21)] <- "MES-like"
my_sc_data$cell_state[my_sc_data$seurat_clusters %in% c(22)] <- "Oligodendrocytes"
my_sc_data$cell_state[my_sc_data$seurat_clusters %in% c(20)] <- "T-Cells"
my_sc_data$cell_state[my_sc_data$seurat_clusters %in% c(15)] <- "Immune"
my_sc_data$cell_state[my_sc_data$seurat_clusters %in% c(2, 24)] <- "Unknown"

p_umap <- DimPlot(my_sc_data, group.by = "cell_state", label = TRUE, repel = TRUE) +
  ggtitle("Cell State Annotation")
ggsave(file.path(fig_dir, "umap_cell_states.png"), p_umap, width = 10, height = 8)

# --- 9. Find Markers ---
cat("Finding markers...\n")
markers <- FindAllMarkers(
  my_sc_data,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5)

# --- 10. Save ---
cat("Saving results...\n")
write.csv(top5, file.path(output_dir, "top5_markers.csv"), row.names = FALSE)
write.csv(markers, file.path(output_dir, "all_markers.csv"), row.names = FALSE)
saveRDS(my_sc_data, file.path(output_dir, "gbm_seurat_object.rds"))

cat("Done!\n")
