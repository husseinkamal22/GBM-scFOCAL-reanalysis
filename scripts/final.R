# ============================================
# GBM scFOCAL Reanalysis - Final Pipeline
# (نسخة مصححة - متسقة مع خطوات الكورس)
# ============================================

library(Seurat)
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scDblFinder)
library(SingleCellExperiment)
library(SingleR)
library(celldex)      # كانت ناقصة - لازمة لـ HumanPrimaryCellAtlasData()
library(BiocSingular)

set.seed(100)          # عشان الكلسترينج يطلع نفس الأرقام كل مرة تشغّل السكريبت

# --- Paths ---
base_dir   <- "D:/GBM-scFOCAL-reanalysis"
data_dir   <- file.path(base_dir, "data")
output_dir <- file.path(base_dir, "output")
fig_dir    <- file.path(base_dir, "figures")
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
rm(counts, counts_mat)

# لازم as.data.frame - fread بترجع data.table ورownames<- بتتجاهل عليها بصمت
meta <- as.data.frame(fread(file.path(data_dir, "GSE229779_cellMetadata.tsv")))
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

p_qc <- VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, pt.size = 0)   # pt.size=0 - من غير كده النقط بتغطي الـ violins بالكامل
ggsave(file.path(fig_dir, "01_qc_violin.png"), p_qc, width = 12, height = 4)

# ============================================
# 4. Filter (قبل الـ doublet detection)
# ============================================
cat("Filtering...\n")
sc <- subset(sc, nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 10)
cat("  After filtering:", ncol(sc), "cells\n")

# ============================================
# 5. Doublet Detection (scDblFinder)
# ============================================
cat("Detecting doublets...\n")

# منشغلش scDblFinder على sc مباشرة - بتحول الـ object لـ SingleCellExperiment
# وتبوظ الـ Seurat object بتاعنا. بنشتغل على نسخة SCE منفصلة بس.
sce_dbl <- as.SingleCellExperiment(sc)
sce_dbl <- scDblFinder(sce_dbl, dbr = 0.01)

sc$doublet_score <- colData(sce_dbl)$scDblFinder.score
sc$doublet_class <- colData(sce_dbl)$scDblFinder.class
rm(sce_dbl)

p_dbl <- DimPlot(sc, group.by = "doublet_class") + ggtitle("Doublet Detection")
ggsave(file.path(fig_dir, "02_doublets.png"), p_dbl, width = 8, height = 6)

cat("  Before:", ncol(sc), "cells\n")
sc <- subset(sc, doublet_class == "singlet")
cat("  After:", ncol(sc), "cells\n")

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
# 9. Manual Annotation (canonical markers + DotPlot + recode)
# ============================================
cat("Manual annotation...\n")

# ماركرز GBM الفعلية (مش PBMC) - Neoplastic/Myeloid/T cells/Oligodendrocytes/Endothelial/Fibroblast
canonical_markers <- c(
  "EGFR", "SOX2", "OLIG1",        # Neoplastic
  "PTPRC", "AIF1", "CD14",        # Myeloid
  "CD3D", "CD3E", "CD2",          # T cells
  "MBP", "PLP1", "MOG",           # Oligodendrocytes
  "PECAM1", "VWF", "CLDN5",       # Endothelial
  "COL1A1", "DCN"                 # Fibroblast
)

marker_to_type <- c(
  EGFR = "Neoplastic", SOX2 = "Neoplastic", OLIG1 = "Neoplastic",
  PTPRC = "Myeloid", AIF1 = "Myeloid", CD14 = "Myeloid",
  CD3D = "T cells", CD3E = "T cells", CD2 = "T cells",
  MBP = "Oligodendrocytes", PLP1 = "Oligodendrocytes", MOG = "Oligodendrocytes",
  PECAM1 = "Endothelial_Fibroblast", VWF = "Endothelial_Fibroblast", CLDN5 = "Endothelial_Fibroblast",
  COL1A1 = "Endothelial_Fibroblast", DCN = "Endothelial_Fibroblast"
)

p_dot <- DotPlot(sc, features = canonical_markers, group.by = "seurat_clusters") +
  RotatedAxis()
ggsave(file.path(fig_dir, "05_dotplot_markers.png"), p_dot, width = 10, height = 6)

# نفس منطق الـ recode اللي عملناه في الكورس، بس مبني على أعلى ماركر فعلي
# لكل كلاستر (مش أرقام كلاسترات مكتوبة بإيدينا) عشان يفضل صح حتى لو
# ترقيم الكلاسترات اتغير مع أي تعديل في الداتا أو الخطوات قبله.
dotplot_data <- p_dot$data

top_marker_per_cluster <- dotplot_data %>%
  group_by(id) %>%
  slice_max(avg.exp.scaled, n = 1, with_ties = FALSE) %>%
  arrange(id)

cluster_to_type <- setNames(
  marker_to_type[as.character(top_marker_per_cluster$features.plot)],
  as.character(top_marker_per_cluster$id)
)

sc$canonical_annotation <- dplyr::recode(
  as.character(sc$seurat_clusters),
  !!!cluster_to_type
)

# ⚠️ خطوة مراجعة يدوية: افتح 05_dotplot_markers.png وقارنه بـ
# canonical_annotation النهائي - لو أي كلاستر إشارته ضعيفة جدًا
# (avg.exp.scaled منخفض)، التصنيف الأوتوماتيكي هنا ممكن يكون غلط
# وهيحتاج مراجعة بإيدك زي ما عملنا مع كلاستر 0 و22 في الكورس.

p_umap_states <- DimPlot(sc, group.by = "canonical_annotation", label = TRUE, repel = TRUE) +
  ggtitle("Manual Annotation")
ggsave(file.path(fig_dir, "06_umap_manual_annotation.png"), p_umap_states, width = 10, height = 8)

# ============================================
# 10. DEG-based annotation (FindAllMarkers) - تأكيد إضافي
# ============================================
cat("Finding markers per cluster (this can take a few minutes)...\n")
Idents(sc) <- "seurat_clusters"

markers_all <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers_all, file.path(output_dir, "all_markers.csv"), row.names = FALSE)

top5 <- markers_all %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5)
write.csv(top5, file.path(output_dir, "top5_markers.csv"), row.names = FALSE)

# ============================================
# 11. Automated Annotation (SingleR)
# ============================================
cat("Running SingleR...\n")
sce <- as.SingleCellExperiment(sc)
hpca_ref <- celldex::HumanPrimaryCellAtlasData()

singler_results <- SingleR(
  test = sce,
  ref = hpca_ref,
  labels = hpca_ref$label.main
)
sc$singler_annotation <- singler_results$labels
rm(sce)

p_singler <- DimPlot(sc, group.by = "singler_annotation", label = TRUE, repel = TRUE) +
  ggtitle("SingleR Annotation")
ggsave(file.path(fig_dir, "07_umap_singler.png"), p_singler, width = 12, height = 8)

# مقارنة الاتنين جنب بعض
p_compare <- p_umap_states + p_singler
ggsave(file.path(fig_dir, "08_manual_vs_singler.png"), p_compare, width = 16, height = 7, dpi = 300)

# ============================================
# 12. تحقق من صحة الـ annotation مقابل الـ ground truth بتاع البيبر
# ============================================
cat("Validating against paper's own CellType column...\n")
if ("CellType" %in% colnames(sc@meta.data)) {
  validation_table <- table(sc$canonical_annotation, sc$CellType)
  write.csv(as.data.frame.matrix(validation_table),
            file.path(output_dir, "annotation_validation.csv"))
  print(validation_table)
}

# ============================================
# 13. Save
# ============================================
cat("Saving...\n")
saveRDS(sc, file.path(output_dir, "gbm_seurat_object.rds"))

cat("Done!\n")
