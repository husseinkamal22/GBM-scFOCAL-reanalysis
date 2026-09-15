# GBM scFOCAL Reanalysis

Single-cell RNA-seq reanalysis of Glioblastoma (GBM) data from the scFOCAL study.

## Datasets

- **GSE229779** — scRNA-seq of GBM patient samples (37,945 cells after QC)

## Pipeline

1. Quality control (mitochondrial %, gene/cell counts)
2. Doublet detection (scDblFinder, dbr = 0.01)
3. Filtering (>200 genes, <4500 genes, <10% MT)
4. Normalization, variable features (2000), scaling
5. PCA (50 PCs) + clustering (resolution 0.5, dims 1:12)
6. Automated annotation (SingleR, HumanPrimaryCellAtlas)
7. Manual cell state annotation (Neoplastic / Myeloid / T-Cells / Oligodendrocytes)

## Requirements

```r
install.packages(c("Seurat", "data.table", "ggplot2", "dplyr", "patchwork"))

# Bioconductor
BiocManager::install(c("scDblFinder", "SingleR", "celldex", "BiocSingular"))
```

## Usage

1. Download data from GEO (GSE229779) into `data/`
2. Run `scripts/final.R`

## Output

- `output/gbm_seurat_object.rds` — Processed Seurat object
- `figures/` — QC, doublet, UMAP, dotplot, and annotation plots
