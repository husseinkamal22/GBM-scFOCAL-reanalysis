# GBM scFOCAL Reanalysis

Single-cell RNA-seq reanalysis of Glioblastoma (GBM) data from the scFOCAL study.

## Datasets

- **GSE229779** — scRNA-seq of GBM patient samples
- **GSE231489** — Additional GBM scRNA-seq data

## Pipeline

1. Quality control (mitochondrial %, gene/cell counts)
2. Filtering (>500 genes, <4500 genes, <10% MT)
3. Normalization, variable features, scaling
4. PCA + clustering (resolution 0.5)
5. Cell state annotation (AC-like, NPC-like, OPC-like, MES-like, etc.)
6. Marker identification

## Requirements

```r
install.packages(c("Seurat", "data.table", "ggplot2", "dplyr", "patchwork"))
```

## Usage

1. Download data from GEO (GSE229779) into `data/`
2. Run `scripts/analysis.R`

## Output

- `output/gbm_seurat_object.rds` — Processed Seurat object
- `output/all_markers.csv` — All marker genes
- `output/top5_markers.csv` — Top 5 markers per cluster
- `figures/` — QC, UMAP, and marker plots
