# GBM scFOCAL Reanalysis

Single-cell RNA-seq reanalysis of Glioblastoma (GBM) data from the scFOCAL study.

## Datasets

- **GSE229779** — scRNA-seq of GBM patient samples

## Pipeline

1. QC (mitochondrial %, gene/cell counts)
2. Filtering (>200 genes, <4500 genes, <10% MT)
3. Doublet detection (scDblFinder, dbr = 0.01)
4. Normalization, variable features (2000), scaling
5. PCA (50 PCs) + clustering (resolution 0.5, dims 1:12)
6. Manual annotation (canonical GBM markers + DotPlot + recode)
7. FindAllMarkers + top5 per cluster
8. Automated annotation (SingleR, HumanPrimaryCellAtlas)
9. Validation against paper's CellType column

## Requirements

```r
install.packages(c("Seurat", "data.table", "ggplot2", "dplyr", "patchwork"))

# Bioconductor
BiocManager::install(c("scDblFinder", "SingleR", "SingleCellExperiment", "celldex", "BiocSingular"))
```

## Usage

1. Download data from GEO (GSE229779) into `data/`
2. Run `scripts/final.R`

## Output

- `output/gbm_seurat_object.rds` — Processed Seurat object
- `output/all_markers.csv` — All marker genes
- `output/top5_markers.csv` — Top 5 markers per cluster
- `output/annotation_validation.csv` — Manual vs paper CellType comparison
- `figures/` — QC, doublet, UMAP, dotplot, and annotation plots
