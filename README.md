# GFAP-dependent heterogeneity in glioma

Code accompanying:

> Hadizadeh M, Mousavi SZ, Ghasemi S.  
> GFAP-dependent transcriptional dynamics and cellular heterogeneity in primary, recurrent, and grade III gliomas.  
> *Cancer Informatics*. 2026;25:11769351261460625.  
> https://doi.org/10.1177/11769351261460625

Single-cell RNA-seq analysis of GFAP-related transcriptional dynamics and cellular heterogeneity across primary GBM, recurrent glioma, and grade III samples from GEO [GSE103224](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103224).

## Repository layout

- `scripts/scRNA_seurat.R` — QC, Harmony integration, clustering, annotation, and differential expression
- `scripts/trajectory.R` — Monocle3 trajectory analysis
- `data/` — input count matrices (not included; see below)
- `results/` — outputs written by the scripts (RDS, CSV, and plots)

## Data

Download the filtered count matrices from GEO series [GSE103224](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103224) (supplementary files on the series page or the individual GSM records). Place them under `data/` with this layout and these exact file names:

```text
data/
├── GBM/
│   ├── GSM2758471_PJ016.filtered.matrix.txt
│   ├── GSM2758472_PJ017.filtered.matrix.txt
│   ├── GSM2758473_PJ018.filtered.matrix.txt
│   ├── GSM2758474_PJ025.filtered.matrix.txt
│   └── GSM2940098_PJ048.filtered.matrix.txt
├── Rec/
│   ├── GSM2758476_PJ032.filtered.matrix.txt
│   └── GSM2758477_PJ035.filtered.matrix.txt
└── grade_III/
    └── GSM2758475_PJ030.filtered.matrix.txt
```

Run both scripts from the repository root so these relative paths resolve correctly.

## Requirements

R ≥ 4.2 and the packages below.

```r
install.packages(c(
  "Seurat", "tidyverse", "harmony", "patchwork", "cowplot",
  "data.table", "BiocManager", "ggridges", "scCustomize", "pheatmap"
))
BiocManager::install(c("metap", "SingleCellExperiment", "monocle3"))

# SeuratWrappers (needed by scripts/trajectory.R)
install.packages("remotes")
remotes::install_github("satijalab/seurat-wrappers")
```

## How to reproduce

From the repository root:

```r
source("scripts/scRNA_seurat.R")  # QC, integration, clustering, annotation, DE
source("scripts/trajectory.R")    # Monocle3 trajectory (loads results/seurat_Rec.rds by default)
```

`scRNA_seurat.R` writes annotated condition subsets to `results/seurat_GBM.rds`, `results/seurat_Rec.rds`, and `results/seurat_gIII.rds`. `trajectory.R` loads one of these (default: Rec). Edit the `condition` setting at the top of `trajectory.R` to run the other subsets. All outputs are written to `results/`.

## Citation

Please cite the paper above.
