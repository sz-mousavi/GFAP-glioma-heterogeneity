# GFAP-dependent heterogeneity in glioma

Code accompanying:

> Hadizadeh M, Mousavi SZ, Ghasemi S.  
> GFAP-dependent transcriptional dynamics and cellular heterogeneity in primary, recurrent, and grade III gliomas.  
> *Cancer Informatics*. 2026;25:11769351261460625.  
> https://doi.org/10.1177/11769351261460625

Single-cell RNA-seq analysis of GFAP-related transcriptional dynamics and cellular heterogeneity across primary GBM, recurrent glioma, and grade III samples from GEO [GSE103224](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103224).

## Repository layout

- `scripts/scRNA_seurat.R` — QC, Harmony integration, clustering, annotation, and differential expression
- `scripts/trajectory.R` — trajectory analysis

## Requirements

R ≥ 4.2 and:

```r
install.packages(c(
  "Seurat", "tidyverse", "harmony", "patchwork",
  "cowplot", "data.table", "BiocManager"
))
BiocManager::install("metap")

## Citation

Please cite the paper above.
