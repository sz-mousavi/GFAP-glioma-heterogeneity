################################################################################
# GSE103224 — Seurat + Harmony + DEG
# Single-cell RNA-seq analysis of high-grade glioma samples from GEO GSE103224
# Dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103224
# ==============================================================================
# CONFIG — edit these two lines for your machine
# ==============================================================================
data_dir    <- "data"      # folder that contains GBM/, Rec/, grade_III/
results_dir <- "results"   # output folder for plots, CSV, RDS

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Libraries
# ==============================================================================
library(Seurat)
library(tidyverse)
library(BiocManager)
library(metap)
library(harmony)
library(patchwork)
library(cowplot)
library(data.table)
library(dplyr)

# ==============================================================================
# Helper: read one filtered matrix (same logic as the original script)
# ==============================================================================
read_filtered_matrix <- function(path) {
  mat <- fread(path, data.table = FALSE)
  print(head(mat[, 1:4]))
  mat <- mat[, -1]
  mat <- mat %>% distinct(V2, .keep_all = TRUE)
  rownames(mat) <- mat[, 1]
  mat <- mat[, -1]
  mat
}

# ==============================================================================
# Load count matrices
# ==============================================================================
GBM1 <- read_filtered_matrix(file.path(data_dir, "GBM/GSM2758471_PJ016.filtered.matrix.txt"))
GBM2 <- read_filtered_matrix(file.path(data_dir, "GBM/GSM2758472_PJ017.filtered.matrix.txt"))
GBM3 <- read_filtered_matrix(file.path(data_dir, "GBM/GSM2758473_PJ018.filtered.matrix.txt"))
GBM4 <- read_filtered_matrix(file.path(data_dir, "GBM/GSM2758474_PJ025.filtered.matrix.txt"))
GBM5 <- read_filtered_matrix(file.path(data_dir, "GBM/GSM2940098_PJ048.filtered.matrix.txt"))

Rec1 <- read_filtered_matrix(file.path(data_dir, "Rec/GSM2758476_PJ032.filtered.matrix.txt"))
Rec2 <- read_filtered_matrix(file.path(data_dir, "Rec/GSM2758477_PJ035.filtered.matrix.txt"))

gIII <- read_filtered_matrix(file.path(data_dir, "grade_III/GSM2758475_PJ030.filtered.matrix.txt"))

# ==============================================================================
# Create Seurat objects + mitochondrial content
# ==============================================================================
GBM1 <- CreateSeuratObject(counts = GBM1, project = "GBM", min.cells = 3, min.features = 200)
GBM1 <- PercentageFeatureSet(GBM1, pattern = "^MT-", col.name = "percent.mt")

GBM2 <- CreateSeuratObject(counts = GBM2, project = "GBM", min.cells = 3, min.features = 200)
GBM2 <- PercentageFeatureSet(GBM2, pattern = "^MT-", col.name = "percent.mt")

GBM3 <- CreateSeuratObject(counts = GBM3, project = "GBM", min.cells = 3, min.features = 200)
GBM3 <- PercentageFeatureSet(GBM3, pattern = "^MT-", col.name = "percent.mt")

GBM4 <- CreateSeuratObject(counts = GBM4, project = "GBM", min.cells = 3, min.features = 200)
GBM4 <- PercentageFeatureSet(GBM4, pattern = "^MT-", col.name = "percent.mt")

GBM5 <- CreateSeuratObject(counts = GBM5, project = "GBM", min.cells = 3, min.features = 200)
GBM5 <- PercentageFeatureSet(GBM5, pattern = "^MT-", col.name = "percent.mt")

Rec1 <- CreateSeuratObject(counts = Rec1, project = "Rec", min.cells = 3, min.features = 200)
Rec1 <- PercentageFeatureSet(Rec1, pattern = "^MT-", col.name = "percent.mt")

Rec2 <- CreateSeuratObject(counts = Rec2, project = "Rec", min.cells = 3, min.features = 200)
Rec2 <- PercentageFeatureSet(Rec2, pattern = "^MT-", col.name = "percent.mt")

gIII <- CreateSeuratObject(counts = gIII, project = "grade_III", min.cells = 3, min.features = 200)
gIII <- PercentageFeatureSet(gIII, pattern = "^MT-", col.name = "percent.mt")

# ==============================================================================
# Filter low-quality cells (thresholds kept as in the original analysis)
# ==============================================================================
VlnPlot(GBM1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
GBM1 <- subset(GBM1, subset = nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 3)

VlnPlot(GBM2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
GBM2 <- subset(GBM2, subset = nFeature_RNA < 4000 & nCount_RNA < 10000 & percent.mt < 4)

VlnPlot(GBM3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
GBM3 <- subset(GBM3, subset = nFeature_RNA < 4000 & nCount_RNA < 10000 & percent.mt < 2)

VlnPlot(GBM4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
GBM4 <- subset(GBM4, subset = nFeature_RNA < 4000 & nCount_RNA < 10000 & percent.mt < 5)

VlnPlot(GBM5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
GBM5 <- subset(GBM5, subset = nFeature_RNA < 3000 & nCount_RNA < 7500 & percent.mt < 1)

VlnPlot(Rec1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Rec1 <- subset(Rec1, subset = nFeature_RNA < 4000 & nCount_RNA < 10000 & percent.mt < 2)

VlnPlot(Rec2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Rec2 <- subset(Rec2, subset = nFeature_RNA < 2000 & nCount_RNA < 4000 & percent.mt < 2)

VlnPlot(gIII, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
gIII <- subset(gIII, subset = nFeature_RNA < 3000 & nCount_RNA < 5000 & percent.mt < 2)

# ==============================================================================
# Merge samples and add status labels
# ==============================================================================
GBM <- merge(
  GBM1,
  c(GBM2, GBM3, GBM4, GBM5),
  add.cell.ids = c("GBM1", "GBM2", "GBM3", "GBM4", "GBM5")
)
Rec <- merge(Rec1, Rec2, add.cell.ids = c("Rec1", "Rec2"))

GBM$status  <- "GBM"
Rec$status  <- "Rec"
gIII$status <- "gIII"

obj <- merge(x = GBM, y = list(Rec, gIII))
rm(GBM1, GBM2, GBM3, GBM4, GBM5, Rec1, Rec2, gIII)

# ==============================================================================
# Normalize, HVG, scale, PCA, UMAP (before integration)
# ==============================================================================
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, verbose = FALSE)

DimPlot(obj, raster = FALSE, group.by = "status", reduction = "umap") +
  plot_annotation(title = "GBM, Rec and Grade_III before integration")

# ==============================================================================
# Harmony integration
# ==============================================================================
obj <- obj %>% RunHarmony("orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(obj, "harmony")
harmony_embeddings[1:5, 1:5]

p1 <- DimPlot(
  object = obj, raster = FALSE, reduction = "harmony",
  pt.size = 0.1, group.by = "orig.ident"
) + NoLegend()
p2 <- VlnPlot(
  object = obj, raster = FALSE, features = "harmony_1",
  group.by = "orig.ident", pt.size = 0.1
) + NoLegend()
plot_grid(p1, p2)

# saveRDS(obj, file.path(results_dir, "Harmony_newObj(27sep).RDS"))
# obj <- readRDS(file.path(results_dir, "Harmony_newObj.RDS"))

# ==============================================================================
# UMAP + clustering on Harmony embeddings
# ==============================================================================
obj <- obj %>%
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:30) %>%
  FindClusters(resolution = 0.1) %>%
  identity()

obj <- SetIdent(obj, value = "orig.ident")
DimPlot(obj, raster = FALSE, group.by = "status", reduction = "umap") +
  plot_annotation(title = "GBM, Rec and Grade_III after integration")
DimPlot(
  obj, raster = FALSE, reduction = "umap",
  group.by = "status", pt.size = 0.1, split.by = "status"
) + NoLegend()

obj <- SetIdent(obj, value = "seurat_clusters")
DimPlot(obj, raster = FALSE, label = TRUE) + NoLegend()
DimPlot(obj, raster = FALSE, reduction = "umap", split.by = "status", label = TRUE)

length(levels(obj))
length(levels(obj$seurat_clusters))

count_table <- table(obj@meta.data$seurat_clusters, obj@meta.data$status)
count_table
write.csv(count_table, file.path(results_dir, "count_table_3C_(res_0.1).csv"))

head(obj)
tail(obj)
saveRDS(obj, file.path(results_dir, "obj(p).RDS"))
obj <- readRDS(file.path(results_dir, "obj(p).RDS"))

# ==============================================================================
# Annotate clusters
# Clusters 4+5 → Cancer cells; 6+7 → Endothelial (same mapping as original)
# ==============================================================================
obj <- RenameIdents(
  obj,
  `0` = "Radial glial",
  `1` = "Immature neurons",
  `2` = "Mature neurons",
  `3` = "Microglial",
  `4` = "Cancer cells",
  `5` = "Cancer cells",
  `6` = "Endothelial",
  `7` = "Endothelial",
  `8` = "Oligodendrocytes",
  `9` = "Immune system"
)

obj_GBM  <- subset(x = obj, subset = status == "GBM")
obj_Rec  <- subset(x = obj, subset = status == "Rec")
obj_gIII <- subset(x = obj, subset = status == "gIII")

saveRDS(obj_GBM,  file.path(results_dir, "Traject&CellChat_GBM.RDS"))
saveRDS(obj_Rec,  file.path(results_dir, "Traject&CellChat_Rec.RDS"))
saveRDS(obj_gIII, file.path(results_dir, "Traject&CellChat_gIII.RDS"))

length(levels(obj))
length(levels(obj$seurat_clusters))
DimPlot(obj, raster = FALSE, label = TRUE)
DimPlot(obj, raster = FALSE, reduction = "umap", split.by = "status", label = TRUE)

# ==============================================================================
# Differential expression
# ==============================================================================
head(obj)
tail(obj)
obj <- JoinLayers(obj)

# Identify DEGs across conditions, per cell type
obj$celltype.status <- paste(Idents(obj), obj$status, sep = "_")
obj$celltype <- Idents(obj)
Idents(obj) <- "celltype.status"

length(levels(obj$celltype))
head(obj)
tail(obj)
saveRDS(obj, file.path(results_dir, "obj_forDEGs.RDS"))
obj <- readRDS(file.path(results_dir, "obj_forDEGs.RDS"))

# --- Rec vs GBM ---------------------------------------------------------------
GBMRec0 <- FindMarkers(obj, ident.1 = "Radial glial_Rec",        ident.2 = "Radial glial_GBM",        verbose = FALSE)
write.csv(GBMRec0, file.path(results_dir, "DEG-Radial glial(Rec-GBM).csv"))

GBMRec1 <- FindMarkers(obj, ident.1 = "Immature neurons_Rec",    ident.2 = "Immature neurons_GBM",    verbose = FALSE)
write.csv(GBMRec1, file.path(results_dir, "DEG-Immature neurons(Rec-GBM).csv"))

GBMRec2 <- FindMarkers(obj, ident.1 = "Mature neurons_Rec",      ident.2 = "Mature neurons_GBM",      verbose = FALSE)
write.csv(GBMRec2, file.path(results_dir, "DEG-Mature neurons(Rec-GBM).csv"))

GBMRec3 <- FindMarkers(obj, ident.1 = "Microglial_Rec",          ident.2 = "Microglial_GBM",          verbose = FALSE)
write.csv(GBMRec3, file.path(results_dir, "DEG-Microglial(Rec-GBM).csv"))

GBMRec4 <- FindMarkers(obj, ident.1 = "Cancer cells_Rec",        ident.2 = "Cancer cells_GBM",        verbose = FALSE)
write.csv(GBMRec4, file.path(results_dir, "DEG-Cancer cells(Rec-GBM).csv"))

GBMRec5 <- FindMarkers(obj, ident.1 = "Endothelial_Rec",         ident.2 = "Endothelial_GBM",         verbose = FALSE)
write.csv(GBMRec5, file.path(results_dir, "DEG-Endothelial(Rec-GBM).csv"))

GBMRec6 <- FindMarkers(obj, ident.1 = "Oligodendrocytes_Rec",    ident.2 = "Oligodendrocytes_GBM",    verbose = FALSE)
write.csv(GBMRec6, file.path(results_dir, "DEG-Oligodendrocytes(Rec-GBM).csv"))

GBMRec7 <- FindMarkers(obj, ident.1 = "Immune system_Rec",       ident.2 = "Immune system_GBM",       verbose = FALSE)
write.csv(GBMRec7, file.path(results_dir, "DEG-Immune system(Rec-GBM).csv"))

# --- GBM vs grade III ---------------------------------------------------------
gIIIGBM0 <- FindMarkers(obj, ident.1 = "Radial glial_GBM",       ident.2 = "Radial glial_gIII",       verbose = FALSE)
write.csv(gIIIGBM0, file.path(results_dir, "DEG-Radial glial(GBM-gIII).csv"))

gIIIGBM1 <- FindMarkers(obj, ident.1 = "Immature neurons_GBM",   ident.2 = "Immature neurons_gIII",   verbose = FALSE)
write.csv(gIIIGBM1, file.path(results_dir, "DEG-Immature neurons(GBM-gIII).csv"))

gIIIGBM2 <- FindMarkers(obj, ident.1 = "Mature neurons_GBM",     ident.2 = "Mature neurons_gIII",     verbose = FALSE)
write.csv(gIIIGBM2, file.path(results_dir, "DEG-Mature neurons(GBM-gIII).csv"))

gIIIGBM3 <- FindMarkers(obj, ident.1 = "Microglial_GBM",         ident.2 = "Microglial_gIII",         verbose = FALSE)
write.csv(gIIIGBM3, file.path(results_dir, "DEG-Microglial(GBM-gIII).csv"))

gIIIGBM4 <- FindMarkers(obj, ident.1 = "Cancer cells_GBM",       ident.2 = "Cancer cells_gIII",       verbose = FALSE)
write.csv(gIIIGBM4, file.path(results_dir, "DEG-Cancer cells(GBM-gIII).csv"))

gIIIGBM5 <- FindMarkers(obj, ident.1 = "Endothelial_GBM",        ident.2 = "Endothelial_gIII",        verbose = FALSE)
write.csv(gIIIGBM5, file.path(results_dir, "DEG-Endothelial(GBM-gIII).csv"))

gIIIGBM6 <- FindMarkers(obj, ident.1 = "Oligodendrocytes_GBM",   ident.2 = "Oligodendrocytes_gIII",   verbose = FALSE)
write.csv(gIIIGBM6, file.path(results_dir, "DEG-Oligodendrocytes(GBM-gIII).csv"))

gIIIGBM7 <- FindMarkers(obj, ident.1 = "Immune system_GBM",      ident.2 = "Immune system_gIII",      verbose = FALSE)
write.csv(gIIIGBM7, file.path(results_dir, "DEG-Immune system(GBM-gIII).csv"))

# --- Rec vs grade III ---------------------------------------------------------
gIIIRec0 <- FindMarkers(obj, ident.1 = "Radial glial_Rec",       ident.2 = "Radial glial_gIII",       verbose = FALSE)
write.csv(gIIIRec0, file.path(results_dir, "DEG-Radial glial(Rec-gIII).csv"))

gIIIRec1 <- FindMarkers(obj, ident.1 = "Immature neurons_Rec",   ident.2 = "Immature neurons_gIII",   verbose = FALSE)
write.csv(gIIIRec1, file.path(results_dir, "DEG-Immature neurons(Rec-gIII).csv"))

gIIIRec2 <- FindMarkers(obj, ident.1 = "Mature neurons_Rec",     ident.2 = "Mature neurons_gIII",     verbose = FALSE)
write.csv(gIIIRec2, file.path(results_dir, "DEG-Mature neurons(Rec-gIII).csv"))

gIIIRec3 <- FindMarkers(obj, ident.1 = "Microglial_Rec",         ident.2 = "Microglial_gIII",         verbose = FALSE)
write.csv(gIIIRec3, file.path(results_dir, "DEG-Microglial(Rec-gIII).csv"))

gIIIRec4 <- FindMarkers(obj, ident.1 = "Cancer cells_Rec",       ident.2 = "Cancer cells_gIII",       verbose = FALSE)
write.csv(gIIIRec4, file.path(results_dir, "DEG-Cancer cells(Rec-gIII).csv"))

gIIIRec5 <- FindMarkers(obj, ident.1 = "Endothelial_Rec",        ident.2 = "Endothelial_gIII",        verbose = FALSE)
write.csv(gIIIRec5, file.path(results_dir, "DEG-Endothelial(Rec-gIII).csv"))

gIIIRec6 <- FindMarkers(obj, ident.1 = "Oligodendrocytes_Rec",   ident.2 = "Oligodendrocytes_gIII",   verbose = FALSE)
write.csv(gIIIRec6, file.path(results_dir, "DEG-Oligodendrocytes(Rec-gIII).csv"))

gIIIRec7 <- FindMarkers(obj, ident.1 = "Immune system_Rec",      ident.2 = "Immune system_gIII",      verbose = FALSE)
write.csv(gIIIRec7, file.path(results_dir, "DEG-Immune system(Rec-gIII).csv"))

# --- Microglial vs Cancer cells (within condition) ----------------------------
Mic_Rec  <- FindMarkers(obj, ident.1 = "Microglial_Rec",  ident.2 = "Cancer cells_Rec",  verbose = FALSE)
write.csv(Mic_Rec,  file.path(results_dir, "DEG-Microglial-CC(Rec).csv"))

Mic_GBM  <- FindMarkers(obj, ident.1 = "Microglial_GBM",  ident.2 = "Cancer cells_GBM",  verbose = FALSE)
write.csv(Mic_GBM,  file.path(results_dir, "DEG-Microglial-CC(GBM).csv"))

Mic_gIII <- FindMarkers(obj, ident.1 = "Microglial_gIII", ident.2 = "Cancer cells_gIII", verbose = FALSE)
write.csv(Mic_gIII, file.path(results_dir, "DEG-Microglial-CC(gIII).csv"))

# --- Radial glial vs Cancer cells (within condition) --------------------------
Rad_Rec  <- FindMarkers(obj, ident.1 = "Radial glial_Rec",  ident.2 = "Cancer cells_Rec",  verbose = FALSE)
write.csv(Rad_Rec,  file.path(results_dir, "DEG-Radial glial-CC(Rec).csv"))

Rad_GBM  <- FindMarkers(obj, ident.1 = "Radial glial_GBM",  ident.2 = "Cancer cells_GBM",  verbose = FALSE)
write.csv(Rad_GBM,  file.path(results_dir, "DEG-Radial glial-CC(GBM).csv"))

Rad_gIII <- FindMarkers(obj, ident.1 = "Radial glial_gIII", ident.2 = "Cancer cells_gIII", verbose = FALSE)
write.csv(Rad_gIII, file.path(results_dir, "DEG-Radial glial-CC(gIII).csv"))
