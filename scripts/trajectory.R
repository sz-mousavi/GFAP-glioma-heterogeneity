# Monocle3 trajectory for GFAP-related glioma subsets
# Cleaned from the lab script used for Radial & cycling cells in recurrence.
# Hadizadeh M, Mousavi SZ, Ghasemi S.
# Cancer Informatics. 2026;25:11769351261460625.
# https://doi.org/10.1177/11769351261460625
# Expected input: a Seurat .rds with UMAP and identities already defined.

library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
library(SingleCellExperiment)

#data<- obj
input_rds <- "Radial&CC in Rec.RDS"
data <- readRDS(input_rds)

# ---------------------------------------------------------------------------
# Delete Cluster(s)
# ---------------------------------------------------------------------------
table(obj$seurat_clusters)
seurat_obj_filtered <- subset(obj, idents = c(7, 8), invert = TRUE)
table(seurat_obj_filtered$seurat_clusters)
data <- sub_obj

#If needed
data <- JoinLayers(data)

# ---------------------------------------------------------------------------
# Add Count Layer to sobj
# ---------------------------------------------------------------------------
count.data <- GetAssayData(object = data, assay = "RNA", layer = "count")
class(count.data)
data[["RNA_counts"]] <- CreateAssayObject(counts = count.data)
DefaultAssay(data) <- "RNA_counts"

# ---------------------------------------------------------------------------
# Converting Seurat object to celldataset object for Monocle3
# ---------------------------------------------------------------------------
cds <- as.cell_data_set(data)
cds <- cluster_cells(cds)

# ---------------------------------------------------------------------------
#cds@colData$ident <- factor(cds@colData$ident,
#                            levels = c("Epithelial_PTC", "Epithelial_ATC"))
#plot_cells(cds, color_cells_by = "celltype.status", group_label_size = 3)

#Get cell metadata
head(colData(cds))
#Get feature/gene metadata
fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))
head(fData(cds))
head(pData(cds))
#Get counts
head(counts(cds))
# ---------------------------------------------------------------------------
# Retrieve clustering information from Surat object
# ---------------------------------------------------------------------------
#Assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

#Assign cluster information
list.cluster <- data@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

#Assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- data@reductions$umap@cell.embeddings

#Plot
cluster.before.traj <- plot_cells(cds,
                                  color_cells_by = "cluster",
                                  label_groups_by_cluster = TRUE,
                                  group_label_size = 5) +
  theme(
    legend.position = "right",
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  )
cluster.before.traj

# ---------------------------------------------------------------------------
# Learn Trajectory
# ---------------------------------------------------------------------------
cds <- learn_graph(cds, use_partition = F)
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster = FALSE,
           label_branch_points = TRUE,
           label_roots = TRUE,
           label_leaves = FALSE,
           group_label_size = 5) +
  theme(
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  )

# ---------------------------------------------------------------------------
# Order cells in Pseudotime
# ---------------------------------------------------------------------------
#cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == 5]))
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = NULL)
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T, graph_label_size = 3,
           label_branch_points = F, label_roots = T, label_leaves = F)

#cds<- saveRDS(cds, "cds.RDS")
#cds<- readRDS("cds.RDS")

#Cells ordered by Monocle3 Pseudotime
head(pseudotime(cds), 10)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, ident, fill = ident)) +
  geom_boxplot() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 14)
  )

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(ident, monocle3_pseudotime), fill = ident)) +
  geom_boxplot() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 14)
  )
# ---------------------------------------------------------------------------
# Find genes that change as a function of pseudotime
# ---------------------------------------------------------------------------
deg <- graph_test(cds, neighbor_graph = "principal_graph")

saveRDS(deg, "deg(PT).RDS")
deg <- readRDS("deg-Rec.RDS")
deg %>% arrange(q_value) %>% filter(status == "OK") %>% head()
write.csv(deg, "deg(PT).csv")
head(deg)

# ---------------------------------------------------------------------------
# Module & HeaatMap
# ---------------------------------------------------------------------------
pr_deg_ids <- row.names(subset(deg, q_value < 0.05 & morans_I > 0.1))
head(pr_deg_ids)
cds <- preprocess_cds(cds)
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution = 1e-2)
head(gene_module_df)
write.csv(gene_module_df, "gene_module_df(PT).csv")
colData(cds)
cell_group_df <- tibble::tibble(cell = row.names(colData(cds)),
                                cell_group = colData(cds)$ident)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale = "column", clustering_method = "ward.D2")

# ---------------------------------------------------------------------------
# Feature & Ridge Plot
# ---------------------------------------------------------------------------
FeaturePlot(data, features = c("APOE", "ID4", "AQP4", "ITM2C", "CLU", "ID3"))
#RidgePlot(data, features = c("PRSS1", "PID1", "MEG8"), sort = T, idents = c("Epithelial", "DC", "Endothelial"))

# ---------------------------------------------------------------------------
# Expression in Pseudotime (ok)
# ---------------------------------------------------------------------------
my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("APOE", "ID4", "AQP4", "ITM2C", "CLU", "ID3")))
cds_subset <- cds[my_genes,]
cds_subset@colData$monocle3_pseudotime <- clusters(cds_subset)
plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime", cell_size = 2)

# ---------------------------------------------------------------------------
# Expression in Pseudotime (not-changed)
# ---------------------------------------------------------------------------
plot_genes_in_pseudotime(subset(cds_subset,
                                rowData(cds_subset)$gene_short_name %in% my_genes),
                         min_expr = 0.5, color_cells_by = "monocle3_pseudotime")
# ---------------------------------------------------------------------------
# bold dot Plot
# ---------------------------------------------------------------------------
plot_cells(cds_subset, label_cell_groups = FALSE, genes = my_genes,
           show_trajectory_graph = FALSE, cell_size = 1)

# ---------------------------------------------------------------------------
# Gene expression of a small number of genes using dot plot
# ---------------------------------------------------------------------------
library(scCustomize)
scCustomize::Clustered_DotPlot(data, features = my_genes,
                               raster = F, exp_color_min = -5, exp_color_max = 5,
                               column_label_size = 16, legend_label_size = 10, legend_title_size = 10,
                               row_label_size = 12, plot_km_elbow = FALSE)
# ---------------------------------------------------------------------------
