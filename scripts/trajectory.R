################################################################################
# Monocle3 trajectory analysis
# Input: an annotated Seurat object written by scripts/scRNA_seurat.R
# Hadizadeh M, Mousavi SZ, Ghasemi S.
# Cancer Informatics. 2026;25:11769351261460625.
# https://doi.org/10.1177/11769351261460625
################################################################################

# ==============================================================================
# CONFIG — run this script from the repository root, after scRNA_seurat.R
# ==============================================================================
results_dir <- "results"

# Which condition subset to use. Files written by scRNA_seurat.R:
#   results/seurat_GBM.rds
#   results/seurat_Rec.rds
#   results/seurat_gIII.rds
condition <- "Rec"

# Optional: restrict to these Seurat identities. Set to NULL to use all cells.
# The original lab script ran this trajectory on radial glial and cancer
# cells in the recurrent-glioma subset. Edit to match the figure you want.
keep_idents <- c("Radial glial", "Cancer cells")

# Cell type used to place the trajectory root (non-interactive).
# This follows the standard Monocle3 vignette pattern (principal-graph node
# closest to the most cells of this identity). Set to NULL to choose the
# root interactively with the Monocle3 shiny app instead.
root_ident <- "Radial glial"

# Example genes plotted along pseudotime in the original analysis.
example_genes <- c("APOE", "ID4", "AQP4", "ITM2C", "CLU", "ID3")

# ==============================================================================
# Libraries
# ==============================================================================
library(Seurat)
library(tidyverse)
library(monocle3)
library(SeuratWrappers)
library(SingleCellExperiment)
library(ggplot2)
library(pheatmap)

# ==============================================================================
# Load the Seurat object written by scRNA_seurat.R
# ==============================================================================
input_rds <- file.path(results_dir, paste0("seurat_", condition, ".rds"))
if (!file.exists(input_rds)) {
  stop(
    "Cannot find ", input_rds, ".\n",
    "Run scripts/scRNA_seurat.R first, from the repository root."
  )
}

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
out_prefix <- file.path(results_dir, paste0("trajectory_", condition))

obj <- readRDS(input_rds)

# Seurat v5: combine layers so Monocle can see a single counts matrix
if (exists("JoinLayers", mode = "function")) {
  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
}

if (!is.null(keep_idents)) {
  present <- levels(Idents(obj))
  missing <- setdiff(keep_idents, present)
  if (length(missing)) {
    stop(
      "These identities are not in the object: ",
      paste(missing, collapse = ", "),
      "\nAvailable identities: ", paste(present, collapse = ", ")
    )
  }
  obj <- subset(obj, idents = keep_idents)
}

# ==============================================================================
# Convert to a Monocle3 cell_data_set; keep Seurat UMAP and cluster labels
# ==============================================================================
cds <- as.cell_data_set(obj)
reducedDim(cds, "UMAP") <- Embeddings(obj, "umap")

# Single partition so learn_graph builds one trajectory
recreate_partitions <- factor(rep(1, ncol(cds)))
names(recreate_partitions) <- colnames(cds)
cds@clusters[["UMAP"]][["partitions"]] <- recreate_partitions

seurat_idents <- Idents(obj)
names(seurat_idents) <- colnames(obj)
cds@clusters[["UMAP"]][["clusters"]] <- seurat_idents

fData(cds)$gene_short_name <- rownames(fData(cds))

p <- plot_cells(
  cds,
  color_cells_by = "cluster",
  label_groups_by_cluster = TRUE,
  group_label_size = 5
) +
  theme(
    legend.position = "right",
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  )
ggsave(paste0(out_prefix, "_clusters.png"), p, width = 8, height = 6)

# ==============================================================================
# Learn the trajectory and order cells in pseudotime
# ==============================================================================
cds <- learn_graph(cds, use_partition = FALSE)

p <- plot_cells(
  cds,
  color_cells_by = "cluster",
  label_groups_by_cluster = FALSE,
  label_branch_points = TRUE,
  label_roots = TRUE,
  label_leaves = FALSE,
  group_label_size = 5
) +
  theme(
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  )
ggsave(paste0(out_prefix, "_graph.png"), p, width = 8, height = 6)

# Principal-graph node closest to the most cells of root_ident
# (Monocle3 vignette pattern), so order_cells() can run non-interactively.
pick_root_pr_nodes <- function(cds, root_ident) {
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), , drop = FALSE])

  cluster_labels <- as.character(clusters(cds))
  names(cluster_labels) <- names(clusters(cds))
  in_root <- cluster_labels[colnames(cds)] == root_ident
  if (!any(in_root) && "ident" %in% colnames(colData(cds))) {
    in_root <- as.character(colData(cds)$ident) == root_ident
  }
  if (!any(in_root)) {
    stop("No cells with identity '", root_ident, "' were found to place the root.")
  }

  root_counts <- table(closest_vertex[in_root, 1])
  igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(root_counts)))
  ]
}

if (is.null(root_ident)) {
  cds <- order_cells(cds, reduction_method = "UMAP")
} else {
  root_pr_nodes <- pick_root_pr_nodes(cds, root_ident)
  cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes = root_pr_nodes)
}

p <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_branch_points = FALSE,
  label_roots = TRUE,
  label_leaves = FALSE,
  graph_label_size = 3
)
ggsave(paste0(out_prefix, "_pseudotime.png"), p, width = 8, height = 6)

colData(cds)$monocle3_pseudotime <- pseudotime(cds)
data_pseudo <- as.data.frame(colData(cds))
data_pseudo$ident <- as.character(clusters(cds)[rownames(data_pseudo)])

p <- ggplot(
  data_pseudo,
  aes(monocle3_pseudotime, reorder(ident, monocle3_pseudotime), fill = ident)
) +
  geom_boxplot() +
  labs(x = "Pseudotime", y = "Identity") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )
ggsave(paste0(out_prefix, "_pseudotime_boxplot.png"), p, width = 8, height = 6)

saveRDS(cds, paste0(out_prefix, "_cds.rds"))

# ==============================================================================
# Genes that change as a function of pseudotime
# ==============================================================================
deg <- graph_test(cds, neighbor_graph = "principal_graph")
deg <- deg[order(deg$q_value), ]
write.csv(deg, paste0(out_prefix, "_pseudotime_degs.csv"))

# ==============================================================================
# Gene modules among significant trajectory genes
# ==============================================================================
pr_deg_ids <- rownames(subset(deg, q_value < 0.05 & morans_I > 0.1))
if (length(pr_deg_ids) >= 2) {
  cds <- preprocess_cds(cds)
  gene_module_df <- find_gene_modules(cds[pr_deg_ids, ], resolution = 1e-2)
  write.csv(gene_module_df, paste0(out_prefix, "_gene_modules.csv"), row.names = FALSE)

  cell_group_df <- tibble::tibble(
    cell = rownames(colData(cds)),
    cell_group = as.character(clusters(cds)[rownames(colData(cds))])
  )
  agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
  rownames(agg_mat) <- paste("Module", rownames(agg_mat))

  png(
    paste0(out_prefix, "_module_heatmap.png"),
    width = 8, height = 6, units = "in", res = 150
  )
  pheatmap::pheatmap(agg_mat, scale = "column", clustering_method = "ward.D2")
  dev.off()
} else {
  message("Not enough significant trajectory genes for module analysis.")
}

# ==============================================================================
# Example gene plots
# ==============================================================================
available <- example_genes[example_genes %in% fData(cds)$gene_short_name]
if (length(available) == 0) {
  message("None of the example genes were found in the object; skipping gene plots.")
} else {
  p <- FeaturePlot(obj, features = available)
  ggsave(paste0(out_prefix, "_example_featureplot.png"), p, width = 10, height = 8)

  cds_subset <- cds[rownames(subset(fData(cds), gene_short_name %in% available)), ]
  p <- plot_genes_in_pseudotime(
    cds_subset,
    color_cells_by = "pseudotime",
    cell_size = 1
  )
  ggsave(paste0(out_prefix, "_example_genes_pseudotime.png"), p, width = 8, height = 8)

  p <- plot_cells(
    cds_subset,
    label_cell_groups = FALSE,
    genes = available,
    show_trajectory_graph = FALSE,
    cell_size = 1
  )
  ggsave(paste0(out_prefix, "_example_genes_umap.png"), p, width = 10, height = 8)

  if (requireNamespace("scCustomize", quietly = TRUE)) {
    png(
      paste0(out_prefix, "_example_dotplot.png"),
      width = 8, height = 6, units = "in", res = 150
    )
    scCustomize::Clustered_DotPlot(
      obj,
      features = available,
      raster = FALSE,
      exp_color_min = -5,
      exp_color_max = 5,
      column_label_size = 16,
      legend_label_size = 10,
      legend_title_size = 10,
      row_label_size = 12,
      plot_km_elbow = FALSE
    )
    dev.off()
  } else {
    message("Package scCustomize is not installed; skipping the clustered dot plot.")
  }
}

message("Trajectory outputs written under ", results_dir, " with prefix trajectory_", condition)
