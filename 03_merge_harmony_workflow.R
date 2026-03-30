#!/usr/bin/env Rscript

# =========================================================
# Multi-sample Seurat merge and Harmony integration workflow
# =========================================================
# Edit the parameters in the "User parameters" section below
# before running the script.
# =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(RColorBrewer)
  library(harmony)
  library(ggplot2)
})

# -----------------------------
# User parameters
# -----------------------------
sample_files <- list(
  sam   = "./sam_Merged_filtered.without_scale_data.rds",     # SAM sample Seurat object
  early = "./early1_early3_Merged_filtered.without_scale_data.rds",  # early-stage sample Seurat object
  mid   = "./mid_Merged_filtered.without_scale_data.rds",     # mid-stage sample Seurat object
  later = "./later_Merged_filtered.without_scale_data.rds"    # later-stage sample Seurat object
)

output_prefix <- "sam_early_mid_later_merged"   # prefix for all output files

mt_pattern  <- "^ZeamMt"    # mitochondrial gene pattern
mir_pattern <- "^zma-MIR"   # miRNA feature pattern
mp_pattern  <- "^ZeamMp"    # plastid gene pattern
mr_pattern  <- "^ZeamMr"    # ribosomal-related feature pattern

nfeatures_use <- 2000       # number of variable features
dims_use <- 1:30            # dimensions used for Harmony, neighbors, and UMAP
cluster_resolution <- 1.5   # clustering resolution
harmony_group <- "orig.ident"  # metadata column used for Harmony batch correction

marker_only_pos <- TRUE     # return only positive marker genes
marker_min_pct <- 0.25      # minimum fraction of cells expressing a marker
marker_logfc <- 0.25        # minimum log fold change threshold for marker detection

umap_width <- 8             # width of UMAP output PDF
umap_height <- 6            # height of UMAP output PDF
qc_width <- 6               # width of QC violin plot PDF
qc_height <- 4              # height of QC violin plot PDF

# -----------------------------
# Helper function
# -----------------------------
stop_if_missing <- function(cond, msg) {
  if (!cond) stop(msg, call. = FALSE)
}

simple_umap_theme <- function() {
  theme_classic() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "plain"),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}

# -----------------------------
# Check input files
# -----------------------------
message("========== Merge workflow started ==========")

for (nm in names(sample_files)) {
  stop_if_missing(
    file.exists(sample_files[[nm]]),
    paste0("Input file not found for sample '", nm, "': ", sample_files[[nm]])
  )
}

# -----------------------------
# Read Seurat objects
# -----------------------------
message("Reading Seurat objects ...")
obj_list <- lapply(sample_files, readRDS)

for (nm in names(obj_list)) {
  stop_if_missing(
    inherits(obj_list[[nm]], "Seurat"),
    paste0("Input object for sample '", nm, "' is not a Seurat object.")
  )
}

# -----------------------------
# Merge objects
# -----------------------------
message("Merging objects ...")
merged.new <- merge(
  x = obj_list[[1]],
  y = obj_list[-1]
)

# -----------------------------
# QC metrics
# -----------------------------
message("Calculating QC metrics ...")
merged.new[["percent.mt"]]  <- PercentageFeatureSet(merged.new, pattern = mt_pattern)
merged.new[["percent.mir"]] <- PercentageFeatureSet(merged.new, pattern = mir_pattern)
merged.new[["percent.mp"]]  <- PercentageFeatureSet(merged.new, pattern = mp_pattern)
merged.new[["percent.mr"]]  <- PercentageFeatureSet(merged.new, pattern = mr_pattern)

if ("orig.ident" %in% colnames(merged.new@meta.data)) {
  my_order <- names(sample_files)
  merged.new$orig.ident <- factor(merged.new$orig.ident, levels = my_order)
}

# -----------------------------
# QC plot
# -----------------------------
message("Saving QC violin plot ...")
pdf(file = paste0(output_prefix, "_QC_nFeature_RNA.pdf"), width = qc_width, height = qc_height)
print(
  VlnPlot(
    merged.new,
    features = c("nFeature_RNA"),
    ncol = 1,
    group.by = "orig.ident",
    pt.size = 0,
    y.max = 15000
  ) + theme_classic()
)
dev.off()

# -----------------------------
# Normalize / HVG / Scale / PCA
# -----------------------------
message("Running normalization and dimensionality reduction ...")
merged.new <- NormalizeData(
  merged.new,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

merged.new <- FindVariableFeatures(
  merged.new,
  selection.method = "vst",
  nfeatures = nfeatures_use
)

merged.new <- ScaleData(
  merged.new,
  features = VariableFeatures(object = merged.new)
)

merged.new <- RunPCA(
  merged.new,
  features = VariableFeatures(object = merged.new)
)

# -----------------------------
# Harmony integration and clustering
# -----------------------------
message("Running Harmony integration ...")
stop_if_missing(
  harmony_group %in% colnames(merged.new@meta.data),
  paste0("Harmony grouping column not found: ", harmony_group)
)

merged.new <- RunHarmony(
  merged.new,
  reduction = "pca",
  group.by.vars = harmony_group,
  reduction.save = "harmony"
)

merged.new <- FindNeighbors(
  merged.new,
  dims = dims_use,
  reduction = "harmony"
)

merged.new <- FindClusters(
  merged.new,
  resolution = cluster_resolution
)

merged.new <- RunUMAP(
  merged.new,
  reduction = "harmony",
  dims = dims_use,
  reduction.name = "umap"
)

# -----------------------------
# UMAP plots
# -----------------------------
message("Saving UMAP plots ...")

# plot by orig.ident if available
pdf(file = paste0(output_prefix, "_UMAP_by_orig.ident.pdf"), width = umap_width, height = umap_height)
print(
  DimPlot(
    merged.new,
    reduction = "umap",
    group.by = "orig.ident",
    label = TRUE,
    label.size = 3
  ) + simple_umap_theme()
)
dev.off()

# plot by orig.ident2 if available
if ("orig.ident2" %in% colnames(merged.new@meta.data)) {
  pdf(file = paste0(output_prefix, "_UMAP_by_orig.ident2.pdf"), width = umap_width, height = umap_height)
  print(
    DimPlot(
      merged.new,
      reduction = "umap",
      group.by = "orig.ident2",
      label = TRUE,
      label.size = 3
    ) + simple_umap_theme()
  )
  dev.off()
}

# cluster colors
cluster_ids <- levels(Idents(merged.new))
n_clusters <- length(cluster_ids)
palette_use <- colorRampPalette(c(
  "#F56867", "#FEB915", "#C798EE", "#7495D3", "#15821E",
  "#3A84E6", "#DB4C6C", "#AF5F3C", "#F9BD3F", "#ed1299",
  "#268785", "#246b93", "#d561dd", "#4aef7b", "#9ed84e",
  "#ff523f", "#03827f", "#931635", "#373bbf", "#59BE86"
))(n_clusters)

all.cells <- colnames(merged.new)

pdf(file = paste0(output_prefix, "_UMAP_by_cluster.pdf"), width = umap_width, height = umap_height)
print(
  DimPlot(
    merged.new,
    reduction = "umap",
    cells = all.cells,
    cols = palette_use,
    label = TRUE,
    label.size = 4
  ) + simple_umap_theme()
)
dev.off()

# -----------------------------
# Per-group UMAP plots
# -----------------------------
message("Saving group-specific UMAP plots ...")
for (grp in names(sample_files)) {
  grp_cells <- grep(paste0(".*_", grp, ".*"), all.cells, value = TRUE)

  if (length(grp_cells) > 0) {
    pdf(file = paste0(output_prefix, "_", grp, "_UMAP_by_cluster.pdf"),
        width = umap_width, height = umap_height)
    print(
      DimPlot(
        merged.new,
        reduction = "umap",
        cells = grp_cells,
        cols = palette_use,
        label = TRUE,
        label.size = 4
      ) + simple_umap_theme()
    )
    dev.off()
  }
}

# -----------------------------
# Marker genes
# -----------------------------
message("Running marker detection ...")
merged.new <- JoinLayers(merged.new)

merged.new.markers <- FindAllMarkers(
  merged.new,
  only.pos = marker_only_pos,
  min.pct = marker_min_pct,
  logfc.threshold = marker_logfc
)

write.csv(
  merged.new.markers,
  file = paste0(output_prefix, "_Cluster.marker.genes.csv"),
  row.names = FALSE
)

# -----------------------------
# Average expression
# -----------------------------
message("Calculating average expression ...")
aver <- as.data.frame(
  AverageExpression(
    merged.new,
    assays = "RNA",
    verbose = FALSE,
    group.by = "seurat_clusters"
  )
)

write.csv(
  aver,
  file = paste0(output_prefix, "_genes.cluster.expre.csv")
)

# -----------------------------
# Run summary
# -----------------------------
message("Saving run summary ...")
run_summary <- data.frame(
  output_prefix = output_prefix,
  sample_names = paste(names(sample_files), collapse = ";"),
  n_samples = length(sample_files),
  total_cells = ncol(merged.new),
  total_genes = nrow(merged.new),
  variable_features = nfeatures_use,
  dims_used = paste(range(dims_use), collapse = ":"),
  cluster_resolution = cluster_resolution,
  n_clusters = length(unique(Idents(merged.new)))
)

write.csv(
  run_summary,
  file = paste0(output_prefix, "_run_summary.csv"),
  row.names = FALSE
)

# cell counts by group
if ("orig.ident" %in% colnames(merged.new@meta.data)) {
  group_cell_counts <- as.data.frame(table(merged.new$orig.ident))
  colnames(group_cell_counts) <- c("group", "cell_number")
  write.csv(
    group_cell_counts,
    file = paste0(output_prefix, "_group_cell_counts.csv"),
    row.names = FALSE
  )
}

# cluster sizes
cluster_sizes <- as.data.frame(table(Idents(merged.new)))
colnames(cluster_sizes) <- c("cluster", "cell_number")
write.csv(
  cluster_sizes,
  file = paste0(output_prefix, "_cluster_cell_counts.csv"),
  row.names = FALSE
)

# -----------------------------
# Save object
# -----------------------------
message("Saving merged Seurat object ...")
merged.new@assays$RNA$scale.data <- matrix(numeric(0), nrow = 0, ncol = 0)

saveRDS(
  merged.new,
  file = paste0(output_prefix, "_Merged_filtered.without_scale_data.rds")
)

message("========== Workflow completed successfully ==========")
message("Main outputs:")
message("1. ", paste0(output_prefix, "_QC_nFeature_RNA.pdf"))
message("2. ", paste0(output_prefix, "_UMAP_by_orig.ident.pdf"))
message("3. ", paste0(output_prefix, "_UMAP_by_cluster.pdf"))
message("4. ", paste0(output_prefix, "_Cluster.marker.genes.csv"))
message("5. ", paste0(output_prefix, "_genes.cluster.expre.csv"))
message("6. ", paste0(output_prefix, "_run_summary.csv"))
message("7. ", paste0(output_prefix, "_Merged_filtered.without_scale_data.rds"))