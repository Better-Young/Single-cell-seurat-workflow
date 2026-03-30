#!/usr/bin/env Rscript

# =========================================================
# SoupX workflow for Seurat object
# =========================================================
# Edit the parameters in the "User parameters" section below
# before running the script.
# =========================================================

suppressPackageStartupMessages({
  library(SoupX)
  library(Seurat)
  library(Matrix)
  library(stringr)
})

# -----------------------------
# User parameters
# -----------------------------
raw_matrix_path <- "./raw_feature_bc_matrix"                         # path to raw 10X matrix directory
filtered_rds <- "./1.early_1.filtered_c3_f200.DoubletFinder_Singlet.rds"  # filtered Seurat object
output_prefix <- "2.early_1.filtered_c3_f200.DoubletFinder_Singlet.SoupX" # prefix for all output files

cluster_col <- "seurat_clusters"   # metadata column containing cluster assignments
assay_name <- "RNA"                # assay used for count matrix extraction
project_name <- "SoupX_Corrected"  # project name for corrected Seurat object

prior_rho <- 0.16                  # prior contamination fraction for SoupX
min_cells <- 3                     # minimum cells per gene for Seurat object creation
min_features <- 200                # minimum detected genes per cell for Seurat object creation

fix_barcode_suffix <- TRUE         # whether to modify barcode suffix in filtered cells
barcode_pattern <- "_early$"       # barcode suffix pattern to remove
barcode_replacement <- ""          # replacement string for barcode correction

# -----------------------------
# Helper function
# -----------------------------
stop_if_missing <- function(cond, msg) {
  if (!cond) stop(msg, call. = FALSE)
}

message("========== SoupX workflow started ==========")
message("Raw matrix path: ", raw_matrix_path)
message("Filtered Seurat object: ", filtered_rds)
message("Output prefix: ", output_prefix)

# -----------------------------
# Check input
# -----------------------------
stop_if_missing(dir.exists(raw_matrix_path),
                paste0("Raw matrix directory not found: ", raw_matrix_path))

stop_if_missing(file.exists(filtered_rds),
                paste0("Filtered Seurat object not found: ", filtered_rds))

# -----------------------------
# Read data
# -----------------------------
message("Reading raw matrix ...")
raw_counts <- Read10X(data.dir = raw_matrix_path)

message("Reading filtered Seurat object ...")
filtered_obj <- readRDS(filtered_rds)

stop_if_missing(inherits(filtered_obj, "Seurat"),
                "The filtered RDS is not a Seurat object.")

stop_if_missing(cluster_col %in% colnames(filtered_obj@meta.data),
                paste0("Metadata column '", cluster_col, "' not found in Seurat object."))

stop_if_missing(assay_name %in% names(filtered_obj@assays),
                paste0("Assay '", assay_name, "' not found in Seurat object."))

filtered_counts <- filtered_obj@assays[[assay_name]]@counts

# -----------------------------
# Fix barcode names if needed
# -----------------------------
if (fix_barcode_suffix) {
  message("Fixing barcode suffix in filtered object ...")
  colnames(filtered_obj) <- str_replace(colnames(filtered_obj), barcode_pattern, barcode_replacement)
  colnames(filtered_counts) <- str_replace(colnames(filtered_counts), barcode_pattern, barcode_replacement)
}

# -----------------------------
# Match genes
# -----------------------------
message("Matching genes between raw and filtered matrices ...")
common_genes <- intersect(rownames(raw_counts), rownames(filtered_counts))

stop_if_missing(length(common_genes) > 0,
                "No common genes found between raw and filtered matrices.")

raw_counts <- raw_counts[common_genes, ]
filtered_counts <- filtered_counts[common_genes, ]

message("Common genes: ", length(common_genes))
message("Raw droplets: ", ncol(raw_counts))
message("Filtered cells: ", ncol(filtered_counts))

# -----------------------------
# Create SoupChannel
# -----------------------------
message("Creating SoupChannel object ...")
sc <- SoupChannel(raw_counts, filtered_counts)

clusters <- filtered_obj@meta.data[[cluster_col]]
names(clusters) <- colnames(filtered_obj)

# make sure cluster vector only contains cells present in filtered_counts
clusters <- clusters[colnames(filtered_counts)]

stop_if_missing(length(clusters) == ncol(filtered_counts),
                "Cluster vector length does not match filtered matrix columns.")

sc <- setClusters(sc, clusters)

# -----------------------------
# Estimate contamination
# -----------------------------
message("Estimating contamination fraction ...")
sc <- autoEstCont(sc, priorRho = prior_rho)

rho_est <- unique(sc$metaData$rho)
message("Estimated contamination fraction (rho): ", paste(round(rho_est, 4), collapse = ", "))

# Save soup profile
write.csv(
  sc$soupProfile,
  paste0(output_prefix, "_soup_profile.csv"),
  row.names = TRUE
)

# -----------------------------
# Adjust counts
# -----------------------------
message("Adjusting counts ...")
corrected_matrix <- adjustCounts(sc, roundToInt = TRUE)

# -----------------------------
# Create corrected Seurat object
# -----------------------------
message("Creating corrected Seurat object ...")
corrected_seurat <- CreateSeuratObject(
  counts = corrected_matrix,
  project = project_name,
  min.cells = min_cells,
  min.features = min_features
)

# Try to keep metadata for overlapping cells
common_cells <- intersect(colnames(corrected_seurat), rownames(filtered_obj@meta.data))
if (length(common_cells) > 0) {
  corrected_seurat <- AddMetaData(
    corrected_seurat,
    metadata = filtered_obj@meta.data[common_cells, , drop = FALSE]
  )
}

# -----------------------------
# Save corrected object
# -----------------------------
saveRDS(corrected_seurat, file = paste0(output_prefix, ".rds"))

# -----------------------------
# Summary statistics
# -----------------------------
median_genes <- median(corrected_seurat$nFeature_RNA)
median_umis <- median(corrected_seurat$nCount_RNA)

run_summary <- data.frame(
  raw_matrix_path = raw_matrix_path,
  filtered_rds = filtered_rds,
  output_prefix = output_prefix,
  assay_name = assay_name,
  cluster_column = cluster_col,
  prior_rho = prior_rho,
  common_genes = length(common_genes),
  raw_droplets = ncol(raw_counts),
  filtered_cells_input = ncol(filtered_counts),
  corrected_cells_output = ncol(corrected_seurat),
  corrected_genes_output = nrow(corrected_seurat),
  median_nFeature_RNA = median_genes,
  median_nCount_RNA = median_umis
)

write.csv(
  run_summary,
  paste0(output_prefix, "_run_summary.csv"),
  row.names = FALSE
)

# Cell-level QC summary
qc_summary <- data.frame(
  cell_id = colnames(corrected_seurat),
  nCount_RNA = corrected_seurat$nCount_RNA,
  nFeature_RNA = corrected_seurat$nFeature_RNA
)

write.csv(
  qc_summary,
  paste0(output_prefix, "_cell_qc_summary.csv"),
  row.names = FALSE
)

message("Median genes: ", median_genes)
message("Median UMIs: ", median_umis)

message("========== Workflow completed successfully ==========")
message("Main outputs:")
message("1. ", paste0(output_prefix, ".rds"))
message("2. ", paste0(output_prefix, "_run_summary.csv"))
message("3. ", paste0(output_prefix, "_cell_qc_summary.csv"))
message("4. ", paste0(output_prefix, "_soup_profile.csv"))