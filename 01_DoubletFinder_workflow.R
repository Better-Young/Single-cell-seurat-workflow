#!/usr/bin/env Rscript

# =========================================================
# DoubletFinder workflow for Seurat object
# =========================================================
# Edit the parameters in the "User parameters" section below
# before running the script.
# =========================================================

suppressPackageStartupMessages({
  library(DoubletFinder)
  library(Seurat)
  library(ggplot2)
})

# -----------------------------
# User parameters
# -----------------------------
input_rds <- "0.early_1.filtered_c3_f200.rds"   # input Seurat object
output_prefix <- "1.early_1.filtered_c3_f200"   # prefix for output files

cluster_col <- "seurat_clusters"   # metadata column used for homotypic modeling
pcs_use <- 1:20                    # PCs used in DoubletFinder
pN_use <- 0.25                     # recommended default for DoubletFinder
doublet_rate <- 0.01               # expected doublet rate
sct_use <- FALSE                   # TRUE if SCTransform workflow was used
remove_scale_data <- TRUE          # remove scale.data before saving to reduce file size

# -----------------------------
# Helper function
# -----------------------------
stop_if_missing <- function(cond, msg) {
  if (!cond) stop(msg, call. = FALSE)
}

message("========== DoubletFinder workflow started ==========")
message("Input file: ", input_rds)
message("Output prefix: ", output_prefix)

# -----------------------------
# Check input
# -----------------------------
stop_if_missing(file.exists(input_rds),
                paste0("Input file not found: ", input_rds))

object <- readRDS(input_rds)

stop_if_missing(inherits(object, "Seurat"),
                "The input RDS is not a Seurat object.")

stop_if_missing("pca" %in% names(object@reductions),
                "PCA reduction not found in Seurat object. Please run RunPCA first.")

stop_if_missing("umap" %in% names(object@reductions),
                "UMAP reduction not found in Seurat object. Please run RunUMAP first.")

stop_if_missing(cluster_col %in% colnames(object@meta.data),
                paste0("Metadata column '", cluster_col, "' not found in Seurat object."))

n_cells <- ncol(object)
n_genes <- nrow(object)

message("Cells: ", n_cells)
message("Genes: ", n_genes)

# -----------------------------
# Parameter sweep
# -----------------------------
message("Running paramSweep ...")
sweep.list <- paramSweep(object, sct = sct_use)

message("Summarizing sweep results ...")
sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)

message("Finding best pK ...")
bcmvn <- find.pK(sweep.stats)

stop_if_missing(nrow(bcmvn) > 0, "find.pK returned no results.")

pK_bcmvn <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
stop_if_missing(!is.na(pK_bcmvn), "Failed to determine optimal pK.")

message("Best pK = ", pK_bcmvn)

write.csv(sweep.stats,
          paste0(output_prefix, "_paramSweep_stats.csv"),
          row.names = FALSE)
write.csv(bcmvn,
          paste0(output_prefix, "_pK_metrics.csv"),
          row.names = FALSE)

# -----------------------------
# Estimate expected doublets
# -----------------------------
message("Estimating homotypic doublet proportion ...")
homotypic.prop <- modelHomotypic(as.vector(object@meta.data[[cluster_col]]))

nExp_poi <- round(doublet_rate * n_cells)
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

message("Estimated doublets (raw): ", nExp_poi)
message("Estimated doublets (adjusted): ", nExp_poi.adj)
message("Homotypic proportion: ", round(homotypic.prop, 4))

# -----------------------------
# Run DoubletFinder
# -----------------------------
message("Running DoubletFinder ...")
object <- doubletFinder(
  object,
  PCs = pcs_use,
  pN = pN_use,
  pK = pK_bcmvn,
  nExp = nExp_poi.adj,
  sct = sct_use
)

group_id <- paste("DF.classifications", pN_use, pK_bcmvn, nExp_poi.adj, sep = "_")

stop_if_missing(group_id %in% colnames(object@meta.data),
                paste0("Doublet classification column not found: ", group_id))

# -----------------------------
# Save metadata
# -----------------------------
meta_out <- object@meta.data
meta_out$cell_id <- rownames(meta_out)
write.csv(meta_out,
          paste0(output_prefix, "_metadata_with_doubletfinder.csv"),
          row.names = FALSE)

# -----------------------------
# Plot UMAP
# -----------------------------
message("Saving UMAP plot ...")
pdf(file = paste0(output_prefix, "_DoubletFinder_UMAP.pdf"), width = 6, height = 5)
print(
  DimPlot(object, reduction = "umap", group.by = group_id, raster = FALSE) +
    theme_classic() +
    labs(title = "DoubletFinder classification")
)
dev.off()

# -----------------------------
# Subset singlets
# -----------------------------
message("Subsetting singlets ...")
singlet_cells <- rownames(object@meta.data)[object@meta.data[[group_id]] == "Singlet"]

stop_if_missing(length(singlet_cells) > 0,
                "No singlet cells detected. Please check the DoubletFinder results.")

object.singlet <- subset(object, cells = singlet_cells)

if (remove_scale_data) {
  assay_name <- DefaultAssay(object.singlet)
  object.singlet@assays[[assay_name]]@scale.data <- matrix(numeric(0), nrow = 0, ncol = 0)
}

saveRDS(object.singlet, file = paste0(output_prefix, "_singlet.rds"))

# -----------------------------
# Summary table
# -----------------------------
classification_table <- as.data.frame(table(object@meta.data[[group_id]]))
colnames(classification_table) <- c("classification", "cell_number")
classification_table$fraction <- classification_table$cell_number / sum(classification_table$cell_number)

write.csv(classification_table,
          paste0(output_prefix, "_classification_summary.csv"),
          row.names = FALSE)

run_summary <- data.frame(
  input_file = input_rds,
  output_prefix = output_prefix,
  cells_total = n_cells,
  genes_total = n_genes,
  cluster_column = cluster_col,
  pcs_used = paste(range(pcs_use), collapse = ":"),
  pN = pN_use,
  best_pK = pK_bcmvn,
  doublet_rate = doublet_rate,
  homotypic_proportion = homotypic.prop,
  expected_doublets_raw = nExp_poi,
  expected_doublets_adjusted = nExp_poi.adj,
  singlet_cells = sum(object@meta.data[[group_id]] == "Singlet"),
  doublet_cells = sum(object@meta.data[[group_id]] == "Doublet")
)

write.csv(run_summary,
          paste0(output_prefix, "_run_summary.csv"),
          row.names = FALSE)

message("========== Workflow completed successfully ==========")
message("Main outputs:")
message("1. ", paste0(output_prefix, "_DoubletFinder_UMAP.pdf"))
message("2. ", paste0(output_prefix, "_metadata_with_doubletfinder.csv"))
message("3. ", paste0(output_prefix, "_classification_summary.csv"))
message("4. ", paste0(output_prefix, "_run_summary.csv"))
message("5. ", paste0(output_prefix, "_singlet.rds"))