# single-cell-seurat-workflow

A collection of Seurat workflows for single-cell RNA-seq data analysis, including:

- DoubletFinder-based doublet detection
- SoupX-based ambient RNA correction
- Multi-sample Seurat merge and Harmony integration

## Files

- `01_DoubletFinder_workflow.R`: doublet detection workflow for Seurat objects
- `02_SoupX_workflow.R`: ambient RNA correction workflow using SoupX
- `03_merge_harmony_workflow.R`: multi-sample merge, Harmony integration, clustering, and marker detection workflow

## Notes

Please edit the **User parameters** section in each script before running.

## Requirements

Main R packages used in this project include:

- Seurat
- DoubletFinder
- SoupX
- harmony
- ggplot2
- dplyr
- patchwork
- RColorBrewer
- Matrix
- stringr
