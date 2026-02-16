# Marker validation in scRNA-seq clusters

## Project description

This analysis was performed to validate the endothelial identity of selected scRNA-seq clusters and to distinguish them from fibroblast populations.

The script evaluates:
- Canonical endothelial markers (CDH5, PECAM1, VWF, etc.)
- LSEC-specific markers (CLEC4G, STAB2, LYVE1, etc.)
- Fibroblast markers (COL1A1, COL3A1, ACTA2, etc.)

## Method

For each cluster:
- avg_log2FC
- pct.1 - pct.2 difference
- adjusted p-value

Genes are classified as:
- Significant upregulated
- Significant downregulated
- Non-significant upregulated
- Non-significant downregulated

Only genes with pct difference ≥ 0.1 are displayed.

## Pipeline to validate endothelial vs fibroblast clusters in scRNA-seq data.

### Steps:

Filter markers per cluster based on significance and expression (filter_markers.R).

Annotate genes with full names (annotate_genes.R).

Visualize expression with heatmaps for log2FC and pct.1
