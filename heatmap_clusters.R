# Heatmap Visualization of Cluster Markers
# Visualize common and cluster-specific genes across clusters


library(readxl)
library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)

input_file <- "results/clusters_filtered_annotated.xlsx"
selected_clusters <- c("cluster1", "cluster2", "cluster4", "cluster10", "cluster13", "cluster15", "cluster20", "cluster21", "cluster24")
output_folder <- "results/heatmaps_pdf/"
output_file   <- "heatmaps_log2FC_pct1.pdf"

if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# --------------------------
# Check available sheets
# --------------------------
all_sheets <- excel_sheets(input_file)
sheets <- intersect(all_sheets, selected_clusters)
if(length(sheets) < length(selected_clusters)) warning("do not exist")

# --------------------------
# Import log2FC and pct.1
# --------------------------
list_log2FC <- lapply(sheets, function(sheet) {
  read_excel(input_file, sheet = sheet) %>%
    mutate(gene = as.character(gene)) %>%
    select(gene, avg_log2FC) %>%
    rename(value = avg_log2FC) %>%
    mutate(cluster = sheet)
})

list_pct1 <- lapply(sheets, function(sheet) {
  read_excel(input_file, sheet = sheet) %>%
    mutate(gene = as.character(gene)) %>%
    select(gene, pct.1) %>%
    rename(value = pct.1) %>%
    mutate(cluster = sheet)
})

# --------------------------
# Convert to matrices
# --------------------------
mat_log2FC <- bind_rows(list_log2FC) %>%
  pivot_wider(names_from = cluster, values_from = value) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames("gene") %>%
  as.matrix()

mat_pct1 <- bind_rows(list_pct1) %>%
  pivot_wider(names_from = cluster, values_from = value) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# --------------------------
# Define color palettes
# --------------------------
palette_log2FC <- colorRampPalette(c("#023858","white","#fd8d3c"))(100)
palette_pct1   <- colorRampPalette(c("#b58900","white","#f6e8c3"))(100)

breaks_log2FC <- seq(min(mat_log2FC), max(mat_log2FC), length.out = 101)
breaks_pct1   <- seq(0, 1, length.out = 101)

# --------------------------
# Split large matrices
# --------------------------
split_matrix <- function(mat, nb_parts = 3) {
  n <- nrow(mat)
  block_size <- ceiling(n / nb_parts)
  lapply(seq_len(nb_parts), function(i) {
    start <- (i - 1) * block_size + 1
    end <- min(i * block_size, n)
    mat[start:end, , drop = FALSE]
  })
}

nb_parts <- 3
list_mat_log2FC <- split_matrix(mat_log2FC, nb_parts)
list_mat_pct1   <- split_matrix(mat_pct1, nb_parts)

# --------------------------
# Synchronized clustering
# --------------------------
get_gene_order <- function(mat) {
  p <- pheatmap(mat, silent=TRUE, cluster_rows=TRUE, cluster_cols=TRUE)
  p$tree_row$order
}

gene_order_list <- list()
for(i in seq_len(nb_parts)){
  order_idx <- get_gene_order(list_mat_log2FC[[i]])
  list_mat_log2FC[[i]] <- list_mat_log2FC[[i]][order_idx,,drop=FALSE]
  gene_order_list[[i]] <- rownames(list_mat_log2FC[[i]])
}

for(i in seq_len(nb_parts)){
  list_mat_pct1[[i]] <- list_mat_pct1[[i]][gene_order_list[[i]],,drop=FALSE]
}

# --------------------------
# Generate PDF heatmaps
# --------------------------
pdf(file.path(output_folder, output_file), width = 15, height = 20)

for(i in seq_len(nb_parts)){
  # Heatmap log2FC
  pheatmap(list_mat_log2FC[[i]], color=palette_log2FC, breaks=breaks_log2FC,
           cluster_rows=FALSE, cluster_cols=TRUE, cellwidth=15, cellheight=15,
           fontsize_row=8, fontsize_col=10,
           main=paste("Heatmap log2FC - part", i))
  
  # Heatmap pct.1
  pheatmap(list_mat_pct1[[i]], color=palette_pct1, breaks=breaks_pct1,
           cluster_rows=FALSE, cluster_cols=TRUE, cellwidth=15, cellheight=15,
           fontsize_row=8, fontsize_col=10,
           main=paste("Heatmap pct.1 - part", i))
}

dev.off()

message("Heatmaps saved in: ", file.path(output_folder, output_file))

