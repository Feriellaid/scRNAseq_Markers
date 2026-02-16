############################################################
# Marker Expression Validation in scRNA-seq Clusters
# NOTE:
# - Cluster identities are anonymized.
# - The original dataset is not included in this repository.
# - The input file must contain one sheet per cluster.
#
############################################################


library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(pheatmap)



# Path to Excel file 
file_path <- "data/markers_by_cluster.xlsx"

# Thresholds
seuil_diff <- 0.1
seuil_p    <- 0.05


# ==========================================================
# Load cluster sheets
# ==========================================================

# The Excel file must contain one sheet per cluster.
# Sheet names are assumed to represent anonymized clusters.

sheet_names <- excel_sheets(file_path)

list_clusters <- lapply(sheet_names, function(sheet) {
  read_excel(file_path, sheet = sheet)
})

names(list_clusters) <- sheet_names


# ==========================================================
# Gene alias mapping
# ==========================================================

gene_map <- tribble(
  ~ref, ~alias,
  
  # Endothelial markers
  "CDH5","CDH5",
  "CDH5","CD144",
  
  "PECAM1","PECAM1",
  "PECAM1","CD31",
  
  "VWF","VWF",
  
  "ERG","ERG",
  "ETS1","ETS1",
  "ETS2","ETS2",
  
  "KLF2","KLF2",
  
  "VCAM1","VCAM1",
  "ICAM1","ICAM1",
  
  "MCAM","MCAM",
  "MCAM","CD146",
  
  "PLVAP","PLVAP",
  
  # LSEC-associated markers
  "CLEC4G","CLEC4G",
  "STAB2","STAB2",
  "LYVE1","LYVE1",
  "FCGR2B","FCGR2B",
  "OIT3","OIT3",
  
  # Fibroblast / mesenchymal markers
  "ACTA2","ACTA2",
  "VIM","VIM",
  "COL1A1","COL1A1",
  "COL3A1","COL3A1",
  "COL4A1","COL4A1",
  "S100A4","S100A4",
  
  # Additional markers
  "CXCR4","CXCR4",
  "GATA4","GATA4",
  "IL33","IL33",
  "CD34","CD34",
  "ACKR1","ACKR1",
  "ACKR3","ACKR3",
  "SEMA3A","SEMA3A",
  "SMAD3","SMAD3"
)


# ==========================================================
# Data cleaning
# ==========================================================

list_clusters_clean <- lapply(list_clusters, function(df) {
  
  df %>%
    select(gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
    mutate(
      gene       = as.character(gene),
      avg_log2FC = as.numeric(avg_log2FC),
      pct.1      = as.numeric(pct.1),
      pct.2      = as.numeric(pct.2),
      p_val_adj  = as.numeric(p_val_adj)
    )
})


# ==========================================================
# Merge clusters and map gene aliases
# ==========================================================

expr_long <- map2_dfr(
  list_clusters_clean,
  names(list_clusters_clean),
  ~ .x %>%
    inner_join(gene_map, by = c("gene" = "alias")) %>%
    mutate(cluster = .y)
)

genes_ref <- unique(gene_map$ref)

expr_state <- expr_long %>%
  group_by(ref, cluster) %>%
  summarise(
    avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
    pct1       = mean(pct.1, na.rm = TRUE),
    pct2       = mean(pct.2, na.rm = TRUE),
    p_adj_min  = min(p_val_adj, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  complete(ref = genes_ref, cluster = names(list_clusters)) %>%
  mutate(
    diff_pct = pct1 - pct2,
    etat = case_when(
      is.na(avg_log2FC) |
        is.na(diff_pct) |
        is.na(p_adj_min) |
        diff_pct < seuil_diff ~ "NA",
      
      p_adj_min < seuil_p & avg_log2FC >= 0 ~ "sig_up",
      p_adj_min < seuil_p & avg_log2FC <  0 ~ "sig_down",
      p_adj_min >= seuil_p & avg_log2FC >= 0 ~ "nonsig_up",
      p_adj_min >= seuil_p & avg_log2FC <  0 ~ "nonsig_down"
    )
  )


# ==========================================================
# Matrix for heatmap categories
# ==========================================================

mat_etat <- expr_state %>%
  mutate(etat_num = case_when(
    etat == "nonsig_down" ~ 0,
    etat == "nonsig_up"   ~ 1,
    etat == "sig_down"    ~ 2,
    etat == "sig_up"      ~ 3,
    TRUE                  ~ NA_real_
  )) %>%
  select(ref, cluster, etat_num) %>%
  pivot_wider(names_from = cluster, values_from = etat_num) %>%
  as.data.frame()

rownames(mat_etat) <- mat_etat$ref
mat_etat$ref <- NULL


# ==========================================================
# Matrix for displayed log2FC values
# ==========================================================

mat_values <- expr_state %>%
  mutate(val_display = ifelse(etat == "NA", NA, avg_log2FC)) %>%
  select(ref, cluster, val_display) %>%
  pivot_wider(names_from = cluster, values_from = val_display) %>%
  as.data.frame()

rownames(mat_values) <- mat_values$ref
mat_values$ref <- NULL

mat_numbers <- apply(mat_values, 2, function(col) {
  ifelse(is.na(col), "", sprintf("%.2f", col))
})

rownames(mat_numbers) <- rownames(mat_values)


# ==========================================================
# Heatmap visualization
# ==========================================================

tes_couleurs <- c(
  "#FFFFCC",   # nonsig_down
  "#90EE90",   # nonsig_up
  "#FFD580",   # sig_down
  "#FF9999"    # sig_up
)

pheatmap(
  mat_etat,
  color           = tes_couleurs,
  breaks          = c(-0.5, 0.5, 1.5, 2.5, 3.5),
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  legend_breaks   = c(0, 1, 2, 3),
  legend_labels   = c("Non-significant down",
                      "Non-significant up",
                      "Significant down",
                      "Significant up"),
  main            = "Marker expression across anonymized clusters",
  fontsize_row    = 9,
  angle_col       = 45,
  na_col          = "white",
  display_numbers = mat_numbers,
  number_color    = "black"
)


writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
