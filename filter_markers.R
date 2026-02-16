#Marker Filtering per Cluster
# This script filters cluster-specific markers from
# differential expression results.


library(readxl)
library(dplyr)
library(purrr)
library(writexl)

input_file  <- "data/markers_by_cluster.xlsx"
output_file <- "results/clusters_filtered.xlsx"

sheet_names <- excel_sheets(input_file)
markers_list <- map(sheet_names, ~ read_excel(input_file, sheet = .x))
names(markers_list) <- sheet_names

markers_filtered_list <- map(markers_list, function(df) {
  df %>%
    mutate(
      p_val_adj = as.numeric(p_val_adj),
      pct.1     = as.numeric(pct.1),
      pct.2     = as.numeric(pct.2),
      avg_log2FC = as.numeric(avg_log2FC)
    ) %>%
    filter(
      p_val_adj < 0.05,
      pct.1 > 0.7,
      abs(pct.1 - pct.2) >= 0.4
    ) %>%
    arrange(desc(avg_log2FC))
})

write_xlsx(markers_filtered_list, output_file)
message("Filtrage terminé : ", output_file)