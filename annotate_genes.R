# Gene Annotation
# ----------------------------------------------------------
# Annotate gene symbols using org.Hs.eg.db
############################################################

library(readxl)
library(dplyr)
library(purrr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(writexl)

input_file  <- "results/clusters_filtered.xlsx"
output_file <- "results/clusters_filtered_annotated.xlsx"

sheet_names <- excel_sheets(input_file)
clusters_list <- map(sheet_names, ~ read_excel(input_file, sheet = .x))
names(clusters_list) <- sheet_names

clusters_list <- map(clusters_list, ~ {
  .x %>%
    select(any_of(c("gene", "avg_log2FC","pct.1", "pct.2","p_val_adj"))) %>%
    mutate(gene = as.character(gene))
})

gene_symbols <- clusters_list %>%
  bind_rows() %>%
  pull(gene) %>%
  unique()

mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = gene_symbols,
  columns  = c("SYMBOL", "GENENAME"),
  keytype  = "SYMBOL"
) %>% distinct(SYMBOL, .keep_all = TRUE)

clusters_annotated <- map(clusters_list, ~ {
  left_join(.x, mapping, by = c("gene" = "SYMBOL"))
})

write_xlsx(clusters_annotated, output_file)
message("Annotation terminée : ", output_file)
