# LOAD REQUIRED LIBRARIES
library(Seurat)
library(dplyr)

# Data retrieved from: https://www.ebi.ac.uk/biostudies/studies/S-BSST1410

# LOAD THE DATASET
mouse_data <- readRDS("C:/Users/gabriel.mecalaguna/Downloads/mm_visium_all_stutility_obj.rds")

# FUNCTION TO COUNT POSITIVE SPOTS
count_positive_spots <- function(seurat_obj, gene_name) {
  gene_data <- FetchData(object = seurat_obj, vars = gene_name)
  seurat_obj[[paste0(gene_name, "_expression")]] <- gene_data[[gene_name]]
  
  gene_by_sample <- seurat_obj@meta.data %>%
    mutate(
      gene_expr = .data[[paste0(gene_name, "_expression")]],
      group = case_when(
        condition == "control" ~ "Control",
        condition == "bleomycin" & day == "d7" ~ "Day 7 - BLM",
        condition == "bleomycin" & day == "d21" ~ "Day 21 - BLM"
      )
    ) %>%
    group_by(sample_name, group) %>%
    summarise(positive_spots = sum(gene_expr > 0)) %>%
    ungroup() %>%
    mutate(gene = gene_name)
  
  return(gene_by_sample)
}

# ANALYZE MULTIPLE GENES
genes_to_analyze <- c("Cdkn2a", "Cdkn1a", "Raet1e", "H60c", "Trdc", "Ulbp1", "Krt8")

results_list <- lapply(genes_to_analyze, function(gene) count_positive_spots(mouse_data, gene))
all_results <- bind_rows(results_list)

# VIEW RESULTS
all_results

openxlsx::write.xlsx(all_results, "gene_positive_spots_results.xlsx")
