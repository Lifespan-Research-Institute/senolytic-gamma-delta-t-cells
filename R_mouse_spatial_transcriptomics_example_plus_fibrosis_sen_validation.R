# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Define file paths
vehicle_path <- "./processed_mouse_seurat_objects/d21_Veh_2_V10A20-052-A1.rds"
bleomycin_path <- "./processed_mouse_seurat_objects/d21_BLM_8_V10A20-053-C1.rds"

# Load Seurat objects
cat("Loading Vehicle sample...\n")
vehicle_obj <- readRDS(vehicle_path)
cat("Loading Bleomycin sample...\n")
bleomycin_obj <- readRDS(bleomycin_path)

# Add sample identifiers before merging
vehicle_obj$sample <- "d21_Veh_2_V10A20-052-A1"
vehicle_obj$condition <- "vehicle"
vehicle_obj$day <- "day21"

bleomycin_obj$sample <- "d21_BLM_8_V10A20-053-C1"
bleomycin_obj$condition <- "bleomycin"
bleomycin_obj$day <- "day21"

# Merge objects (vehicle first, then bleomycin)
cat("Merging samples...\n")
merged_obj <- merge(vehicle_obj, bleomycin_obj, 
                    add.cell.ids = c("Vehicle", "Bleomycin"),
                    project = "SpatialTranscriptomics")

# Rename the slices to more meaningful names
cat("Renaming slices...\n")
current_images <- names(merged_obj@images)
cat("Current slice names:", paste(current_images, collapse = ", "), "\n")

# Rename slices to Vehicle and Bleomycin
names(merged_obj@images) <- c("Vehicle", "Bleomycin")
cat("New slice names:", paste(names(merged_obj@images), collapse = ", "), "\n")

# Check merged object
cat("Merged object summary:\n")
print(merged_obj)
cat("Number of spots per sample:\n")
table(merged_obj$sample)

# Define all γδ T cell markers
genes <- c("Cd3e", "Trdc", "Tcrg-C1", "Tcrg-C2", "Tcrg-C4", "Prf1")
available_genes <- genes[genes %in% rownames(merged_obj)]

cat("Available genes:", paste(available_genes, collapse = ", "), "\n")

# Check if essential Cd3e is present
if (!"Cd3e" %in% available_genes) {
  stop("Cd3e is not available in the dataset")
}

# Get expression data for all available markers
expr_data <- FetchData(merged_obj, vars = available_genes)

# Define CD3+ spots
cd3_positive <- expr_data[["Cd3e"]] > 0

# Check for gamma-delta TCR markers
gd_markers <- c("Trdc", "Tcrg-C1", "Tcrg-C2", "Tcrg-C4")
available_gd_markers <- gd_markers[gd_markers %in% available_genes]

cat("Available γδ TCR markers:", paste(available_gd_markers, collapse = ", "), "\n")

# Create logical OR across all available gamma-delta markers
gd_marker_positive <- rep(FALSE, nrow(expr_data))
for (marker in available_gd_markers) {
  gd_marker_positive <- gd_marker_positive | (expr_data[[marker]] > 0)
}

# Define γδ T cells: Cd3e+ AND (at least one γδ TCR marker+)
gamma_delta_positive <- cd3_positive & gd_marker_positive

# Create classification for all γδ T cells
merged_obj$GD_T_cell <- ifelse(
  gamma_delta_positive,
  "γδ T cells",
  "Other"
)

# Count γδ T cells per sample
cat("\nγδ T cell counts (Cd3e+ AND any γδ TCR marker+):\n")
table(merged_obj$GD_T_cell, merged_obj$sample)

# If Prf1 is available, also create cytotoxic γδ T cell classification
if ("Prf1" %in% available_genes) {
  prf1_positive <- expr_data[["Prf1"]] > 0
  cytotoxic_gamma_delta_positive <- gamma_delta_positive & prf1_positive
  
  merged_obj$Cytotoxic_GD_T_cell <- ifelse(
    cytotoxic_gamma_delta_positive,
    "Cytotoxic γδ T cells",
    "Other"
  )
  
  cat("\nCytotoxic γδ T cell counts (γδ T cells AND Prf1+):\n")
  table(merged_obj$Cytotoxic_GD_T_cell, merged_obj$sample)
}

# Plot γδ T cells spatial distribution for Vehicle
p_vehicle_gd <- SpatialPlot(
  object = merged_obj,
  group.by = "GD_T_cell",
  images = "Vehicle",
  pt.size = 2.5,
  cols = c("γδ T cells" = "red", "Other" = "#E0E0E0"),
  alpha = c("γδ T cells" = 0.7)
) + 
  labs(fill = "γδ T cells")
p_vehicle_gd

# Plot γδ T cells spatial distribution for Bleomycin
p_blm_gd <- SpatialPlot(
  object = merged_obj,
  group.by = "GD_T_cell",
  images = "Bleomycin",
  pt.size = 2.5,
  cols = c("γδ T cells" = "red", "Other" = "#E0E0E0"),
  alpha = c("γδ T cells" = 0.7)
) + 
  labs(fill = "γδ T cells")
p_blm_gd

# If Prf1 is available, also plot cytotoxic γδ T cells
if ("Prf1" %in% available_genes) {
  # Plot cytotoxic γδ T cells for Vehicle
  p_vehicle_cytotoxic <- SpatialPlot(
    object = merged_obj,
    group.by = "Cytotoxic_GD_T_cell",
    images = "Vehicle",
    pt.size = 2.5,
    cols = c("Cytotoxic γδ T cells" = "red", "Other" = "#E0E0E0"),
    alpha = c("Cytotoxic γδ T cells" = 0.7)
  ) + 
    labs(fill = "Cytotoxic γδ T cells")
  p_vehicle_cytotoxic
  
  # Plot cytotoxic γδ T cells for Bleomycin
  p_blm_cytotoxic <- SpatialPlot(
    object = merged_obj,
    group.by = "Cytotoxic_GD_T_cell",
    images = "Bleomycin",
    pt.size = 2.5,
    cols = c("Cytotoxic γδ T cells" = "red", "Other" = "#E0E0E0"),
    alpha = c("Cytotoxic γδ T cells" = 0.7)
  ) + 
    labs(fill = "Cytotoxic γδ T cells")
  p_blm_cytotoxic
}

# Save plots (uncomment to save)
# ggsave(
#   filename = "vehicle_gamma_delta_T_cells.tiff",
#   plot = p_vehicle_gd,
#   device = "tiff",
#   dpi = 600,
#   width = 8,
#   height = 6
# )

# ggsave(
#   filename = "blm_gamma_delta_T_cells.tiff",
#   plot = p_blm_gd,
#   device = "tiff",
#   dpi = 600,
#   width = 8,
#   height = 6
# )

# Plot senescence and fibrosis markers
SpatialFeaturePlot(
  object = merged_obj,
  features = "Cdkn1a",
  pt.size = 2.5,
  alpha = c(0.7, 1)
)

SpatialFeaturePlot(
  object = merged_obj,
  features = "Cdkn2a",
  pt.size = 2.5,
  alpha = c(0.7, 1)
)

SpatialFeaturePlot(
  object = merged_obj,
  features = "Col1a1",
  pt.size = 2.5,
  alpha = c(0.7, 1)
)

# Define the features you want to plot
features <- c("Cdkn1a", "Cdkn2a", "Col1a1")

# Loop through each feature
for (feature in features) {
  # Create the spatial feature plot
  p <- SpatialFeaturePlot(
    object = merged_obj,
    features = feature,
    pt.size = 2.5,
    alpha = c(0.7, 1)
  )
  
  # Save as TIFF with 600 DPI
  ggsave(
    filename = paste0(feature, "_spatial_plot.tiff"),
    plot = p,
    device = "tiff",
    dpi = 600,
    width = 8,
    height = 6
  )
  
  # Optional: print progress
  cat("Saved:", feature, "_spatial_plot.tiff\n")
}

