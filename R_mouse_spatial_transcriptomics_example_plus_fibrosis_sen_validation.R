# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Define file paths
vehicle_path <- "C:/Users/gabriel.mecalaguna/Documents/processed_mouse_seurat_objects/d21_Veh_2_V10A20-052-A1.rds"
bleomycin_path <- "C:/Users/gabriel.mecalaguna/Documents/processed_mouse_seurat_objects/d21_BLM_8_V10A20-053-C1.rds"

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


# Get expression data for Cd3e and Trdc
expr_data <- FetchData(merged_obj, vars = c("Cd3e", "Trdc", "Prf1", "Cdkn1a"))

# Create γδ T cell classification
merged_obj$GD_T_cell <- ifelse(
  expr_data[["Trdc"]] > 0 & expr_data[["Cd3e"]] > 0 & expr_data[["Prf1"]] > 0, 
  "γδ T cells", 
  "Other"
)

# Plot γδ T cells spatial distribution

# Plot γδ T cells spatial distribution with custom legend title
p_vehicle_prf <- SpatialPlot(
  object = merged_obj,
  group.by = "GD_T_cell",
  images = "Vehicle",
  pt.size = 2.5,
  cols = c("γδ T cells" = "red", "Other" = "#E0E0E0"),
  alpha = c("γδ T cells" = 0.7)
) + 
  labs(fill = "γδ T cells")

p_vehicle_prf

# Save as TIFF with 600 DPI
#ggsave(
 # filename = "vehicle_prf_gamma_delta_T_cells.tiff",
#  plot = p_vehicle_prf,
#  device = "tiff",
#  dpi = 600,
#  width = 8,      # width in inches
#  height = 6     # height in inches
  
#)

# Plot γδ T cells spatial distribution with custom legend title
p_blm_prf <- SpatialPlot(
  object = merged_obj,
  group.by = "GD_T_cell",
  images = "Bleomycin",
  pt.size = 2.5,
  cols = c("γδ T cells" = "red", "Other" = "#E0E0E0"),
  alpha = c("γδ T cells" = 0.7)
) + 
  labs(fill = "γδ T cells")

p_blm_prf

# Save as TIFF with 600 DPI
#ggsave(
#  filename = "blm_prf_gamma_delta_T_cells.tiff",
#  plot = p_blm_prf,
#  device = "tiff",
#  dpi = 600,
#  width = 8,      # width in inches
#  height = 6     # height in inches
  
#)


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
