# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: CHM HEIGHT EXTRACTION & VISUALIZATION PIPELINE ####
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
# Load required libraries for point cloud processing, spatial operations, and plotting
library(lidR)
library(RCSF)
library(RMCC)
library(gstat)
library(sf)
library(tictoc)
library(geometry)
library(dplyr)
library(future)
library(sp)
library(terra)
library(exactextractr) 

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration & Path Management ####
# ──────────────────────────────────────────────────────────────────────────────
# === CONFIGURE PATHS ===
# Change this single variable for each new batch!
date_folder <- "03. 30 October 2025"

# Extract the date part and create a safe filename format
# (e.g., "17. 02 March 2026" -> "02_March_2026")
file_date <- sub("^\\d+\\.\\s*", "", date_folder)
file_date_safe <- gsub(" ", "_", file_date)

# ──────────────────────────────────────────────────────────────────────────────
# 3. Spatial Data Loading & Preprocessing ####
# ──────────────────────────────────────────────────────────────────────────────
# --- LOAD TREE POLYGONS ---
path_trees <- paste0("E:/Remote Sensing Media/", date_folder, "/08. Crown Polygons/Crown_Polygons_", file_date_safe, ".shp")
trees <- st_read(path_trees)

# Clean tree geometry dataset by removing unnecessary metadata columns
# This prevents duplicate name errors and ESRI shapefile driver issues during export
trees <- trees %>%
  select(-any_of(c("grop_ld", "group_ulid", "N_GM", "id", "N_FG", "N_BG", "BBox")))

# --- LOAD CANOPY HEIGHT MODEL (CHM) ---
# Read in the pre-processed and clamped Canopy Height Model
ctg_chm <- rast(paste0("E:/Remote Sensing Media/", date_folder, "/07. Canopy Height Models/Master_Site_CHM_Single_", file_date_safe, ".tif"))

# # --- LOAD CANOPY HEIGHT MODEL (CHM) ---
# # Read in the pre-processed Canopy Height Model (using the virtual raster reference)
# ctg_chm <- rast(paste0("E:/Remote Sensing Media/", date_folder, "/07. Canopy Height Models/rasterize_canopy.vrt"))

# ──────────────────────────────────────────────────────────────────────────────
# 4. CHM & Species Visualization ####
# ──────────────────────────────────────────────────────────────────────────────
# Define a consistent color palette mapped to specific Eucalyptus species
species_colors <- c(
  "Cladocalyx"    = "#336998",
  "Grandis"       = "#97dde3",
  "Cloeziana"     = "#ffffff",
  "Urophylla"     = "#e3acff",
  "Grandis clone" = "#ff7da0"
)

# Plot the base Canopy Height Model with a capped height range for better contrast
plot(ctg_chm, range = c(0, 5), main = "Canopy Height Model with Tree Species")

# Overlay the individual tree crown polygons, mapped to their respective species colors
# Note: 'border = "black"' adds a black outline ensuring white/light polygons remain visible
plot(trees$geometry, 
     add = TRUE, 
     col = species_colors[trees$Species])

# ──────────────────────────────────────────────────────────────────────────────
# 5. Metric Extraction & Export ####
# ──────────────────────────────────────────────────────────────────────────────
tic()

# Calculate exact maximum tree height within each delineated crown polygon
# (exact_extract outputs directly as a vector, simplifying alignment)
trees$Tree_Height <- exact_extract(ctg_chm, trees, 'max')

# Define final output paths for the metrics
output_path_shp <- paste0("E:/Remote Sensing Media/", date_folder, "/09. Crown Metrics/Crown_Metrics_", file_date_safe, ".shp")
output_path_csv <- paste0("E:/Remote Sensing Media/", date_folder, "/09. Crown Metrics/Crown_Metrics_", file_date_safe, ".csv")

# Export updated geospatial shapefile (overwriting old files to prevent schema errors)
st_write(trees, output_path_shp, delete_dsn = TRUE)

# Export lightweight tabular CSV data (drops the messy spatial geometry text)
write.csv(st_drop_geometry(trees), output_path_csv, row.names = FALSE)

message("✅ Extraction complete. Metrics saved to: ", output_path_csv)
toc()