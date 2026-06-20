# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: SPATIAL MAPPING & CHM VISUALIZATION PIPELINE ####
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────
# Description: Functions as a targeted visual inspection and cartographic mapping 
#              tool. It loads normalized point clouds and Canopy Height Models (CHMs), 
#              applies focal smoothing matrices to mitigate raster noise, and generates 
#              interactive 3D point cloud renders and publication-ready 2D spatial maps.
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
# Load required libraries for point cloud visualization and spatial mapping
library(lidR)
library(RCSF)
library(RMCC)
library(gstat)
library(sf)
library(tictoc)
library(geometry)
library(dplyr)
library(terra)
library(rgl)
library(ggplot2)
library(ggspatial)
library(tidyterra)

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration & Path Management ####
# ──────────────────────────────────────────────────────────────────────────────
# === PIPELINE TOGGLES ===
process_whole_site <- TRUE    # Set to TRUE to process the Master Site file, FALSE for individual plots
save_smoothed_output <- TRUE  # Set to TRUE to export the smoothed CHM back to the disk

# === CONFIGURE BATCH AND PLOT ===
# Change this single variable for each new batch!
date_folder <- "24. 23 April 2026"

# Define the specific plot number to visualize (Used if process_whole_site == FALSE)
Number <- 28

# --- DYNAMIC DATE EXTRACTION ---
# Extract the date part and create a safe filename format
# (e.g., "26. 08 May 2026" -> "08_May_2026")
file_date <- sub("^\\d+\\.\\s*", "", date_folder)
file_date_safe <- gsub(" ", "_", file_date)

# Base path to the remote sensing media dataset
myPath <- paste0("E:/Remote Sensing Media/", date_folder, "/")

# --- PATH CONSTRUCTION ---
if (process_whole_site) {
  # Dynamic path for site-wide processing
  name_chm <- paste0("Master_Site_CHM_Single_", file_date_safe)
  path_chm <- paste0(myPath, "07. Canopy Height Models/", name_chm, ".tif")
  
  # Set point cloud path to NULL to prevent rendering site-wide data (hardware safety)
  path_normalised <- NULL 
  
  # Set plot subtitle for cartography
  map_subtitle <- paste0("Master Site Level - ", file_date)
  
} else {
  # Dynamic paths for individual plot processing
  name_clipped <- paste0("Plot_", Number, "_", file_date_safe)
  name_classified <- paste0(name_clipped, "_classified")
  name_normalised <- paste0(name_classified, "_normalised")
  name_chm <- paste0(name_normalised, "_chm")
  
  path_normalised <- paste0(myPath, "06. Point Clouds Normalised/", name_normalised, ".las")
  path_chm <- paste0(myPath, "07. Canopy Height Models/", name_chm, ".tif")
  
  # Set plot subtitle for cartography
  map_subtitle <- paste0("Plot ", Number, " - ", file_date)
}


# ──────────────────────────────────────────────────────────────────────────────
# 3. Data Loading & Preprocessing ####
# ──────────────────────────────────────────────────────────────────────────────
# --- SAFELY LOAD & FILTER DATA ---

# 3.1 Normalised Point Cloud (Skipped if process_whole_site is TRUE)
if (!is.null(path_normalised) && file.exists(path_normalised)) {
  las_normalised <- readLAS(path_normalised)
  message("Loaded: ", name_normalised, ".las")
} else {
  if(process_whole_site) {
    message("Skipped: Point cloud loading disabled for master site to conserve memory.")
  } else {
    message("Skipped: Could not find folder or file for 06. Point Clouds Normalised")
  }
  las_normalised <- NULL
}

# 3.2 Canopy Height Model (CHM)
if (file.exists(path_chm)) {
  las_chm <- rast(path_chm)
  
  # Filter out rogue photogrammetry noise
  las_chm[las_chm > 15] <- NA  # Replace 15 with a realistic upper biological threshold
  las_chm[las_chm < 0] <- 0    # Cleans up any sub-surface negative noise
  
  # --- Smooth the CHM ---
  # Define the moving window 'w' (A 5x5 matrix is used here for smoothing)
  w <- matrix(1, nrow = 5, ncol = 5) 
  
  # Apply the focal smoothing function to the filtered CHM to remove micro-noise
  smoothed_chm <- terra::focal(las_chm, w = w, fun = max, na.rm = TRUE)
  message("Loaded and smoothed: ", name_chm, ".tif")
  
  # --- Export Smoothed CHM ---
  if (save_smoothed_output) {
    # Dynamically build output path in the same directory as the input
    out_dir <- dirname(path_chm)
    original_filename <- basename(path_chm)
    smoothed_filename <- paste0("A_Smoothed_", original_filename)
    out_path <- file.path(out_dir, smoothed_filename)
    
    # Write raster to disk
    terra::writeRaster(smoothed_chm, filename = out_path, overwrite = TRUE)
    message("Successfully exported smoothed CHM to: ", out_path)
  }
  
} else {
  message("Error: Could not find CHM file at specified path:\n", path_chm)
  las_chm <- NULL
  smoothed_chm <- NULL 
}

# ──────────────────────────────────────────────────────────────────────────────
# 4. Standard Visualization Plots ####
# ──────────────────────────────────────────────────────────────────────────────
# Plot height-normalised point cloud (Removed "axis = TRUE" to prevent rendering warnings)
if (!is.null(las_normalised)) {
  plot(las_normalised, color = "RGB", size = 2, bg = "white")
  plot(las_normalised, size = 2, bg = "white")
}

# ──────────────────────────────────────────────────────────────────────────────
# 5. Spatial Mapping (ggplot2) ####
# ──────────────────────────────────────────────────────────────────────────────
# Only generate the 2D map if the smoothed CHM was successfully created
if (!is.null(smoothed_chm)) {
  chm_map <- ggplot() +
    # Add the Smoothed Canopy Height Model raster layer
    geom_spatraster(data = smoothed_chm) +
    
    # (Optional) Add the tree locations - uncomment if PlotTrees is loaded
    # geom_sf(data = PlotTrees, color = "red", size = 1.5, shape = 16) +
    
    # Define the continuous color palette for canopy height
    scale_fill_viridis_c(
      name = "Canopy\nHeight (m)", 
      option = "viridis", 
      na.value = "transparent"
    ) +
    
    # Add the spatial scale bar
    annotation_scale(
      location = "br", 
      width_hint = 0.4, 
      text_cex = 0.8
    ) +
    
    # Add the North Arrow (Positioned inside to prevent title collisions)
    annotation_north_arrow(
      location = "tr", 
      which_north = "false",
      style = north_arrow_fancy_orienteering()
    ) +
    
    # Map formatting and dynamic titles
    labs(
      title = "Smoothed Canopy Height Model",
      subtitle = map_subtitle, # Dynamically set based on pipeline toggle
      x = "Longitude",
      y = "Latitude"
    ) +
    
    # Lock spatial coordinates
    coord_sf() +
    
    # Theme and aesthetic adjustments
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
      # Angle the longitude text so coordinates never overlap
      axis.text.x = element_text(angle = 45, hjust = 1),
      # Normalized margins to keep elements balanced
      plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "in")
    )
  
  # Render the final map
  print(chm_map)
}