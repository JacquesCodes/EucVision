# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: DATA VISUALIZATION & ITD (INDIVIDUAL TREE DETECTION) PIPELINE ####
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────
# Description: Serves as an interactive visual inspection and diagnostic tool 
#              for the broader point cloud processing pipeline. It selectively 
#              loads specific temporal plots to render multi-stage 3D point cloud 
#              layers (raw, ground-classified, height-normalized) alongside 2D 
#              Canopy Height Models. It also features modular, optional routines 
#              for Individual Tree Detection (ITD), 3D Dalponte canopy segmentation, 
#              and dynamic 3D rgl animation exports.
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
# Load required libraries for 3D point cloud visualization, spatial processing, and raster manipulation
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

# (Reverted 3D plotting back to default external rgl window)

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration & Path Management ####
# ──────────────────────────────────────────────────────────────────────────────
# === CONFIGURE BATCH AND PLOT ===
# Change this single variable for each new batch!
date_folder <- "20. 23 March 2026"

# Define the specific plot number to visualize
Number <- "28"

# Extract the date part and create a safe filename format
# (e.g., "17. 02 March 2026" -> "02_March_2026")
file_date <- sub("^\\d+\\.\\s*", "", date_folder)
file_date_safe <- gsub(" ", "_", file_date)

# Base path to the remote sensing media dataset
myPath <- paste0("E:/Remote Sensing Media/",date_folder,"/")

# Dynamically construct the file prefixes based on the targeted plot and date
name_clipped <- paste0("Plot_", Number, "_", file_date_safe)
name_classified <- paste0(name_clipped, "_classified")
name_normalised <- paste0(name_classified, "_normalised")
name_chm <- paste0(name_normalised, "_chm")

# Build the full absolute file paths
path_clipped <- paste0(myPath, "04. Point Clouds Clipped/", name_clipped, ".las")
path_classified <- paste0(myPath, "05. Point Clouds Ground Classified/", name_classified, ".las")
path_normalised <- paste0(myPath, "06. Point Clouds Normalised/", name_normalised, ".las")
path_chm <- paste0(myPath, "07. Canopy Height Models/", name_chm, ".tif")

# ──────────────────────────────────────────────────────────────────────────────
# 3. Data Loading & Preprocessing ####
# ──────────────────────────────────────────────────────────────────────────────
# --- SAFELY LOAD POINT CLOUDS & RASTERS ---

# 3.1 Cropped Point Cloud
if (file.exists(path_clipped)) {
  las <- readLAS(path_clipped)
  message("Loaded: ", name_clipped, ".las")
} else {
  message("Skipped: Could not find folder or file for 04. Point Clouds Clipped")
  las <- NULL
}

# 3.2 Classified Point Cloud
if (file.exists(path_classified)) {
  las_classified <- readLAS(path_classified)
  message("Loaded: ", name_classified, ".las")
} else {
  message("Skipped: Could not find folder or file for 05. Point Clouds Ground Classified")
  las_classified <- NULL
}

# 3.3 Normalised Point Cloud
if (file.exists(path_normalised)) {
  las_normalised <- readLAS(path_normalised)
  message("Loaded: ", name_normalised, ".las")
} else {
  message("Skipped: Could not find folder or file for 06. Point Clouds Normalised")
  las_normalised <- NULL
}

# 3.4 Canopy Height Model (CHM)
if (file.exists(path_chm)) {
  las_chm <- rast(path_chm)
  
  # Filter out rogue photogrammetry noise to ensure clean visualization
  las_chm[las_chm > 15] <- NA  # Replace 15 with a realistic upper biological threshold
  las_chm[las_chm < 0] <- 0    # Cleans up any sub-surface negative noise
  
  message("Loaded: ", name_chm, ".tif")
} else {
  message("Skipped: Could not find folder or file for 07. Canopy Height Models")
  las_chm <- NULL
}

# --- SAFELY LOAD & TRANSFORM SHAPEFILES ---
# Load and transform the Crown Polygons shapefile
path_trees <- paste0(myPath, "08. Crown Polygons/Crown_Polygons_", file_date_safe, ".shp")

if (file.exists(path_trees)) {
  trees <- st_read(path_trees, quiet = TRUE)
  
  # The CRS Fix: Assign raster metadata to the shapefile
  if (!is.null(las_chm)) {
    st_crs(trees) <- st_crs(las_chm)
  }
  
  # Filter for just the target plot we are currently visualizing
  PlotTrees <- trees[trees$Plot == Number | trees$Plot == as.numeric(Number),]
  message("Loaded and filtered Crown Polygons for Plot ", Number)
  
} else {
  message("Skipped: Could not find folder or file for 08. Crown Polygons")
  PlotTrees <- NULL
}

# ──────────────────────────────────────────────────────────────────────────────
# 4. Standard Visualization Plots ####
# ──────────────────────────────────────────────────────────────────────────────

# 4.1 Plot Base Cropped Point Cloud
if (!is.null(las)) {
  plot(las, size = 4, bg = "white")
}

# 4.2 Plot Ground Classified Points & Digital Terrain Model (DTM)
if (!is.null(las_classified)) {
  gnd <- filter_ground(las_classified)
  plot(gnd, size = 4, bg = "white")
  
  # Generate a temporary 1m DTM using TIN for visualization
  dtm_tin_0 <- rasterize_terrain(las_classified, res = 0.05, algorithm = tin())
  plot_dtm3d(dtm_tin_0, bg = "white")
}

# 4.3 Plot Height-Normalised Point Cloud
if (!is.null(las_normalised)) {
  plot(las_normalised, color = "RGB", size = 3, bg = "white")
}

# 4.4 Plot 2D Canopy Height Model (CHM)
if (!is.null(las_chm)) {
  if (!is.null(PlotTrees) && nrow(PlotTrees) > 0) {
    
    # Plot Base Image using the full extent of the CHM
    plot(las_chm, col = height.colors(50), main = paste("Full CHM with Plot", Number, "Crown Overlay"))
    
    # Add Polygons (They will render accurately over the CHM extent)
    plot(st_geometry(PlotTrees), add = TRUE, border = "red", lwd = 2)
  } else {
    # Fallback if no polygons exist for this plot
    plot(las_chm, col = height.colors(50), main = paste("Full Site CHM (Plot", Number, "Polygons Not Found)"))
  }
}

# ──────────────────────────────────────────────────────────────────────────────
# 5. Individual Tree Detection (ITD) & Segmentation (ITS) ####
# ──────────────────────────────────────────────────────────────────────────────

MinimumTreeHeight <- 0.5

# Initialize the Local Maximum Filter (LMF) algorithm for treetops
lmf_algorithm <- lmf(ws = 3, hmin = MinimumTreeHeight, shape = "circular")

# Execute ITD and Segmentation if the required spatial data is available
if (!is.null(las_chm) && !is.null(las_normalised)) {
  
  # --- 1. Locate Individual Treetops (ITD) ---
  ttops <- locate_trees(las = las_chm, algorithm = lmf_algorithm)
  
  # --- 2. Segment the 3D Point Cloud (ITS) ---
  algo_dalponte <- dalponte2016(chm = las_chm, treetops = ttops)
  las_segmented <- segment_trees(las = las_normalised, algorithm = algo_dalponte)
  las_trees_only <- filter_poi(las_segmented, !is.na(treeID))
  
  # --- 3. Optional: Generate Dynamic Crown Polygons from Points ---
  delineated_crowns <- crown_metrics(las_trees_only, func = .stdtreemetrics, geom = "convex")
  
  # --- 4. Visualizations ---
  # Visualization A: 2D CHM with Treetops and Delineated Crowns
  plot(las_chm, col = height.colors(50), main = "CHM with Treetops & Crown Boundaries")
  plot(sf::st_geometry(ttops), add = TRUE, pch = 3, col = "black")
  if (!is.null(delineated_crowns)) {
    plot(sf::st_geometry(delineated_crowns), add = TRUE, border = "white", lwd = 2)
  }
  
  # Visualization B: 3D Segmented Point Cloud
  plot(las_trees_only,
       color = "treeID",
       colorPalette = pastel.colors(200),
       bg = "white",
       size = 3)
  
} else {
  message("Cannot perform ITD or Segmentation: Both CHM and Normalised point clouds are required.")
}

# ──────────────────────────────────────────────────────────────────────────────
# 6. Optional: 3D Animation Export ####
# ──────────────────────────────────────────────────────────────────────────────

# # --- Define spin motion ---
# spin <- spin3d(axis = c(0, 0, 1), rpm = 6)
# 
# # --- Save the animation ---
# movie3d(
#   movie = "Tree_Tops_animation",   # Base filename for output
#   dir = getwd(),                   # Output directory
#   spin,                            # The animation function defined above
#   duration = 10,                   # Animation length in seconds
#   fps = 25,                        # Frames per second
#   clean = TRUE,                    # Remove individual frame images after compiling
#   type = "gif"                     # Export format
# )