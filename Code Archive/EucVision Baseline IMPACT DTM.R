# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: BASELINE STAND-LEVEL DTM GENERATION PIPELINE 
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# Description: Processes juvenile Eucalyptus SfM point clouds to classify ground
#              points and generate a 0.05m high-resolution Digital Terrain Model.
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Environment
# ──────────────────────────────────────────────────────────────────────────────
library(lidR)
library(terra)
library(future)

# Enable parallel processing for heavy point cloud operations.
# Constrained to 4 workers to safely manage RAM limits during processing.
plan(multisession, workers = 4)

# ──────────────────────────────────────────────────────────────────────────────
# 2. Path Configuration
# ──────────────────────────────────────────────────────────────────────────────
input_las_dir  <- "E:/Remote Sensing Media/03. 30 October 2025/03. Point Clouds"
chunk_out_dir  <- "E:/Remote Sensing Media/00. Baseline DTM and plot cropping/Classified_Chunks"
output_dtm     <- "E:/Remote Sensing Media/00. Baseline DTM and plot cropping/Baseline_DTM_0.05m.tif"

# Load the baseline point cloud catalog
baseline_ctg <- readLAScatalog(input_las_dir)

# ──────────────────────────────────────────────────────────────────────────────
# 3. Engine & Memory Optimization
# ──────────────────────────────────────────────────────────────────────────────
# Configure chunk sizes to balance processing speed and RAM limitations
opt_chunk_size(baseline_ctg)   <- 50 
opt_chunk_buffer(baseline_ctg) <- 5

# Stream chunks directly to disk to prevent memory overflow
opt_output_files(baseline_ctg) <- paste0(chunk_out_dir, "/Tile_{XLEFT}_{YBOTTOM}_ground")

# SPEED OPTIMIZATION: Load only spatial coordinates (X, Y, Z) into RAM.
# Automatically drops heavy RGB color data which is unnecessary for classification.
opt_select(baseline_ctg) <- "xyz"

# ──────────────────────────────────────────────────────────────────────────────
# 4. Ground Classification (Cloth Simulation Filter)
# ──────────────────────────────────────────────────────────────────────────────
# Apply CSF to invert the point cloud and drape a simulated cloth to find the ground.
baseline_classified_ctg <- classify_ground(baseline_ctg, csf(
  sloop_smooth     = TRUE,  # Adjusts cloth to adhere naturally to sloped terrain
  class_threshold  = 0.05,  # 5cm threshold captures the vertical "fuzziness" of SfM ground points
  cloth_resolution = 1,     # 1m node spacing suited for 1x1m plots
  rigidness        = 3,     # CRITICAL: Max stiffness prevents the cloth from sagging into juvenile canopy gaps
  time_step        = 1      # Standard gravity simulation speed
))

# ──────────────────────────────────────────────────────────────────────────────
# 5. Rasterize Master Digital Terrain Model (DTM)
# ──────────────────────────────────────────────────────────────────────────────
# Reload the newly classified chunks as a fresh catalog for rasterization
baseline_classified_ctg <- readLAScatalog(chunk_out_dir)

# Generate a continuous surface model using the Triangulated Irregular Network (TIN) algorithm.
# Resolution kept at 0.05m to capture fine micro-topography beneath small trees.
master_dtm <- rasterize_terrain(baseline_classified_ctg,
                                res = 0.05, 
                                algorithm = tin())

# ──────────────────────────────────────────────────────────────────────────────
# 6. Export Final Raster
# ──────────────────────────────────────────────────────────────────────────────
# Save the DTM to disk (overwriting any previous iterations)
writeRaster(master_dtm, output_dtm, overwrite = TRUE)

print(paste("Processing complete. Master DTM saved to:", output_dtm))