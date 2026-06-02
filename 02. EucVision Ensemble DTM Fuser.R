# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: ULTIMATE ENSEMBLE BASELINE DTM FUSION 
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────
# Description: Fuses multiple temporal baseline Digital Terrain Models (DTMs) 
#              into a single, robust "Ultimate" DTM. It uses a pixel-wise 
#              temporal maximum approach to eliminate transient SfM sinkholes, 
#              ensuring the most reliable solid ground model is preserved.
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
library(terra)   # Core package for spatial raster operations
library(tictoc)  # For tracking script execution time

tic()
print("Fusing multiple temporal DTMs into an Ultimate Baseline DTM...")

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration & Path Management ####
# ──────────────────────────────────────────────────────────────────────────────
# Define the paths to your individual, smoothed baseline DTMs
# These represent different flight dates where ground visibility may have varied
# path_dtm_1 <- "E:/Remote Sensing Media/03. 30 October 2025/05b. Baseline Plot DTMs/Master_Baseline_DTM_Smoothed_30_October_2025.tif"
# path_dtm_2 <- "E:/Remote Sensing Media/04. 07 November 2025/05b. Baseline Plot DTMs/Master_Baseline_DTM_Smoothed_07_November_2025.tif"
# path_dtm_3 <- "E:/Remote Sensing Media/05. 14 November 2025/05b. Baseline Plot DTMs/Master_Baseline_DTM_Smoothed_14_November_2025.tif"
path_dtm_1 <- "E:/Remote Sensing Media/03. 30 October 2025/05b. Baseline Plot DTMs/Master_Baseline_DTM_Single_30_October_2025.tif"
path_dtm_2 <- "E:/Remote Sensing Media/04. 07 November 2025/05b. Baseline Plot DTMs/Master_Baseline_DTM_Single_07_November_2025.tif"
path_dtm_3 <- "E:/Remote Sensing Media/05. 14 November 2025/05b. Baseline Plot DTMs/Master_Baseline_DTM_Single_14_November_2025.tif"

# Define the final output destination
# output_path <- "E:/Remote Sensing Media/00. Baseline DTM/Ultimate_Ensemble_Baseline_DTM.tif"
output_path <- "E:/Remote Sensing Media/00. Baseline DTM/IMPACT_OAL_Baseline_DTM.tif"

# Ensure the output directory exists before attempting to write out
if (!dir.exists(dirname(output_path))) dir.create(dirname(output_path), recursive = TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# 3. Spatial Data Loading & Alignment ####
# ──────────────────────────────────────────────────────────────────────────────
print("Loading individual DTMs and enforcing spatial alignment...")

# Read the raster files into memory as SpatRasters
dtm_1 <- terra::rast(path_dtm_1)
dtm_2 <- terra::rast(path_dtm_2)
dtm_3 <- terra::rast(path_dtm_3)

# Ensure they all perfectly align geometrically
# DTM 1 acts as the master reference grid. If extents or resolutions differ, 
# we force them to match using bilinear interpolation.
if (!ext(dtm_2) == ext(dtm_1)) dtm_2 <- terra::resample(dtm_2, dtm_1, method = "bilinear")
if (!ext(dtm_3) == ext(dtm_1)) dtm_3 <- terra::resample(dtm_3, dtm_1, method = "bilinear")

# ──────────────────────────────────────────────────────────────────────────────
# 4. Temporal Fusion ####
# ──────────────────────────────────────────────────────────────────────────────
print("Applying pixel-wise temporal maximum across all dates...")

# Stack them into a single multi-layer SpatRaster
dtm_stack <- c(dtm_1, dtm_2, dtm_3)

# Apply the Pixel-wise Temporal MAXIMUM
# Because PTD SfM errors are inherently biased downward (sinkholes), taking the 
# absolute highest elevation value across all 3 dates ensures that if even ONE 
# flight caught the true solid ground, it overwrites the sinkholes from the others.
ultimate_dtm <- terra::app(dtm_stack, fun = "max", na.rm = TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# 5. Export Final DTM ####
# ──────────────────────────────────────────────────────────────────────────────
print("Writing out final surface model to disk...")

# Export the fused Baseline DTM
terra::writeRaster(ultimate_dtm, filename = output_path, overwrite = TRUE)

print(paste("Ultimate Ensemble DTM saved to:", output_path))
toc()