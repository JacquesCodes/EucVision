# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: END-TO-END BASELINE SITE-LEVEL DTM GENERATION 
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────
# Description: Processes raw baseline point clouds across the entire site using 
#              optimal geometric chunks for 32GB of RAM & 6 Workers, classifies ground, 
#              generates a Master DTM VRT, and crops the final raster to the 
#              study boundary.
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
library(lidR)
library(RCSF)
library(RMCC)
library(sf)
library(tictoc)
library(dplyr)
library(future)
library(terra)

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration & Path Management ####
# ──────────────────────────────────────────────────────────────────────────────
# Define your "Golden Baseline" date where ground is most visible
date_folder <- "30. 30 June 2026 (ALS)"

# Extract the date part and create a safe filename format
# (e.g., "17. 02 March 2026" -> "02_March_2026")
file_date <- sub("^\\d+\\.\\s*", "", date_folder)
file_date_safe <- gsub(" ", "_", file_date)

las_folder <- paste0("E:/Remote Sensing Media/", date_folder, "/03. Point Clouds/")

cropped_dir <- paste0("E:/Remote Sensing Media/", date_folder, "/03b. Point Clouds Cropped/")
denoised_dir <- paste0("E:/Remote Sensing Media/", date_folder, "/04b. Point Clouds Denoised/")
classified_dir <- paste0("E:/Remote Sensing Media/", date_folder, "/05. Point Clouds Ground Classified/")
dtm_dir <- paste0("E:/Remote Sensing Media/", date_folder, "/05b. Baseline Plot DTMs/")

# Ensure all output directories exist
for (dir in c(cropped_dir, denoised_dir, classified_dir, dtm_dir)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Enable parallel processing for lidR operations to run across 6 CPU logical processors.
# Note: Users should adjust 'workers' based on their available CPU logical processors and RAM.
plan(multisession, workers = 6)

# ──────────────────────────────────────────────────────────────────────────────
# 3. Spatial Data Loading & Boundary Definition ####
# ──────────────────────────────────────────────────────────────────────────────
# Load your final mask boundary to cookie-cut the site at the very end
boundary <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. QGIS Shapefiles/04. IMPACT Plot & Compartment Boundaries/IMPACT_Plot_&_Compartment_Boundaries_EPSG_2048.shp")

# ──────────────────────────────────────────────────────────────────────────────
# 3.5. Point Cloud Cropping (Fixing Sparse Edges) ####
# ──────────────────────────────────────────────────────────────────────────────

# Load the new OAL boundary to clip the raw data
oal_boundary <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. QGIS Shapefiles/05. IMPACT OAL Boundaries/IMPACT_Boundaries_EPSG_2048.shp")

# Buffer the boundary slightly (5m) to ensure we don't create artificial edge artifacts inside the study plots
oal_boundary_buffered <- st_buffer(oal_boundary, 5)

# Load the unclipped raw point clouds
ctg_raw_unclipped <- readLAScatalog(las_folder)

# Setup engine to write cropped tiles
opt_output_files(ctg_raw_unclipped) <- paste0(cropped_dir, "Tile_{XLEFT}_{YBOTTOM}_cropped")
opt_chunk_size(ctg_raw_unclipped) <- 100

tic()
print("Cropping raw point clouds to OAL Boundary to remove sparse edges...")
# clip_roi extracts only the points inside the shapefile and writes them to the cropped_dir
invisible(clip_roi(ctg_raw_unclipped, oal_boundary_buffered))
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 4. Noise Removal (Statistical Outlier Removal) ####
# ──────────────────────────────────────────────────────────────────────────────

# NEW: Load the dense, cropped point clouds instead of the raw ones
ctg_raw <- readLAScatalog(cropped_dir)

# --- CRITICAL SITE-LEVEL ENGINE SETTINGS ---
opt_independent_files(ctg_raw) <- FALSE   # Ignores original file overlaps
opt_chunk_size(ctg_raw) <- 100            # Breaks the site into 100x100m geometric tiles
opt_chunk_buffer(ctg_raw) <- 5            # 5m buffer is plenty for noise removal
opt_select(ctg_raw) <- "xyz"

# Name the outputs dynamically by their spatial coordinates
opt_output_files(ctg_raw) <- paste0(denoised_dir, "Tile_{XLEFT}_{YBOTTOM}_denoised")

tic()
print("Applying Statistical Outlier Removal (SOR) to the entire site...")

# Apply STRICT SOR algorithm to ruthlessly prune sub-surface pit noise
ctg_denoised <- classify_noise(ctg_raw, sor(k = 25, m = 1.2))
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 5. Ground Classification (Progressive TIN Densification) ####
# ──────────────────────────────────────────────────────────────────────────────

# Reload the explicitly denoised 100x100m chunks from disk
ctg_denoised <- readLAScatalog(denoised_dir)

# Tell the engine to drop the noise points (Class 7) in the next step
opt_filter(ctg_denoised) <- "-drop_class 7"

# --- CRITICAL SITE-LEVEL ENGINE SETTINGS ---
opt_independent_files(ctg_denoised) <- FALSE
opt_chunk_size(ctg_denoised) <- 100
opt_chunk_buffer(ctg_denoised) <- 10      # Massive 10m buffer allows PTD to reach the roads!
opt_select(ctg_denoised) <- "xyz"
opt_output_files(ctg_denoised) <- paste0(classified_dir, "Tile_{XLEFT}_{YBOTTOM}_classified")

tic()
print("Applying Progressive TIN Densification (PTD) tuned for SfM Canopy Shroud...")

# Tuned specifically for dense SfM canopy crusts using Progressive TIN Densification (PTD)
ctg_classified <- classify_ground(
  ctg_denoised, 
  ptd(res = 12,          # 12m guarantees the search grid overhangs canopy widths
      angle = 6,         # Strict 6-degree limit prevents the TIN from climbing canopy walls
      distance = 1.5)    # Keeps the vertical step distance tight
)
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 6. DTM Generation (Site Level) ####
# ──────────────────────────────────────────────────────────────────────────────

# Reload the explicitly classified chunks from disk
ctg_classified <- readLAScatalog(classified_dir)

# --- CRITICAL SITE-LEVEL ENGINE SETTINGS ---
opt_independent_files(ctg_classified) <- FALSE
opt_chunk_size(ctg_classified) <- 100          
opt_chunk_buffer(ctg_classified) <- 2          # TIN only needs a 2m buffer to stitch seams

opt_select(ctg_classified) <- "xyzc"
opt_filter(ctg_classified) <- "-keep_class 2"  # Only use ground points

opt_output_files(ctg_classified) <- paste0(dtm_dir, "Tile_{XLEFT}_{YBOTTOM}_DTM")

tic()
print("Generating continuous site DTMs at 0.05m resolution...")
site_dtms <- rasterize_terrain(ctg_classified, 
                               res = 0.05, 
                               algorithm = tin())
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 7. Master DTM Consolidation (VRT & Single File) ####
# ──────────────────────────────────────────────────────────────────────────────
tic()
print("Building Master Baseline VRT and consolidating single physical DTM...")

dtm_files <- list.files(dtm_dir, pattern = "\\.tif$", full.names = TRUE)
vrt_path <- paste0(dtm_dir, "Master_Baseline_DTM_", file_date_safe, ".vrt")
single_dtm_path <- paste0(dtm_dir, "Master_Baseline_DTM_Single_", file_date_safe, ".tif")

# Create the VRT and write it out to a single physical .tif file
site_dtm_vrt <- terra::vrt(dtm_files, vrt_path, overwrite = TRUE)
terra::writeRaster(site_dtm_vrt, filename = single_dtm_path, overwrite = TRUE)

print(paste("- Master VRT successfully saved to:", vrt_path))
print(paste("- Raw Physical DTM successfully saved to:", single_dtm_path))
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 8. Spatial DTM Sinkhole Filling & Final Boundary Masking ####
# ──────────────────────────────────────────────────────────────────────────────
tic()
print("Applying Morphological Closing to destroy SfM sinkholes...")

smoothed_dtm_path <- paste0(dtm_dir, "Master_Baseline_DTM_Smoothed_", file_date_safe, ".tif")
raw_dtm <- terra::rast(single_dtm_path)

# 1. Morphological Closing for <= 2m holes
print("Step 1/5: Dilating terrain over small sinkholes...")
dtm_dilated <- terra::focal(raw_dtm, w = 41, fun = "max", na.rm = TRUE)

print("Step 2/5: Eroding terrain back to true natural slopes...")
dtm_closed <- terra::focal(dtm_dilated, w = 41, fun = "min", na.rm = TRUE)

# 2. Targeted Surgical Patch for Massive (4-5m) Sinkholes
print("Step 3/5: Targeted surgical fill for massive 5m sinkholes...")
reference_surface <- terra::focal(dtm_closed, w = 101, fun = "mean", na.rm = TRUE)

dtm_patched <- dtm_closed
dtm_patched[dtm_patched < (reference_surface - 0.5)] <- NA
dtm_filled <- terra::focal(dtm_patched, w = 111, fun = "mean", na.rm = TRUE, na.only = TRUE)

# 3. Final Light Blend
print("Step 4/5: Final surface blending...")
final_dtm <- terra::focal(dtm_filled, w = 5, fun = "median", na.rm = TRUE)

# 4. Boundary Masking (Shaving the Edge Artifacts and Excess Site Data)
print("Step 5/5: Cropping continuous DTM to the specific study boundary...")

boundary_vect <- terra::vect(boundary)
if (terra::crs(boundary_vect) != terra::crs(final_dtm)) {
  boundary_vect <- terra::project(boundary_vect, terra::crs(final_dtm))
}

# Crop reduces the processing extent, Mask cookie-cuts the exact shape
final_dtm_clipped <- terra::crop(final_dtm, boundary_vect)
final_dtm_clipped <- terra::mask(final_dtm_clipped, boundary_vect)

# Export the finalized, absolutely perfect DTM
terra::writeRaster(final_dtm_clipped, filename = smoothed_dtm_path, overwrite = TRUE)

print(paste("- SINKHOLE-FREE & CLIPPED DTM successfully saved to:", smoothed_dtm_path))
print("Pipeline Complete!")
toc()