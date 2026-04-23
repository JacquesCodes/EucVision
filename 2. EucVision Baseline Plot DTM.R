# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: END-TO-END BASELINE PLOT-LEVEL DTM GENERATION 
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/) 
# Description: Processes raw baseline point clouds (crops to plots, classifies 
#              ground) and generates a Master DTM VRT for future normalization.
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
date_folder <- "03. 30 October 2025"
file_date_safe <- gsub(" ", "_", sub("^\\d+\\.\\s*", "", date_folder))

las_folder <- paste0("E:/Remote Sensing Media/", date_folder, "/03. Point Clouds")
clipped_dir <- paste0("E:/Remote Sensing Media/", date_folder, "/04. Point Clouds Clipped/")
denoised_dir <- paste0("E:/Remote Sensing Media/", date_folder, "/04b. Point Clouds Denoised/")
classified_dir <- paste0("E:/Remote Sensing Media/", date_folder, "/05. Point Clouds Ground Classified/")
dtm_dir <- paste0("E:/Remote Sensing Media/", date_folder, "/05b. Baseline Plot DTMs/")

# Ensure output directories exist
for (dir in c(clipped_dir, denoised_dir, classified_dir, dtm_dir)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# ──────────────────────────────────────────────────────────────────────────────
# 3. Spatial Data Loading & Boundary Definition ####
# ──────────────────────────────────────────────────────────────────────────────
# Load individual plot boundary shapefiles and sort them chronologically by ID
plots_buffered_unsorted <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LAScatalog Boundaries/LAScatalog Plot Boundaries.shp")
plots <- plots_buffered_unsorted[order(plots_buffered_unsorted$id), ]

# Load your new mask boundary
boundary <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LAScatalog Boundaries/IMPACT Plot & Compartment Boundaries 2048.shp")

all_las_files <- list.files(las_folder, pattern = "\\.(las|laz)$", full.names = TRUE, ignore.case = TRUE)
top_files <- all_las_files[grepl("Top", basename(all_las_files), ignore.case = TRUE)]
bot_files <- all_las_files[grepl("Bottom", basename(all_las_files), ignore.case = TRUE)]

# ──────────────────────────────────────────────────────────────────────────────
# 4. Point Cloud Cropping (Plot Level) ####
# ──────────────────────────────────────────────────────────────────────────────
plan(multisession) # Enable parallel processing
tic()              

if (length(top_files) > 0 && length(bot_files) > 0) {
  print("Top and Bottom point clouds detected. Splitting processing by Plot ID...")
  
  plots_top <- plots %>% filter(id <= 21)
  plots_bot <- plots %>% filter(id >= 22)
  
  ctg_top <- readLAScatalog(top_files)
  ctg_bot <- readLAScatalog(bot_files)
  
  # Configure TOP cropping
  opt_independent_files(ctg_top) <- FALSE
  opt_select(ctg_top) <- "xyz"
  opt_output_files(ctg_top) <- paste0(clipped_dir, "Plot_{id}_", file_date_safe)
  
  # Configure BOTTOM cropping
  opt_independent_files(ctg_bot) <- FALSE
  opt_select(ctg_bot) <- "xyz"
  opt_output_files(ctg_bot) <- paste0(clipped_dir, "Plot_{id}_", file_date_safe)
  
  print("Cropping Top plots (1-21)...")
  ctg_clipped_top <- clip_roi(ctg_top, plots_top)
  
  print("Cropping Bottom plots (22-75)...")
  ctg_clipped_bot <- clip_roi(ctg_bot, plots_bot)
  
} else {
  print("Single point cloud detected. Processing entirely...")
  ctg <- readLAScatalog(las_folder)
  opt_independent_files(ctg) <- FALSE
  opt_select(ctg) <- "xyz"
  opt_output_files(ctg) <- paste0(clipped_dir, "Plot_{id}_", file_date_safe)
  ctg_clipped <- clip_roi(ctg, plots)
}
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 5. Noise Removal (Statistical Outlier Removal) ####
# ──────────────────────────────────────────────────────────────────────────────
plan(multisession)

# Reload the explicitly cropped chunks from disk
ctg_clipped <- readLAScatalog(clipped_dir)

opt_independent_files(ctg_clipped) <- FALSE   # Suppresses overlap warning
opt_chunk_buffer(ctg_clipped) <- 2            # Prevents edge artifacts during noise detection
opt_output_files(ctg_clipped) <- paste0(denoised_dir, "{ORIGINALFILENAME}_denoised")

tic()
print("Applying Statistical Outlier Removal (SOR)...")

# 2. Apply a STRICTER SOR algorithm to ruthlessly prune sub-surface pit noise
# k = 25 (forces the algorithm to look at a much wider neighborhood of points)
# m = 1.2 (drops points that are just 1.2 standard deviations away from the average, down from 2.0)
ctg_denoised <- classify_noise(ctg_clipped, sor(k = 25, m = 1.2))
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 6. Ground Classification (Progressive TIN Densification) ####
# ──────────────────────────────────────────────────────────────────────────────
plan(multisession)

# Reload the explicitly denoised chunks from disk
ctg_denoised <- readLAScatalog(denoised_dir)

# Tell the engine to drop the noise points (Class 7) in the next step
opt_filter(ctg_denoised) <- "-drop_class 7"

# Setup engine for Ground Classification
opt_independent_files(ctg_denoised) <- FALSE
opt_chunk_size(ctg_denoised) <- 0
opt_chunk_buffer(ctg_denoised) <- 2
opt_select(ctg_denoised) <- "xyz"
opt_output_files(ctg_denoised) <- paste0(classified_dir, "{ORIGINALFILENAME}_classified")

tic()
print("Applying Progressive TIN Densification (PTD) tuned for SfM Canopy Shroud...")

# Tuned specifically for dense SfM canopy crusts
ctg_classified <- classify_ground(
  ctg_denoised, 
  ptd(res = 12,          # 12m guarantees the search grid overhangs canopy widths
      angle = 6,         # Strict 6-degree limit prevents the TIN from climbing canopy walls
      distance = 1.5)    # Keeps the vertical step distance tight
)
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 7. DTM Generation (Plot Level) ####
# ──────────────────────────────────────────────────────────────────────────────
plan(multisession, workers = 6) 

# Reload the explicitly classified chunks from disk
ctg_classified <- readLAScatalog(classified_dir)

# --- CRITICAL ENGINE SETTINGS FOR SEAMLESS 1:1 PLOT EXPORT ---
opt_independent_files(ctg_classified) <- FALSE # Allow access to neighbors
opt_chunk_size(ctg_classified) <- 0            # Keep 1-to-1 plot output
opt_chunk_buffer(ctg_classified) <- 2          # TIN needs buffer to anchor triangles seamlessly

opt_select(ctg_classified) <- "xyzc"
opt_filter(ctg_classified) <- "-keep_class 2"  # Only use ground points

# Stream .tif outputs directly to disk
opt_output_files(ctg_classified) <- paste0(dtm_dir, "{ORIGINALFILENAME}_DTM")

tic()
print("Generating baseline DTMs at 0.05m resolution...")
plot_dtms <- rasterize_terrain(ctg_classified, 
                               res = 0.05, 
                               algorithm = tin())
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 8. Master DTM Consolidation (VRT & Single File) ####
# ──────────────────────────────────────────────────────────────────────────────
tic()
print("Building Master Baseline VRT and consolidating single physical DTM...")

dtm_files <- list.files(dtm_dir, pattern = "\\.tif$", full.names = TRUE)
vrt_path <- paste0(dtm_dir, "Master_Baseline_DTM_", file_date_safe, ".vrt")
single_dtm_path <- paste0(dtm_dir, "Master_Baseline_DTM_Single_", file_date_safe, ".tif")

# 1. Create the VRT using terra
site_dtm_vrt <- terra::vrt(dtm_files, vrt_path, overwrite = TRUE)

# 2. Write the VRT out to a single physical .tif file
terra::writeRaster(site_dtm_vrt, filename = single_dtm_path, overwrite = TRUE)

print(paste("- Master VRT successfully saved to:", vrt_path))
print(paste("- Raw Physical DTM successfully saved to:", single_dtm_path))
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 9. Spatial DTM Sinkhole Filling (Morphological & Targeted Deep Fill) ####
# ──────────────────────────────────────────────────────────────────────────────
tic()
print("Applying Morphological Closing to destroy small SfM sinkholes...")

smoothed_dtm_path <- paste0(dtm_dir, "Master_Baseline_DTM_Smoothed_", file_date_safe, ".tif")
raw_dtm <- terra::rast(single_dtm_path)

# 1. Morphological Closing for <= 2m holes
print("Step 1/5: Dilating terrain over small sinkholes...")
dtm_dilated <- terra::focal(raw_dtm, w = 41, fun = "max", na.rm = TRUE)

print("Step 2/5: Eroding terrain back to true natural slopes...")
dtm_closed <- terra::focal(dtm_dilated, w = 41, fun = "min", na.rm = TRUE)

# 2. Targeted Surgical Patch for Massive (4-5m) Sinkholes
print("Step 3/5: Targeted surgical fill for massive 5m sinkholes...")

# Create a heavily smoothed 5m reference surface to identify the massive drops
# Swapped to "mean" to prevent the laptop from freezing for hours
reference_surface <- terra::focal(dtm_closed, w = 101, fun = "mean", na.rm = TRUE)

# If the ground suddenly drops more than 0.5m below the general local area, 
# it's a massive SfM artifact. Turn those pixels into NA (a blank void).
dtm_patched <- dtm_closed
dtm_patched[dtm_patched < (reference_surface - 0.5)] <- NA

# Fill ONLY the NA voids using a large moving average of the healthy edges.
# w = 111 covers 5.55m, allowing it to easily bridge across the massive holes.
# na.only = TRUE ensures that the rest of your beautiful terrain is untouched!
dtm_filled <- terra::focal(dtm_patched, w = 111, fun = "mean", na.rm = TRUE, na.only = TRUE)

# 3. Final Light Blend
print("Step 4/5: Final surface blending...")
final_dtm <- terra::focal(dtm_filled, w = 5, fun = "median", na.rm = TRUE)

# 4. Boundary Masking (Shaving the Edge Artifacts)
print("Step 5/5: Shaving boundary edge artifacts...")

# Convert the sf boundary object to terra's vector format
boundary_vect <- terra::vect(boundary)

# Ensure Coordinate Reference Systems match perfectly
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