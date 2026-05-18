# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: PLOT-LEVEL LiDAR POINT CLOUD PROCESSING PIPELINE USING DTM ####
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
# Load required libraries for point cloud processing, spatial operations, and parallel computing
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
date_folder <- "20. 23 March 2026"

# Extract the date part and create a safe filename format
file_date <- sub("^\\d+\\.\\s*", "", date_folder)
file_date_safe <- gsub(" ", "_", file_date)

# Define processing directories
las_folder     <- paste0("E:/Remote Sensing Media/", date_folder, "/03. Point Clouds/")
clipped_dir    <- paste0("E:/Remote Sensing Media/", date_folder, "/04. Point Clouds Clipped/")
normalised_dir <- paste0("E:/Remote Sensing Media/", date_folder, "/06. Point Clouds Normalised/")
chm_dir        <- paste0("E:/Remote Sensing Media/", date_folder, "/07. Canopy Height Models/")
polygons_dir   <- paste0("E:/Remote Sensing Media/", date_folder, "/08. Crown Polygons/")
metrics_dir    <- paste0("E:/Remote Sensing Media/", date_folder, "/09. Crown Metrics/")

# Define the absolute path to your Baseline DTM (Stays constant across batches)
baseline_dtm_path <- "E:/Remote Sensing Media/03. 30 October 2025/05b. Baseline Plot DTMs/Master_Baseline_DTM_Smoothed_30_October_2025.tif"

# Ensure all output directories exist
for (dir in c(clipped_dir, normalised_dir, chm_dir, metrics_dir)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# ──────────────────────────────────────────────────────────────────────────────
# 3. Spatial Data Loading & Boundary Definition ####
# ──────────────────────────────────────────────────────────────────────────────
# Load individual plot boundary shapefiles and sort them chronologically/spatially by ID
plots_buffered_unsorted <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LAScatalog Boundaries/LAScatalog Plot Boundaries.shp")
plots <- plots_buffered_unsorted[order(plots_buffered_unsorted$id), ]

# --- BATCH PROCESSING & CATALOG SETUP ---
all_las_files <- list.files(las_folder, pattern = "\\.(las|laz)$", full.names = TRUE, ignore.case = TRUE)

# Identify Top and Bottom point cloud files based on filenames
top_files <- all_las_files[grepl("Top", basename(all_las_files), ignore.case = TRUE)]
bot_files <- all_las_files[grepl("Bottom", basename(all_las_files), ignore.case = TRUE)]

# Load individual tree crown polygons for final height extraction
trees <- st_read(paste0(polygons_dir, "Crown_Polygons_", file_date_safe, ".shp"))

# Validate and enforce the correct Coordinate Reference System (EPSG: 2048)
if (is.na(st_crs(trees)$epsg) || st_crs(trees)$epsg != 2048) {
  trees <- st_transform(trees, 2048)
  print("Transformed CRS to 2048 successfully.")
}

# Clean tree geometry dataset by removing unnecessary metadata columns
trees <- trees %>%
  select(-any_of(c("group_ulid", "N_GM", "id", "N_FG", "N_BG", "BBox")))

# ──────────────────────────────────────────────────────────────────────────────
# 4. Point Cloud Cropping (Plot Level) ####
# ──────────────────────────────────────────────────────────────────────────────
plan(multisession) # Enable parallel processing
tic()              # Track execution time

# Check if we need to split processing based on Top/Bottom files
if (length(top_files) > 0 && length(bot_files) > 0) {
  print("Top and Bottom point clouds detected. Splitting processing by Plot ID...")
  
  plots_top <- plots %>% filter(id <= 21)
  plots_bot <- plots %>% filter(id >= 22)
  
  ctg_top <- readLAScatalog(top_files)
  ctg_bot <- readLAScatalog(bot_files)
  
  # Configure engine settings for TOP cropping
  opt_independent_files(ctg_top) <- FALSE
  opt_select(ctg_top) <- "xyz"
  opt_output_files(ctg_top) <- paste0(clipped_dir, "Plot_{id}_", file_date_safe)
  
  # Configure engine settings for BOTTOM cropping
  opt_independent_files(ctg_bot) <- FALSE
  opt_select(ctg_bot) <- "xyz"
  opt_output_files(ctg_bot) <- paste0(clipped_dir, "Plot_{id}_", file_date_safe)
  
  print("Cropping Top plots (1-21)...")
  ctg_clipped_top <- clip_roi(ctg_top, plots_top)
  
  print("Cropping Bottom plots (22-75)...")
  ctg_clipped_bot <- clip_roi(ctg_bot, plots_bot)
  
  ctg_clipped <- readLAScatalog(clipped_dir)
  
} else {
  print("Single point cloud or no Top/Bottom distinction detected. Processing entirely...")
  
  ctg <- readLAScatalog(las_folder)
  opt_independent_files(ctg) <- FALSE
  opt_select(ctg) <- "xyz"
  opt_output_files(ctg) <- paste0(clipped_dir, "Plot_{id}_", file_date_safe)
  
  ctg_clipped <- clip_roi(ctg, plots)
}
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 5. Ground Classification (Not needed) ####
# ──────────────────────────────────────────────────────────────────────────────

# Not needed when using baseline DTM for Normalization

# ──────────────────────────────────────────────────────────────────────────────
# 6. Height Normalization (USING BASELINE DTM) ####
# ──────────────────────────────────────────────────────────────────────────────
ctg_clipped <- readLAScatalog(clipped_dir)

plan(multisession, workers = 6)

opt_independent_files(ctg_clipped) <- TRUE
opt_select(ctg_clipped) <- "xyz"
opt_output_files(ctg_clipped) <- paste0(normalised_dir, "{ORIGINALFILENAME}_classified_normalised")

# Load the Master Baseline Smoothed DTM generated in October 2025
baseline_dtm <- rast(baseline_dtm_path)

tic()
print("Normalizing point clouds against the October baseline DTM...")

# Perform raster-based elevation normalization
ctg_normalised <- normalize_height(las = ctg_clipped, algorithm = baseline_dtm)
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 7. Canopy Height Model (CHM) Generation ####
# ──────────────────────────────────────────────────────────────────────────────
# RAM Optimization Strategies (Useful for heavy workloads):
# 1. Use smaller chunk/plot sizes.
# 2. Use opt_select() to load only needed fields into memory.
# 3. Decrease active workers (threads) since each worker consumes its own RAM footprint.
# 4. Exclude ground points and sub-surface noise.

ctg_normalised <- readLAScatalog(normalised_dir)

# Limit the amount of workers (threads) if RAM constrained
plan(multisession, workers = 6)
opt_independent_files(ctg_normalised) <- TRUE
opt_select(ctg_normalised) <- "xyz"

# Filter noise: Drop ground points and sub-surface/extreme-height noise
opt_filter(ctg_normalised) <- "-drop_z_below 0 -drop_z_above 30"

tic()
opt_output_files(ctg_normalised) <- paste0(chm_dir, "{*}_chm")

# Rasterize highest points into a continuous CHM grid with interpolation
ctg_chm <- rasterize_canopy(ctg_normalised,
                            res = 0.05,
                            algorithm = p2r(na.fill = tin()))
print("Rasterize canopy time:")
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 8. Master CHM Consolidation & Dynamic Metric Extraction ####
# ──────────────────────────────────────────────────────────────────────────────
tic()
print("Extracting metrics and consolidating dynamic master CHM...")

# 1. Gather all individual CHM files and create a Virtual Raster (VRT)
chm_files <- list.files(chm_dir, pattern = "\\.tif$", full.names = TRUE)
site_chm_vrt <- terra::vrt(chm_files)

print("Step 1/2: Extracting max heights from VRT to calculate dynamic cap...")

# 2. Extract exact maximum tree heights directly from the VRT in memory
trees$Tree_Height <- exact_extract(site_chm_vrt, trees, 'max')

# 3. Find the absolute tallest valid tree to use as our dynamic ceiling.
# We use is.finite() to ignore any Inf edge-artifacts that might have sneaked 
# into polygons sitting directly on the plot boundary.
max_tree_height <- max(trees$Tree_Height[is.finite(trees$Tree_Height)], na.rm = TRUE)

# Round up to the nearest whole meter for a clean QGIS color ramp (e.g., 22.4m -> 23m)
dynamic_cap <- ceiling(max_tree_height)
print(paste("-> Dynamic CHM cap safely set to:", dynamic_cap, "meters"))

# 4. CLAMP THE ARTIFACTS dynamically
site_chm_clamped <- terra::clamp(site_chm_vrt, lower = 0, upper = dynamic_cap)

print("Step 2/2: Exporting perfectly capped Master CHM and metrics...")

# 5. Write the perfectly capped raster out to a single physical .tif file
single_chm_path <- paste0(chm_dir, "Master_Site_CHM_Single_", file_date_safe, ".tif")
terra::writeRaster(site_chm_clamped, filename = single_chm_path, overwrite = TRUE)

# 6. Save to shapefile (completely overwriting old files to prevent schema errors)
st_write(trees, paste0(metrics_dir, "Crown_Metrics_", file_date_safe, ".shp"), delete_dsn = TRUE)

# 7. Save lightweight tabular CSV data (drops the messy spatial geometry text)
write.csv(st_drop_geometry(trees), paste0(metrics_dir, "Crown_Metrics_", file_date_safe, ".csv"), row.names = FALSE)

print(paste("Master CHM successfully saved to:", single_chm_path))
toc()