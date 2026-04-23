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
date_folder <- "07. 28 November 2025 DTM"

# Extract the date part and create a safe filename format
# (e.g., "17. 02 March 2026" -> "02_March_2026")
file_date <- sub("^\\d+\\.\\s*", "", date_folder)
file_date_safe <- gsub(" ", "_", file_date)

# ──────────────────────────────────────────────────────────────────────────────
# 3. Spatial Data Loading & Boundary Definition ####
# ──────────────────────────────────────────────────────────────────────────────
# Load individual plot boundary shapefiles and sort them chronologically/spatially by ID
plots_buffered_unsorted <- st_read(paste0("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LAScatalog Boundaries/LAScatalog Plot Boundaries.shp"))
plots <- plots_buffered_unsorted[order(plots_buffered_unsorted$id), ]

# --- BATCH PROCESSING & CATALOG SETUP ---
las_folder <- paste0("E:/Remote Sensing Media/", date_folder, "/03. Point Clouds")
all_las_files <- list.files(las_folder, pattern = "\\.(las|laz)$", full.names = TRUE, ignore.case = TRUE)

# Identify Top and Bottom point cloud files based on filenames
top_files <- all_las_files[grepl("Top", basename(all_las_files), ignore.case = TRUE)]
bot_files <- all_las_files[grepl("Bottom", basename(all_las_files), ignore.case = TRUE)]

# Load individual tree crown polygons for final height extraction
trees <- st_read(paste0("E:/Remote Sensing Media/", date_folder, "/08. Crown Polygons/Crown_Polygons_", file_date_safe, ".shp"))

# Validate and enforce the correct Coordinate Reference System (EPSG: 2048)
if (is.na(st_crs(trees)$epsg) || st_crs(trees)$epsg != 2048) {
  trees <- st_transform(trees, 2048)
  print("Transformed CRS to 2048 successfully.")
}

# Clean tree geometry dataset by removing unnecessary metadata columns
trees <- trees %>%
  select(-any_of(c("group_ulid", "N_GM", "id", "N_FG", "N_BG", "BBox")))

# Optional Debugging Visualization (Uncomment if needed)
# plot(ctg)
# plot(plots$geometry, add = TRUE)
# plot(trees$geometry, add = TRUE, col = "red")

# ──────────────────────────────────────────────────────────────────────────────
# 4. Point Cloud Cropping (Plot Level) ####
# ──────────────────────────────────────────────────────────────────────────────
plan(multisession) # Enable parallel processing
tic()              # Track execution time

# Check if we need to split processing based on Top/Bottom files
if (length(top_files) > 0 && length(bot_files) > 0) {
  print("Top and Bottom point clouds detected. Splitting processing by Plot ID...")
  
  # Split plot boundaries by ID to align with the physical flight paths
  plots_top <- plots %>% filter(id <= 21)
  plots_bot <- plots %>% filter(id >= 22)
  
  # Initialize separate LiDAR catalogs
  ctg_top <- readLAScatalog(top_files)
  ctg_bot <- readLAScatalog(bot_files)
  
  # Configure engine settings for TOP cropping
  opt_independent_files(ctg_top) <- FALSE
  opt_select(ctg_top) <- "xyz"
  # Note: using {id} instead of {ID} to match the dataframe's column name casing
  opt_output_files(ctg_top) <- paste0("E:/Remote Sensing Media/", date_folder, "/04. Point Clouds Clipped/", "Plot_{id}_", file_date_safe)
  
  # Configure engine settings for BOTTOM cropping
  opt_independent_files(ctg_bot) <- FALSE
  opt_select(ctg_bot) <- "xyz"
  opt_output_files(ctg_bot) <- paste0("E:/Remote Sensing Media/", date_folder, "/04. Point Clouds Clipped/", "Plot_{id}_", file_date_safe)
  
  # Execute spatial clipping operations
  print("Cropping Top plots (1-21)...")
  ctg_clipped_top <- clip_roi(ctg_top, plots_top)
  
  print("Cropping Bottom plots (22-75)...")
  ctg_clipped_bot <- clip_roi(ctg_bot, plots_bot)
  
  # Re-read the fully clipped directory as a single unified catalog for downstream processing
  ctg_clipped <- readLAScatalog(paste0("E:/Remote Sensing Media/", date_folder, "/04. Point Clouds Clipped"))
  
} else {
  print("Single point cloud or no Top/Bottom distinction detected. Processing entirely...")
  
  # Standard processing for a single or unstructured catalog
  ctg <- readLAScatalog(las_folder)
  opt_independent_files(ctg) <- FALSE
  opt_select(ctg) <- "xyz"
  opt_output_files(ctg) <- paste0("E:/Remote Sensing Media/", date_folder, "/04. Point Clouds Clipped/", "Plot_{id}_", file_date_safe)
  
  # Execute spatial clipping operations
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
plan(multisession, workers = 6)

opt_independent_files(ctg_clipped) <- TRUE
opt_select(ctg_clipped) <- "xyz"
opt_output_files(ctg_clipped) <- paste0("E:/Remote Sensing Media/", date_folder,"/06. Point Clouds Normalised/", "{ORIGINALFILENAME}_classified_normalised")

# 2. Load the Master Baseline DTM VRT generated in October 2025
baseline_dtm <- rast("E:/Remote Sensing Media/03. 30 October 2025/05. Baseline Plot DTMs/Master_Baseline_DTM_03._30_October_2025.vrt")

tic()
print("Normalizing point clouds against the October baseline DTM...")

# 3. Perform raster-based elevation normalization
# lidR automatically queries the DTM raster to subtract ground elevation from Z coordinates
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

ctg_normalised <- readLAScatalog(paste0("E:/Remote Sensing Media/",date_folder,"/06. Point Clouds Normalised"))

# Limit the amount of workers (threads) if RAM constrained
plan(multisession, workers = 6)
opt_independent_files(ctg_normalised) <- TRUE
opt_select(ctg_normalised) <- "xyz"

# Filter noise: Drop ground points and sub-surface/extreme-height noise
opt_filter(ctg_normalised) <- "-drop_z_below 0 -drop_z_above 30"

tic()
# Stream outputs directly to disk rather than holding in memory
opt_output_files(ctg_normalised) <- paste0("E:/Remote Sensing Media/",date_folder,"/07. Canopy Height Models/", "{*}_chm")

# Rasterize highest points into a continuous CHM grid with interpolation
ctg_chm <- rasterize_canopy(ctg_normalised,
                            res = 0.05,
                            algorithm = p2r(na.fill = tin()))
print("Rasterize canopy time:")
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 8. Master CHM Consolidation & Metric Extraction ####
# ──────────────────────────────────────────────────────────────────────────────
tic()
print("Consolidating individual CHM plots into a master single file...")

# 1. Gather all individual CHM files generated in Step 7
chm_dir <- paste0("E:/Remote Sensing Media/", date_folder, "/07. Canopy Height Models/")
chm_files <- list.files(chm_dir, pattern = "\\.tif$", full.names = TRUE)

# 2. Create a Virtual Raster (VRT) to stitch them together
site_chm_vrt <- terra::vrt(chm_files)

# 3. Write the VRT out to a single physical .tif file
single_chm_path <- paste0(chm_dir, "Master_Site_CHM_Single_", file_date_safe, ".tif")
terra::writeRaster(site_chm_vrt, filename = single_chm_path, overwrite = TRUE)

print(paste("Master CHM successfully saved to:", single_chm_path))

# 4. Re-read the singular physical file for the extraction step
site_chm_single <- terra::rast(single_chm_path)

# 5. Calculate exact maximum tree height within each delineated crown polygon
trees$Tree_Height <- exact_extract(site_chm_single, trees, 'max')

# 6. Save to shapefile (completely overwriting old files to prevent schema errors)
st_write(trees, paste0("E:/Remote Sensing Media/",date_folder,"/09. Crown Metrics/Crown_Metrics_", file_date_safe, ".shp"), delete_dsn = TRUE)

# 7. Save lightweight tabular CSV data (drops the messy spatial geometry text)
write.csv(st_drop_geometry(trees), paste0("E:/Remote Sensing Media/",date_folder,"/09. Crown Metrics/Crown_Metrics_", file_date_safe, ".csv"), row.names = FALSE)
toc()