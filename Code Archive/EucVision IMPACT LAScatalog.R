# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: LiDAR POINT CLOUD PROCESSING & CANOPY HEIGHT EXTRACTION PIPELINE ####
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
date_folder <- "22. 08 April 2026"

# Extract the date part and create a safe filename format
# (e.g., "17. 02 March 2026" -> "02_March_2026")
file_date <- sub("^\\d+\\.\\s*", "", date_folder)
file_date_safe <- gsub(" ", "_", file_date)

# --- CREATE MISSING IMPACT DIRECTORIES ---
# Define all required output directories for the IMPACT workflow
impact_dirs <- c(
  paste0("E:/Remote Sensing Media/", date_folder, "/04. Point Clouds Clipped IMPACT"),
  paste0("E:/Remote Sensing Media/", date_folder, "/05. Point Clouds Ground Classified IMPACT"),
  paste0("E:/Remote Sensing Media/", date_folder, "/06. Point Clouds Normalised IMPACT"),
  paste0("E:/Remote Sensing Media/", date_folder, "/07. Canopy Height Models IMPACT"),
  paste0("E:/Remote Sensing Media/", date_folder, "/09. Crown Metrics IMPACT"),
  paste0("E:/Remote Sensing Media/", date_folder, "/10. Digital Terrain Models")
)

# Loop through and create any directories that do not already exist
for (dir in impact_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("Created missing directory:", dir, "\n")
  }
}

# ──────────────────────────────────────────────────────────────────────────────
# 3. Spatial Data Loading & Boundary Definition ####
# ──────────────────────────────────────────────────────────────────────────────
# --- BATCH PROCESSING SETUP ---
las_folder <- paste0("E:/Remote Sensing Media/", date_folder, "/03. Point Clouds")
all_las_files <- list.files(las_folder, pattern = "\\.(las|laz)$", full.names = TRUE, ignore.case = TRUE)

# Identify Top and Bottom point cloud files based on filenames
top_files <- all_las_files[grepl("Top", basename(all_las_files), ignore.case = TRUE)]
bot_files <- all_las_files[grepl("Bottom", basename(all_las_files), ignore.case = TRUE)]

# Load boundary shapefiles for clipping the Top and Bottom point clouds
IMPACT_Top <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LAScatalog Boundaries/Top IMPACT Boundaries.shp")
IMPACT_Bottom <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LAScatalog Boundaries/Bottom IMPACT Boundaries.shp")

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

# ──────────────────────────────────────────────────────────────────────────────
# 4. Point Cloud Clipping (IMPACT Site) ####
# ──────────────────────────────────────────────────────────────────────────────
plan(multisession) # Enable parallel processing
tic()

if (length(top_files) > 0 && length(bot_files) > 0) {
  print("Top and Bottom point clouds detected. Cropping separately...")
  
  # Initialize separate LiDAR catalogs
  ctg_top <- readLAScatalog(top_files)
  ctg_bot <- readLAScatalog(bot_files)
  
  # Configure engine settings for TOP cropping
  opt_independent_files(ctg_top) <- FALSE
  opt_select(ctg_top) <- "xyz"
  opt_output_files(ctg_top) <- paste0("E:/Remote Sensing Media/",date_folder,"/04. Point Clouds Clipped IMPACT/IMPACT_Site_Top_", file_date_safe)
  
  # Configure engine settings for BOTTOM cropping
  opt_independent_files(ctg_bot) <- FALSE
  opt_select(ctg_bot) <- "xyz"
  opt_output_files(ctg_bot) <- paste0("E:/Remote Sensing Media/",date_folder,"/04. Point Clouds Clipped IMPACT/IMPACT_Site_Bottom_", file_date_safe)
  
  # Execute spatial clipping operations
  print("Cropping Top IMPACT Boundary...")
  ctg_clipped_top <- clip_roi(ctg_top, IMPACT_Top)
  
  print("Cropping Bottom IMPACT Boundary...")
  ctg_clipped_bot <- clip_roi(ctg_bot, IMPACT_Bottom)
  
  # Re-read the fully clipped directory as a single unified catalog for downstream processing
  ctg_clipped <- readLAScatalog(paste0("E:/Remote Sensing Media/", date_folder, "/04. Point Clouds Clipped IMPACT"))
  
} else {
  print("Single point cloud or no Top/Bottom distinction detected. Make sure your files contain 'Top' and 'Bottom' in the names.")
}

toc()

# ──────────────────────────────────────────────────────────────────────────────
# 5. Ground Classification ####
# ──────────────────────────────────────────────────────────────────────────────
plan(multisession, workers = 6) # Scale up to 6 workers for heavy processing
opt_independent_files(ctg_clipped) <- TRUE
opt_select(ctg_clipped) <- "xyz"

# CRITICAL FOR WHOLE SITE PROCESSING: Define spatial chunks and buffers
# This splits the large catalog across CPU workers and prevents edge artifacts
opt_chunk_size(ctg_clipped) <- 200  # Process in 200m x 200m chunks
opt_chunk_buffer(ctg_clipped) <- 10 # 10m overlap buffer around chunks

tic()
opt_output_files(ctg_clipped) <- paste0("E:/Remote Sensing Media/",date_folder,"/05. Point Clouds Ground Classified IMPACT/", "IMPACT_Tile_{XLEFT}_{YBOTTOM}_classified", file_date_safe)

# Apply Cloth Simulation Filter (CSF) to identify ground points
ctg_classified <- classify_ground(
  ctg_clipped, 
  csf(sloop_smooth = TRUE, 
      class_threshold = 0.15, 
      cloth_resolution = 1.5, 
      rigidness = 3,
      time_step = 1)
)
toc()

# ──────────────────────────────────────────────────────────────────────────────
# Optional: Generate Digital Terrain Model (DTM) ####
# ──────────────────────────────────────────────────────────────────────────────

# Ensure output folder exists (e.g., ".../10. Digital Terrain Models/") before running.

ctg_classified <- readLAScatalog(paste0("E:/Remote Sensing Media/",date_folder,"/05. Point Clouds Ground Classified IMPACT"))
plan(multisession, workers = 6)
opt_independent_files(ctg_classified) <- TRUE
opt_select(ctg_classified) <- "xyzc"

# Maintain exact chunking/buffer parameters to match classification
opt_chunk_size(ctg_classified) <- 200
opt_chunk_buffer(ctg_classified) <- 10

tic()
opt_output_files(ctg_classified) <- paste0("E:/Remote Sensing Media/",date_folder,"/10. Digital Terrain Models/", "IMPACT_Tile_{XLEFT}_{YBOTTOM}_dtm_", file_date_safe)

# Rasterize ground points into DTM using TIN algorithm.
# Default 0.5m resolution smooths micro-noise; adjust to 0.05m to match CHM pixel scaling.
ctg_dtm <- rasterize_terrain(ctg_classified,
                             res = 0.5,
                             algorithm = tin())
print("Rasterize DTM time:")
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 6. Height Normalization ####
# ──────────────────────────────────────────────────────────────────────────────
plan(multisession, workers = 6)
opt_independent_files(ctg_classified) <- TRUE
opt_select(ctg_classified) <- "xyzc"

# Maintain consistent chunking logic
opt_chunk_size(ctg_classified) <- 200 
opt_chunk_buffer(ctg_classified) <- 10 

tic()
opt_output_files(ctg_classified) <- paste0("E:/Remote Sensing Media/",date_folder,"/06. Point Clouds Normalised IMPACT/", "IMPACT_Tile_{XLEFT}_{YBOTTOM}_normalised", file_date_safe)

# Normalize point elevations to calculate absolute tree heights above ground
ctg_normalised <- normalize_height(las = ctg_classified, algorithm = tin())
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 7. Canopy Height Model (CHM) Generation ####
# ──────────────────────────────────────────────────────────────────────────────
# Re-load normalized catalog to ensure a clean state
ctg_normalised <- readLAScatalog(paste0("E:/Remote Sensing Media/",date_folder,"/06. Point Clouds Normalised IMPACT"))

plan(multisession, workers = 4)
opt_independent_files(ctg_normalised) <- TRUE
opt_select(ctg_normalised) <- "xyz"

# Filter noise: Drop extreme outlier points (below ground or impossibly high)
opt_filter(ctg_normalised) <- "-drop_z_below 0 -drop_z_above 30"

# Maintain consistent chunking logic
opt_chunk_size(ctg_normalised) <- 200 
opt_chunk_buffer(ctg_normalised) <- 10 

tic()
opt_output_files(ctg_normalised) <- paste0("E:/Remote Sensing Media/",date_folder,"/07. Canopy Height Models IMPACT/", "IMPACT_Tile_{XLEFT}_{YBOTTOM}_chm", file_date_safe)

# Rasterize highest points into a continuous CHM grid
ctg_chm <- rasterize_canopy(ctg_normalised,
                            res = 0.05,
                            algorithm = p2r(na.fill = tin()))
print("Rasterize canopy time:")
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 8. Metric Extraction & Export ####
# ──────────────────────────────────────────────────────────────────────────────
tic()

# Consolidate individual chunked CHM tiles via Virtual Raster (VRT)
chm_files <- list.files(paste0("E:/Remote Sensing Media/",date_folder,"/07. Canopy Height Models IMPACT/"), 
                        pattern = "\\.tif$", full.names = TRUE)
site_chm_vrt <- terra::vrt(chm_files)

# Write the VRT out to a single physical .tif file for exact_extract
single_chm_path <- paste0("E:/Remote Sensing Media/",date_folder,"/07. Canopy Height Models IMPACT/IMPACT_Site_CHM_Single_", file_date_safe, ".tif")
terra::writeRaster(site_chm_vrt, filename = single_chm_path, overwrite = TRUE)

# Optional: Cleanup chunk files to save SSD space
file.remove(chm_files)

# Re-read the singular physical file for geospatial extraction
site_chm_single <- terra::rast(single_chm_path)

# Calculate exact maximum tree height within each delineated crown polygon
trees$Tree_Height <- exact_extract(site_chm_single, trees, 'max')

# Export updated geospatial shapefile
st_write(trees, paste0("E:/Remote Sensing Media/",date_folder,"/09. Crown Metrics IMPACT/Crown_Metrics_", file_date_safe, ".shp"), delete_dsn = TRUE)

# Export lightweight tabular CSV data (drops heavy geometry)
write.csv(st_drop_geometry(trees), paste0("E:/Remote Sensing Media/",date_folder,"/09. Crown Metrics IMPACT/Crown_Metrics_", file_date_safe, ".csv"), row.names = FALSE)
toc()