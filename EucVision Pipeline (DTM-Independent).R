# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: PLOT-LEVEL LiDAR POINT CLOUD PROCESSING PIPELINE USING LAS ####
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────
# Description: Automates the plot-level processing of raw point clouds.
#              Unlike the SfM pipeline, this script is DTM-independent.
#              The pipeline handles spatial clipping, advanced ground 
#              classification using Progressive TIN Densification (PTD) tuned 
#              for TLS/SfM, and point-cloud-based height normalization. It 
#              subsequently generates high-resolution Canopy Height Models 
#              (CHMs) and extracts absolute maximum tree heights within 
#              pre-defined crown polygons for downstream metric analysis.
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
date_folder <- "07. December 2025 (TLS)"

# Extract the date part and create a safe filename format
# (e.g., "17. 02 March 2026" -> "02_March_2026")
file_date <- sub("^\\d+\\.\\s*", "", date_folder)
file_date_safe <- gsub(" ", "_", file_date)

# Define all processing directories dynamically based on the date
las_folder     <- paste0("E:/Remote Sensing Media/", date_folder, "/03. Point Clouds/")
clipped_dir    <- paste0("E:/Remote Sensing Media/", date_folder, "/04. Point Clouds Clipped/")
classified_dir <- paste0("E:/Remote Sensing Media/", date_folder, "/05. Point Clouds Ground Classified/")
normalised_dir <- paste0("E:/Remote Sensing Media/", date_folder, "/06. Point Clouds Normalised/")
chm_dir        <- paste0("E:/Remote Sensing Media/", date_folder, "/07. Canopy Height Models/")
polygons_dir   <- paste0("E:/Remote Sensing Media/", date_folder, "/08. Crown Polygons/")
metrics_dir    <- paste0("E:/Remote Sensing Media/", date_folder, "/09. Crown Metrics/")

# Ensure all output directories exist before processing begins
for (dir in c(clipped_dir, classified_dir, normalised_dir, chm_dir, metrics_dir)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# ──────────────────────────────────────────────────────────────────────────────
# 3. Spatial Data Loading & Boundary Definition ####
# ──────────────────────────────────────────────────────────────────────────────
# Load individual plot boundary shapefiles and sort them chronologically/spatially by ID
plots_buffered_unsorted <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. QGIS Shapefiles/1. LAScatalog Plot Boundaries/LAScatalog Plot Boundaries.shp")
plots <- plots_buffered_unsorted[order(plots_buffered_unsorted$id), ]

plan(multisession, workers = 2)

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

  # Split plot boundaries by ID to align with the physical flight paths
  plots_top <- plots %>% filter(id <= 21)
  plots_bot <- plots %>% filter(id >= 22)

  # Initialize separate LiDAR catalogs
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

  # Execute spatial clipping operations
  print("Cropping Top plots (1-21)...")
  ctg_clipped_top <- clip_roi(ctg_top, plots_top)

  print("Cropping Bottom plots (22-75)...")
  ctg_clipped_bot <- clip_roi(ctg_bot, plots_bot)

  # Re-read the fully clipped directory as a single unified catalog for downstream processing
  ctg_clipped <- readLAScatalog(clipped_dir)

} else {
  print("Single point cloud or no Top/Bottom distinction detected. Processing entirely...")

  # Standard processing for a single or unstructured catalog
  ctg <- readLAScatalog(las_folder)
  opt_independent_files(ctg) <- FALSE
  opt_select(ctg) <- "xyz"
  opt_output_files(ctg) <- paste0(clipped_dir, "Plot_{id}_", file_date_safe)

  # Execute spatial clipping operations
  ctg_clipped <- clip_roi(ctg, plots)
}

toc()

# ──────────────────────────────────────────────────────────────────────────────
# 5. Ground Classification ####
# ──────────────────────────────────────────────────────────────────────────────

ctg_clipped <- readLAScatalog(clipped_dir)

opt_independent_files(ctg_clipped) <- TRUE
opt_select(ctg_clipped) <- "xyz"

tic()
# Stream outputs directly to disk rather than holding in memory
opt_output_files(ctg_clipped) <- paste0(classified_dir, "{*}_classified")

# PTD: Tuned specifically for dense TLS point clouds in clipped small plots
ctg_classified <- classify_ground(
  ctg_clipped,            
  ptd(res = 1.5,         # Small enough to get plenty of seeds across a 6m width, large enough to bypass thick trunks.
      angle = 20,        # Relaxed to allow the TIN to capture natural terrain slopes and micro-topography.
      distance = 0.2,    # Tightened significantly to prevent the TIN from jumping into low shrubs or tall grass.
      spacing = 0.1)     # (Optional) Halts densification when triangles hit 10cm, preserving fine TLS terrain details.
)

# # PTD: Tuned specifically for dense SfM canopy crusts for whole IMPACT OAL site
# ctg_classified <- classify_ground(
#   ctg_clipped,           
#   ptd(res = 12,          # 12m guarantees the search grid overhangs canopy widths
#       angle = 6,         # Strict 6-degree limit prevents the TIN from climbing canopy walls
#       distance = 1.5)    # Keeps the vertical step distance tight
# )

# LidR creator recommends using the PTD algorithm exclusively from March 2026.

# # SCF: Tuned Simulated cloth function for SfM individual plots 
# ctg_classified <- classify_ground(
#   ctg_clipped,           
#   csf(sloop_smooth = TRUE, 
#       class_threshold = 0.01, 
#       cloth_resolution = 1, 
#       time_step = 1))

# SCF Parameter Configuration Notes:
# Class_threshold: The distance to the simulated cloth to classify a point into ground/non-ground. 
#  - Default is 0.5. Must be smaller than the smallest tree. 0.01 preferred for best height estimations.
#  - The higher the value, the higher the ground classifications become.
# Cloth_resolution: The distance between particles in the simulated cloth.
#  - Default is 0.5. PREFERRED = 1 for 1m x 1m plots. 
#  - Do not lower below 1; otherwise, the cloth falls between points under closed canopy and classifies trees as ground.
#  - Above 1m seems to cut off half of a plot

toc()

# ──────────────────────────────────────────────────────────────────────────────
# 6. Height Normalization (Using Classified LAS) ####
# ──────────────────────────────────────────────────────────────────────────────

opt_independent_files(ctg_classified) <- TRUE
# Load only x-, y-, z- coordinates and the classification ("c") values into RAM
opt_select(ctg_classified) <- "xyzc"

tic()
# Stream outputs directly to disk rather than holding in memory
opt_output_files(ctg_classified) <- paste0(normalised_dir, "{*}_normalised")

# Perform point cloud-based elevation normalization without a pre-existing raster
ctg_normalised <- normalize_height(las = ctg_classified, algorithm = tin())
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 7. Canopy Height Model (CHM) Generation ####
# ──────────────────────────────────────────────────────────────────────────────
ctg_normalised <- readLAScatalog(normalised_dir)

# Limit the amount of workers (threads) if RAM constrained
# plan(multisession, workers = 6)
opt_independent_files(ctg_normalised) <- TRUE
opt_select(ctg_normalised) <- "xyz"

# Filter noise: Drop ground points and sub-surface/extreme-height noise
opt_filter(ctg_normalised) <- "-drop_z_below 0 -drop_z_above 30"

tic()
# Stream outputs directly to disk rather than holding in memory
opt_output_files(ctg_normalised) <- paste0(chm_dir, "{*}_chm")

# Rasterize highest points into a continuous CHM grid with interpolation
ctg_chm <- rasterize_canopy(ctg_normalised,
                            res = 0.05,
                            algorithm = p2r(na.fill = tin()))
print("Rasterize canopy time:")
toc()

# ──────────────────────────────────────────────────────────────────────────────
# 8. Metric Extraction & Export ####
# ──────────────────────────────────────────────────────────────────────────────
tic()

# Calculate exact maximum tree height within each delineated crown polygon
trees$Tree_Height <- exact_extract(ctg_chm, trees, 'max')

# Save to shapefile (completely overwriting old files to prevent schema errors)
st_write(trees, paste0(metrics_dir, "Crown_Metrics_", file_date_safe, ".shp"), delete_dsn = TRUE)

# Save lightweight tabular CSV data (drops the messy spatial geometry text)
write.csv(st_drop_geometry(trees), paste0(metrics_dir, "Crown_Metrics_", file_date_safe, ".csv"), row.names = FALSE)
toc()