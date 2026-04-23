# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: BATCH PLOT-LEVEL LiDAR PROCESSING PIPELINE USING DTM ####
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
library(stringr)

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration & Batch Management ####
# ──────────────────────────────────────────────────────────────────────────────
base_dir <- "E:/Remote Sensing Media"

# Define the absolute path to your Baseline DTM (Stays constant across batches)
baseline_dtm_path <- "E:/Remote Sensing Media/03. 30 October 2025/05b. Baseline Plot DTMs/Master_Baseline_DTM_Smoothed_30_October_2025.tif"

# Load the baseline DTM once into memory before the loop starts to save time
baseline_dtm <- rast(baseline_dtm_path)

# --- EXCLUDE LIST ---
# Clearly mark any folders you want the batch processor to completely ignore
exclude_list <- c("000. Projects",
                  "00. Dataset template", 
                  "01. 25 February 2025", 
                  "17. 03 March 2026 (Multispectral)",
                  "20. 24 March 2026 (Multispectral)")

# Fetch all folders and filter for standard dated folders not on the exclude list
folders <- list.dirs(base_dir, recursive = FALSE)
dataset_folders <- folders[grepl("^\\d{2}\\.", basename(folders)) & !basename(folders) %in% exclude_list]

# ──────────────────────────────────────────────────────────────────────────────
# 3. Static Spatial Data Loading ####
# ──────────────────────────────────────────────────────────────────────────────
# Load plot boundary shapefiles once (Static across all batches)
plots_buffered_unsorted <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LAScatalog Boundaries/LAScatalog Plot Boundaries.shp")
plots <- plots_buffered_unsorted[order(plots_buffered_unsorted$id), ]

# Enable parallel processing globally for the session
plan(multisession, workers = 6) 

print("Starting batch processing pipeline...")

# ──────────────────────────────────────────────────────────────────────────────
# 4. BATCH PROCESSING LOOP ####
# ──────────────────────────────────────────────────────────────────────────────
for (folder_path in dataset_folders) {
  
  # Extract dynamic date strings for this specific iteration
  date_folder <- basename(folder_path)
  file_date <- sub("^\\d+\\.\\s*", "", date_folder)
  file_date_safe <- gsub(" ", "_", file_date)
  
  print(paste("================================================================"))
  print(paste("PROCESSING DATASET:", date_folder))
  print(paste("================================================================"))
  
  # Define processing directories dynamically
  las_folder     <- file.path(folder_path, "03. Point Clouds")
  clipped_dir    <- file.path(folder_path, "04. Point Clouds Clipped")
  normalised_dir <- file.path(folder_path, "06. Point Clouds Normalised")
  chm_dir        <- file.path(folder_path, "07. Canopy Height Models")
  polygons_dir   <- file.path(folder_path, "08. Crown Polygons")
  metrics_dir    <- file.path(folder_path, "09. Crown Metrics")
  
  # Ensure all output directories exist for this batch
  for (dir in c(clipped_dir, normalised_dir, chm_dir, metrics_dir)) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  }
  
  # --- PRE-FLIGHT CHECKS ---
  # Skip this batch if no point clouds exist yet
  all_las_files <- list.files(las_folder, pattern = "\\.(las|laz)$", full.names = TRUE, ignore.case = TRUE)
  if (length(all_las_files) == 0) {
    print(paste("-> SKIPPED: No .las/.laz files found in", las_folder))
    next
  }
  
  # Skip this batch if the Crown Polygons shapefile hasn't been created yet
  crown_shp_path <- file.path(polygons_dir, paste0("Crown_Polygons_", file_date_safe, ".shp"))
  if (!file.exists(crown_shp_path)) {
    print(paste("-> SKIPPED: No Crown Polygons found for", date_folder))
    next
  }
  
  # Load and validate dynamic crown polygons
  trees <- st_read(crown_shp_path, quiet = TRUE)
  if (is.na(st_crs(trees)$epsg) || st_crs(trees)$epsg != 2048) {
    trees <- st_transform(trees, 2048)
  }
  trees <- trees %>% select(-any_of(c("group_ulid", "N_GM", "id", "N_FG", "N_BG", "BBox")))
  
  # ────────────────────────────────────────────────────────────────────────────
  # 4.1. Point Cloud Cropping
  # ────────────────────────────────────────────────────────────────────────────
  tic("Cropping complete")
  top_files <- all_las_files[grepl("Top", basename(all_las_files), ignore.case = TRUE)]
  bot_files <- all_las_files[grepl("Bottom", basename(all_las_files), ignore.case = TRUE)]
  
  if (length(top_files) > 0 && length(bot_files) > 0) {
    print("Splitting processing by Plot ID (Top/Bottom)...")
    plots_top <- plots %>% filter(id <= 21)
    plots_bot <- plots %>% filter(id >= 22)
    
    ctg_top <- readLAScatalog(top_files)
    ctg_bot <- readLAScatalog(bot_files)
    
    opt_independent_files(ctg_top) <- FALSE
    opt_select(ctg_top) <- "xyz"
    opt_output_files(ctg_top) <- file.path(clipped_dir, paste0("Plot_{id}_", file_date_safe))
    
    opt_independent_files(ctg_bot) <- FALSE
    opt_select(ctg_bot) <- "xyz"
    opt_output_files(ctg_bot) <- file.path(clipped_dir, paste0("Plot_{id}_", file_date_safe))
    
    suppressMessages(clip_roi(ctg_top, plots_top))
    suppressMessages(clip_roi(ctg_bot, plots_bot))
    
  } else {
    print("Processing single catalog entirely...")
    ctg <- readLAScatalog(las_folder)
    opt_independent_files(ctg) <- FALSE
    opt_select(ctg) <- "xyz"
    opt_output_files(ctg) <- file.path(clipped_dir, paste0("Plot_{id}_", file_date_safe))
    
    suppressMessages(clip_roi(ctg, plots))
  }
  ctg_clipped <- readLAScatalog(clipped_dir)
  toc()
  
  # ────────────────────────────────────────────────────────────────────────────
  # 4.2. Height Normalization (USING BASELINE DTM)
  # ────────────────────────────────────────────────────────────────────────────
  tic("Normalization complete")
  print("Normalizing point clouds against the October baseline DTM...")
  opt_independent_files(ctg_clipped) <- TRUE
  opt_select(ctg_clipped) <- "xyz"
  opt_output_files(ctg_clipped) <- file.path(normalised_dir, "{ORIGINALFILENAME}_classified_normalised")
  
  ctg_normalised <- normalize_height(las = ctg_clipped, algorithm = baseline_dtm)
  toc()
  
  # ────────────────────────────────────────────────────────────────────────────
  # 4.3. Canopy Height Model (CHM) Generation
  # ────────────────────────────────────────────────────────────────────────────
  tic("CHM Rasterization complete")
  print("Generating individual Canopy Height Models...")
  ctg_normalised <- readLAScatalog(normalised_dir)
  opt_independent_files(ctg_normalised) <- TRUE
  opt_select(ctg_normalised) <- "xyz"
  opt_filter(ctg_normalised) <- "-drop_z_below 0 -drop_z_above 30"
  opt_output_files(ctg_normalised) <- file.path(chm_dir, "{*}_chm")
  
  ctg_chm <- rasterize_canopy(ctg_normalised, res = 0.05, algorithm = p2r(na.fill = tin()))
  toc()
  
  # ────────────────────────────────────────────────────────────────────────────
  # 4.4. Master CHM Consolidation & Dynamic Metric Extraction
  # ────────────────────────────────────────────────────────────────────────────
  tic("Metric extraction complete")
  print("Extracting metrics and consolidating dynamic master CHM...")
  
  chm_files <- list.files(chm_dir, pattern = "\\.tif$", full.names = TRUE)
  site_chm_vrt <- terra::vrt(chm_files)
  
  # Extract exact maximum tree heights directly from the VRT
  trees$Tree_Height <- exact_extract(site_chm_vrt, trees, 'max')
  
  # Calculate dynamic cap safely
  max_tree_height <- max(trees$Tree_Height[is.finite(trees$Tree_Height)], na.rm = TRUE)
  dynamic_cap <- ceiling(max_tree_height)
  print(paste("-> Dynamic CHM cap safely set to:", dynamic_cap, "meters"))
  
  # Clamp the artifacts dynamically
  site_chm_clamped <- terra::clamp(site_chm_vrt, lower = 0, upper = dynamic_cap)
  
  # Write the perfectly capped raster out
  single_chm_path <- file.path(chm_dir, paste0("Master_Site_CHM_Single_", file_date_safe, ".tif"))
  terra::writeRaster(site_chm_clamped, filename = single_chm_path, overwrite = TRUE)
  
  # Save Shapefile and CSV
  st_write(trees, file.path(metrics_dir, paste0("Crown_Metrics_", file_date_safe, ".shp")), delete_dsn = TRUE, quiet = TRUE)
  write.csv(st_drop_geometry(trees), file.path(metrics_dir, paste0("Crown_Metrics_", file_date_safe, ".csv")), row.names = FALSE)
  toc()
  
  # --- GARBAGE COLLECTION ---
  # Clear large catalog objects from RAM before the next iteration starts
  rm(ctg_clipped, ctg_normalised, ctg_chm, site_chm_vrt, site_chm_clamped, trees, chm_files, all_las_files)
  gc()
}

print("================================================================")
print("BATCH PIPELINE COMPLETE! All datasets processed successfully.")
print("================================================================")