# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: BATCH PLOT-LEVEL SfM PROCESSING PIPELINE ####
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
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
baseline_dtm_path <- "E:/Remote Sensing Media/00. Baseline DTM/Ultimate_Ensemble_Baseline_DTM.tif"

baseline_dtm <- rast(baseline_dtm_path)

# --- RUN CONTROLS ---
# Set to a specific folder name to run only that dataset (e.g., "01. 25 February 2025") 
# Set to NULL to run the full batch process.
target_date_override <- NULL 

# Set to TRUE to keep intermediate directories (Clipped, Normalised). FALSE deletes them.
keep_intermediate_dirs <- FALSE 

# --- EXCLUDE LIST ---
exclude_list <- c("000. Projects",
                  "00. Baseline DTM",
                  "00. Dataset Template", 
                  "01. 25 February 2025",
                  "17. 03 March 2026 (Multispectral)",
                  "20. 24 March 2026 (Multispectral)")

folders <- list.dirs(base_dir, recursive = FALSE)
dataset_folders <- folders[grepl("^\\d{2}\\.", basename(folders)) & !basename(folders) %in% exclude_list]

# Apply the target date override if one is provided
if (!is.null(target_date_override)) {
  dataset_folders <- dataset_folders[basename(dataset_folders) == target_date_override]
  if (length(dataset_folders) == 0) {
    stop("Target date folder not found! Please check the spelling and try again.")
  }
}

# ──────────────────────────────────────────────────────────────────────────────
# 3. Static Spatial Data Loading ####
# ──────────────────────────────────────────────────────────────────────────────
plots_buffered_unsorted <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LAScatalog Boundaries/LAScatalog Plot Boundaries.shp")
plots <- plots_buffered_unsorted[order(plots_buffered_unsorted$id), ]

plan(multisession, workers = 6) 

print("Starting processing pipeline...")

# ──────────────────────────────────────────────────────────────────────────────
# 4. PROCESSING LOOP ####
# ──────────────────────────────────────────────────────────────────────────────
for (folder_path in dataset_folders) {
  
  date_folder <- basename(folder_path)
  file_date <- sub("^\\d+\\.\\s*", "", date_folder)
  file_date_safe <- gsub(" ", "_", file_date)
  
  print(paste("================================================================"))
  print(paste("PROCESSING DATASET:", date_folder))
  print(paste("================================================================"))
  
  las_folder     <- file.path(folder_path, "03. Point Clouds")
  clipped_dir    <- file.path(folder_path, "04. Point Clouds Clipped")
  normalised_dir <- file.path(folder_path, "06. Point Clouds Normalised")
  chm_dir        <- file.path(folder_path, "07. Canopy Height Models")
  polygons_dir   <- file.path(folder_path, "08. Crown Polygons")
  metrics_dir    <- file.path(folder_path, "09. Crown Metrics")
  
  single_chm_path <- file.path(chm_dir, paste0("Master_Site_CHM_Single_", file_date_safe, ".tif"))
  
  for (dir in c(clipped_dir, normalised_dir, chm_dir, metrics_dir)) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  }
  
  if (file.exists(single_chm_path)) {
    print(paste("-> SKIPPED: Master CHM already exists for", date_folder))
    next
  }
  
  all_las_files <- list.files(las_folder, pattern = "\\.(las|laz)$", full.names = TRUE, ignore.case = TRUE)
  if (length(all_las_files) == 0) {
    print(paste("-> SKIPPED: No .las/.laz files found in", las_folder))
    next
  }
  
  crown_shp_path <- file.path(polygons_dir, paste0("Crown_Polygons_", file_date_safe, ".shp"))
  if (!file.exists(crown_shp_path)) {
    print(paste("-> SKIPPED: No Crown Polygons found for", date_folder))
    next
  }
  
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
  # 4.2. Height Normalization 
  # ────────────────────────────────────────────────────────────────────────────
  tic("Normalization complete")
  print("Normalizing point clouds against the Ensemble Baseline DTM...")
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
  
  trees$Tree_Height <- exact_extract(site_chm_vrt, trees, 'max')
  
  max_tree_height <- max(trees$Tree_Height[is.finite(trees$Tree_Height)], na.rm = TRUE)
  dynamic_cap <- ceiling(max_tree_height)
  print(paste("-> Dynamic CHM cap safely set to:", dynamic_cap, "meters"))
  
  site_chm_clamped <- terra::clamp(site_chm_vrt, lower = 0, upper = dynamic_cap)
  
  terra::writeRaster(site_chm_clamped, filename = single_chm_path, overwrite = TRUE)
  
  st_write(trees, file.path(metrics_dir, paste0("Crown_Metrics_", file_date_safe, ".shp")), delete_dsn = TRUE, quiet = TRUE)
  write.csv(st_drop_geometry(trees), file.path(metrics_dir, paste0("Crown_Metrics_", file_date_safe, ".csv")), row.names = FALSE)
  toc()
  
  # --- FOLDER CLEANUP ---
  if (!keep_intermediate_dirs) {
    print("Cleaning up intermediate directories to save drive space...")
    unlink(clipped_dir, recursive = TRUE)
    unlink(normalised_dir, recursive = TRUE)
  } else {
    print("Retaining intermediate directories as requested...")
  }
  
  # --- GARBAGE COLLECTION ---
  rm(ctg_clipped, ctg_normalised, ctg_chm, site_chm_vrt, site_chm_clamped, trees, chm_files, all_las_files)
  gc()
}

print("================================================================")
print("PIPELINE COMPLETE! All designated datasets processed successfully.")
print("================================================================")