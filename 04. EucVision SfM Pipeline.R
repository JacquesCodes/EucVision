# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: BATCH PLOT-LEVEL SfM PROCESSING PIPELINE
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/) 
# ──────────────────────────────────────────────────────────────────────────────
# Description: Automates the plot-level processing of multiple temporal SfM 
#              point cloud datasets. It iteratively clips raw point clouds to 
#              plot boundaries, normalizes them against an ultimate baseline DTM, 
#              generates Canopy Height Models (CHMs), extracts maximum tree 
#              heights, and enforces strict EPSG:2048 CRS on metric outputs.
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
library(lidR)            # Core point cloud manipulation and processing engine
library(sf)              # Spatial vector data handling (reading/writing shapefiles)
library(terra)           # Spatial raster data handling (DTMs, CHMs, clamping)
library(dplyr)           # Data wrangling and piping logic for dataframes
library(future)          # Parallel processing framework for chunked lidR operations
library(tictoc)          # Script execution timing 
library(exactextractr)   # Fast raster metric extraction (zonal statistics)
library(stringr)         # String manipulation for file paths
library(gstat)           # Spatial interpolation (often used in DTM/CHM generation)
library(geometry)        # Mesh and polygon geometry functions
library(sp)              # Legacy spatial framework (required by some older lidR functions)

# The strict, exact OGC Well-Known Text (WKT) for EPSG:2048 to force the South/West axis fix
pure_epsg_2048_wkt <- 'PROJCS["Hartebeesthoek94 / Lo19",GEOGCS["Hartebeesthoek94",DATUM["Hartebeesthoek94",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6148"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4148"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",19],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Southing",SOUTH],AXIS["Westing",WEST],AUTHORITY["EPSG","2048"]]'

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration & Batch Management ####
# ──────────────────────────────────────────────────────────────────────────────
base_dir <- "E:/Remote Sensing Media"
baseline_dtm_path <- "E:/Remote Sensing Media/00. Baseline DTM/Ultimate_Ensemble_Baseline_DTM.tif"

# Load the ultimate baseline surface model into memory as a SpatRaster
# (We will use this as the MASTER CRS truth for all vector layers)
baseline_dtm <- rast(baseline_dtm_path)

# --- RUN CONTROLS ---
# Set to a specific folder name to run only that dataset (e.g., "01. 25 February 2025") 
# Set to NULL to run the full batch process.
# target_date_override <- NULL
target_date_override <- "24. 23 April 2026"

# Disk Space Management: TRUE retains intermediate point clouds, FALSE deletes them.
keep_intermediate_dirs <- FALSE

# --- EXCLUDE LIST ---
# Folders to ignore during the batch processing loop
exclude_list <- c("000. Projects",
                  "00. Baseline DTM",
                  "00. Dataset Template", 
                  "01. 25 February 2025",
                  "07. December 2025 (TLS)",
                  "17. 03 March 2026 (Multispectral)",
                  "20. 24 March 2026 (Multispectral)")

# Scan the base directory and filter for valid date folders
folders <- list.dirs(base_dir, recursive = FALSE)
dataset_folders <- folders[grepl("^\\d{2}\\.", basename(folders)) & !basename(folders) %in% exclude_list]

# Apply the target date override if one was explicitly provided
if (!is.null(target_date_override)) {
  dataset_folders <- dataset_folders[basename(dataset_folders) == target_date_override]
  if (length(dataset_folders) == 0) {
    stop("CRITICAL ERROR: Target date folder not found! Please check spelling.")
  }
}

# ──────────────────────────────────────────────────────────────────────────────
# 3. Static Spatial Data Loading ####
# ──────────────────────────────────────────────────────────────────────────────
print("Loading static spatial data and initializing parallel processing...")

# st_read loads the vector data. quiet = TRUE suppresses messy console output.
plots_buffered_unsorted <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. QGIS Shapefiles/1. LAScatalog Plot Boundaries/LAScatalog Plot Boundaries.shp", quiet = TRUE)

# Force the clipping polygons to perfectly inherit the master DTM's spatial grid (CRS)
st_crs(plots_buffered_unsorted) <- st_crs(baseline_dtm)

# Enforce numeric sorting by the 'id' column to guarantee consistent processing order
plots <- plots_buffered_unsorted[order(plots_buffered_unsorted$id), ]

# Enable parallel processing for lidR operations to run across 6 CPU logical processors.
# Note: Users should adjust 'workers' based on their available CPU logical processors and RAM.
plan(multisession, workers = 6)

print("Starting processing pipeline...")

# ──────────────────────────────────────────────────────────────────────────────
# MASTER BATCH LOOP START ####
# ──────────────────────────────────────────────────────────────────────────────
for (folder_path in dataset_folders) {
  
  # ────────────────────────────────────────────────────────────────────────────
  # 4. Dataset Initialization & Checking ####
  # ────────────────────────────────────────────────────────────────────────────
  # Dynamically extract and format the date from the folder name
  date_folder <- basename(folder_path)
  
  # Extract the date part and create a safe filename format
  # (e.g., "17. 02 March 2026" -> "02_March_2026")
  file_date <- sub("^\\d+\\.\\s*", "", date_folder)
  file_date_safe <- gsub(" ", "_", file_date)
  
  print(paste("================================================================"))
  print(paste("PROCESSING DATASET:", date_folder))
  print(paste("================================================================"))
  
  # Map out all input/output directories for the current specific date
  las_folder     <- file.path(folder_path, "03. Point Clouds")
  clipped_dir    <- file.path(folder_path, "04. Point Clouds Clipped")
  normalised_dir <- file.path(folder_path, "06. Point Clouds Normalised")
  chm_dir        <- file.path(folder_path, "07. Canopy Height Models")
  polygons_dir   <- file.path(folder_path, "08. Crown Polygons")
  metrics_dir    <- file.path(folder_path, "09. Crown Metrics")
  
  # Define the target path for the final output of this iteration
  single_chm_path <- file.path(chm_dir, paste0("Master_Site_CHM_Single_", file_date_safe, ".tif"))
  
  # dir.create safely generates the folder structure if it doesn't exist yet
  for (dir in c(clipped_dir, normalised_dir, chm_dir, metrics_dir)) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  }
  
  # --- Logic Gates for Skipping Iterations ---
  # If the final CHM already exists, skip to the next folder to save time
  if (file.exists(single_chm_path)) {
    print(paste("-> SKIPPED: Master CHM already exists for", date_folder))
    next
  }
  
  # list.files looks for anything ending in .las or .laz. If empty, skip.
  all_las_files <- list.files(las_folder, pattern = "\\.(las|laz)$", full.names = TRUE, ignore.case = TRUE)
  if (length(all_las_files) == 0) {
    print(paste("-> SKIPPED: No .las/.laz files found in", las_folder))
    next
  }
  
  # We require the Crown Polygons shapefile to extract metrics later. If missing, skip.
  crown_shp_path <- file.path(polygons_dir, paste0("Crown_Polygons_", file_date_safe, ".shp"))
  if (!file.exists(crown_shp_path)) {
    print(paste("-> SKIPPED: No Crown Polygons found for", date_folder))
    next
  }
  
  # select(-any_of()) strips out legacy tracking columns to keep the final dataframe clean
  trees <- st_read(crown_shp_path, quiet = TRUE) %>% 
    mutate(Tree = round(Tree, 2)) %>%  # Force 2 decimal precision ---
    select(-any_of(c("group_ulid", "N_GM", "id", "N_FG", "N_BG", "BBox")))
  
  # --- APPLIED CRS FIX: THE CROWN POLYGONS ---
  # Force the tree polygons to inherit the master grid immediately upon loading
  st_crs(trees) <- st_crs(baseline_dtm)
  
  # ────────────────────────────────────────────────────────────────────────────
  # 5. Point Cloud Cropping ####
  # ────────────────────────────────────────────────────────────────────────────
  tic("Cropping complete")
  
  # grepl searches the file names for "Top" or "Bottom" to split heavy site workloads
  top_files <- all_las_files[grepl("Top", basename(all_las_files), ignore.case = TRUE)]
  bot_files <- all_las_files[grepl("Bottom", basename(all_las_files), ignore.case = TRUE)]
  
  if (length(top_files) > 0 && length(bot_files) > 0) {
    print("Splitting processing by Plot ID (Top/Bottom)...")
    
    # filter() subsets our spatial boundaries into top (<=21) and bottom (>=22) plots
    plots_top <- plots %>% filter(id <= 21)
    plots_bot <- plots %>% filter(id >= 22)
    
    # readLAScatalog treats the files as a continuous virtual dataset without loading them into RAM
    ctg_top <- readLAScatalog(top_files)
    ctg_bot <- readLAScatalog(bot_files)
    
    # opt_independent_files = FALSE tells lidR the files are parts of a whole, so it should resolve overlapping points
    opt_independent_files(ctg_top) <- FALSE
    # opt_select = "xyz" strips out color, intensity, and classification data to drastically reduce memory usage
    opt_select(ctg_top) <- "xyz"
    # Defines the naming convention for the clipped outputs using the 'id' attribute from the plot polygons
    opt_output_files(ctg_top) <- file.path(clipped_dir, paste0("Plot_{id}_", file_date_safe))
    
    opt_independent_files(ctg_bot) <- FALSE
    opt_select(ctg_bot) <- "xyz"
    opt_output_files(ctg_bot) <- file.path(clipped_dir, paste0("Plot_{id}_", file_date_safe))
    
    # clip_roi cookie-cuts the point cloud using the polygon boundaries. suppressMessages hides lidR's verbose output.
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
  
  # Load the newly clipped directory as our working catalog for the next step
  ctg_clipped <- readLAScatalog(clipped_dir)
  toc()
  
  # ────────────────────────────────────────────────────────────────────────────
  # 6. Height Normalization ####
  # ────────────────────────────────────────────────────────────────────────────
  tic("Normalization complete")
  print("Normalizing point clouds against the Ultimate Baseline DTM...")
  
  # opt_independent_files = TRUE treats each clipped plot as isolated. No need to look at overlaps anymore.
  opt_independent_files(ctg_clipped) <- TRUE
  opt_select(ctg_clipped) <- "xyz"
  # {ORIGINALFILENAME} preserves the "Plot_X_Date" structure generated in step 5
  opt_output_files(ctg_clipped) <- file.path(normalised_dir, "{ORIGINALFILENAME}_classified_normalised")
  
  # normalize_height subtracts the Z-value of the baseline DTM from the Z-value of every point, 
  # effectively flattening the topography so ground is at Z = 0 and canopy is absolute height.
  ctg_normalised <- normalize_height(las = ctg_clipped, algorithm = baseline_dtm)
  toc()
  
  # ────────────────────────────────────────────────────────────────────────────
  # 7. Canopy Height Model (CHM) Generation ####
  # ────────────────────────────────────────────────────────────────────────────
  tic("CHM Rasterization complete")
  print("Generating individual Canopy Height Models...")
  
  ctg_normalised <- readLAScatalog(normalised_dir)
  opt_independent_files(ctg_normalised) <- TRUE
  opt_select(ctg_normalised) <- "xyz"
  
  # opt_filter applies a strict height threshold *while* loading data. 
  # Drops points below 0 (sub-surface noise) and above 30m (birds, extreme noise).
  opt_filter(ctg_normalised) <- "-drop_z_below 0 -drop_z_above 30"
  opt_output_files(ctg_normalised) <- file.path(chm_dir, "{*}_chm")
  
  # rasterize_canopy builds the 2.5D surface. 
  # res = 0.05 creates 5cm pixels. 
  # algorithm = p2r(na.fill = tin()) uses a point-to-raster approach and interpolates empty pixels using a Triangular Irregular Network.
  ctg_chm <- rasterize_canopy(ctg_normalised, res = 0.05, algorithm = p2r(na.fill = tin()))
  toc()
  
  # ────────────────────────────────────────────────────────────────────────────
  # 8. Metric Extraction & Consolidation ####
  # ────────────────────────────────────────────────────────────────────────────
  tic("Metric extraction complete")
  print("Extracting metrics and consolidating dynamic master CHM...")
  
  # list.files gathers the paths to all the individual plot CHMs generated in step 7
  chm_files <- list.files(chm_dir, pattern = "\\.tif$", full.names = TRUE)
  # terra::vrt stitches them together logically as a Virtual Raster without duplicating files on disk
  site_chm_vrt <- terra::vrt(chm_files)
  
  # exact_extract calculates the "max" pixel value (tallest canopy point) falling inside each tree polygon
  # It is vastly faster and more memory-efficient than standard GIS extract functions.
  # (No st_crs fix needed here because we safely applied it when the trees were loaded in Section 4!)
  trees$Tree_Height <- exact_extract(site_chm_vrt, trees, 'max')
  
  # is.finite() strips out completely empty polygons that might have returned NA or Inf.
  # max() finds the absolute tallest tree across the whole site to determine our clamping ceiling.
  max_tree_height <- max(trees$Tree_Height[is.finite(trees$Tree_Height)], na.rm = TRUE)
  dynamic_cap <- ceiling(max_tree_height)
  print(paste("-> Dynamic CHM cap safely set to:", dynamic_cap, "meters"))
  
  # terra::clamp forcefully cuts off any remaining spike artifacts in the CHM above our dynamic cap
  site_chm_clamped <- terra::clamp(site_chm_vrt, lower = 0, upper = dynamic_cap)
  
  # writeRaster writes out the final, stitched, clamped CHM image
  terra::writeRaster(site_chm_clamped, filename = single_chm_path, overwrite = TRUE)
  
  # Define output paths for metric files
  out_shp_path <- file.path(metrics_dir, paste0("Crown_Metrics_", file_date_safe, ".shp"))
  out_csv_path <- file.path(metrics_dir, paste0("Crown_Metrics_", file_date_safe, ".csv"))
  
  # st_write exports the updated polygons (now containing Tree_Height). delete_dsn = TRUE allows overwriting.
  st_write(trees, out_shp_path, delete_dsn = TRUE, quiet = TRUE)
  
  # --- INJECT PURE EPSG:2048 WKT INTO .PRJ FILE ---
  prj_path <- sub("\\.shp$", ".prj", out_shp_path, ignore.case = TRUE)
  writeLines(pure_epsg_2048_wkt, prj_path)
  
  # st_drop_geometry strips the spatial mapping data so the dataframe can be written to a clean, lightweight CSV
  write.csv(st_drop_geometry(trees), out_csv_path, row.names = FALSE)
  toc()
  
  # ────────────────────────────────────────────────────────────────────────────
  # 9. Housekeeping & Memory Management ####
  # ────────────────────────────────────────────────────────────────────────────
  if (!keep_intermediate_dirs) {
    print("Cleaning up intermediate directories to save drive space...")
    # unlink(..., recursive = TRUE) is R's equivalent of "delete folder and all contents"
    unlink(clipped_dir, recursive = TRUE)
    unlink(normalised_dir, recursive = TRUE)
  } else {
    print("Retaining intermediate directories as requested...")
  }
  
  # rm() removes heavy R objects from the environment workspace.
  # gc() (Garbage Collection) forces the operating system to reclaim that RAM before the next loop iteration.
  rm(ctg_clipped, ctg_normalised, ctg_chm, site_chm_vrt, site_chm_clamped, trees, chm_files, all_las_files)
  gc()
}

print("================================================================")
print("PIPELINE COMPLETE! All designated datasets processed successfully.")
print("================================================================")