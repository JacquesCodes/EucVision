# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: BATCH CHM HEIGHT EXTRACTION PIPELINE
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Project: EucXylo (https://eucxylo.sun.ac.za/) 
# ──────────────────────────────────────────────────────────────────────────────
# Description: Automates the extraction of maximum canopy heights from temporal
#              CHMs using individual tree crown polygons. It loops through all
#              valid datasets, extracts metrics, cleans attributes, and ensures
#              the output shapefiles retain the strict EPSG:2048 CRS label.
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
library(sf)            # Spatial vector data handling
library(terra)         # Spatial raster data handling
library(dplyr)         # Data wrangling and piping logic
library(exactextractr) # Fast raster metric extraction (zonal statistics)
library(tictoc)        # Script execution timing 

# The strict, exact OGC Well-Known Text (WKT) for EPSG:2048 to force the South/West axis fix
pure_epsg_2048_wkt <- 'PROJCS["Hartebeesthoek94 / Lo19",GEOGCS["Hartebeesthoek94",DATUM["Hartebeesthoek94",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6148"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4148"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",19],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Southing",SOUTH],AXIS["Westing",WEST],AUTHORITY["EPSG","2048"]]'

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration & Batch Management ####
# ──────────────────────────────────────────────────────────────────────────────
base_dir <- "E:/Remote Sensing Media"

# --- RUN CONTROLS ---
# Set to a specific folder name to run only that dataset (e.g., "01. 25 February 2025") 
# Set to NULL to run the full batch process.
target_date_override <- "31. 26 June 2026 Oblique"

# Folders to ignore during the batch processing loop
exclude_list <- c("000. Projects",
                  "00. Baseline DTM",
                  "00. Dataset Template",
                  "01. 25 February 2025", 
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

print("Starting Batch Metric Extraction Pipeline...")

# ──────────────────────────────────────────────────────────────────────────────
# MASTER BATCH LOOP START ####
# ──────────────────────────────────────────────────────────────────────────────
for (folder_path in dataset_folders) {
  
  # Dynamically extract and format the date from the folder name
  date_folder <- basename(folder_path)
  
  # Extract the date part and create a safe filename format
  # (e.g., "17. 02 March 2026" -> "02_March_2026")
  file_date_safe <- gsub(" ", "_", sub("^\\d+\\.\\s*", "", date_folder))
  
  print(paste("================================================================"))
  print(paste("EXTRACTING METRICS FOR:", date_folder))
  print(paste("================================================================"))
  
  # Define paths
  polygons_dir <- file.path(folder_path, "08. Crown Polygons")
  chm_dir      <- file.path(folder_path, "07. Canopy Height Models")
  metrics_dir  <- file.path(folder_path, "09. Crown Metrics")
  
  # Ensure the output directory exists
  if (!dir.exists(metrics_dir)) dir.create(metrics_dir, recursive = TRUE)
  
  # Define specific input files
  path_trees <- file.path(polygons_dir, paste0("Crown_Polygons_", file_date_safe, ".shp"))
  path_chm   <- file.path(chm_dir, paste0("Master_Site_CHM_Single_", file_date_safe, ".tif"))
  
  # Define final output paths for the metrics
  output_path_shp <- file.path(metrics_dir, paste0("Crown_Metrics_", file_date_safe, ".shp"))
  output_path_csv <- file.path(metrics_dir, paste0("Crown_Metrics_", file_date_safe, ".csv"))
  
  # --- Logic Gates for Skipping Iterations ---
  if (!file.exists(path_trees)) {
    print(paste("-> SKIPPED: No Crown Polygons found for", date_folder))
    next
  }
  if (!file.exists(path_chm)) {
    print(paste("-> SKIPPED: No Master CHM found for", date_folder))
    next
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # 3. Spatial Data Loading & Metric Extraction ####
  # ────────────────────────────────────────────────────────────────────────────
  tic("Metric extraction complete")
  
  # Load polygons and clean columns
  trees <- st_read(path_trees, quiet = TRUE)
  
  # Force the Tree column to exactly 2 decimal places and drop unnecessary columns
  trees <- trees %>%
    mutate(Tree = round(Tree, 2)) %>%
    select(-any_of(c("grop_ld", "group_ulid", "N_GM", "id", "N_FG", "N_BG", "BBox")))
  
  # Load CHM
  ctg_chm <- rast(path_chm)
  
  # Ensure Shapefiles' CRS are seen as EPSG:2048 in R
  st_crs(trees) <- st_crs(ctg_chm)
  
  # Extract exact max tree height
  print("Running exact_extract spatial analysis...")
  trees$Tree_Height <- exact_extract(ctg_chm, trees, 'max')
  
  # ────────────────────────────────────────────────────────────────────────────
  # 4. Export & EPSG:2048 CRS Injection ####
  # ────────────────────────────────────────────────────────────────────────────
  print("Exporting data and enforcing EPSG:2048...")
  
  # Export lightweight tabular CSV data
  write.csv(st_drop_geometry(trees), output_path_csv, row.names = FALSE)
  
  # Export geospatial shapefile
  st_write(trees, output_path_shp, delete_dsn = TRUE, quiet = TRUE)
  
  # Bypassing the GDAL driver to forcefully inject the EPSG:2048 WKT
  prj_path <- sub("\\.shp$", ".prj", output_path_shp, ignore.case = TRUE)
  writeLines(pure_epsg_2048_wkt, prj_path)
  
  toc()
  
  # Clean memory for the next iteration
  rm(trees, ctg_chm)
  gc()
}

print("================================================================")
print("BATCH EXTRACTION COMPLETE!")
print("================================================================")