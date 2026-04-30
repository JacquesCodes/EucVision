# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: TLS POINT CLOUD PREPROCESSING & SPATIAL CORRECTION PIPELINE ####
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────
# Description: Automates the preprocessing and spatial correction of raw 
#              Terrestrial Laser Scanning (TLS) point clouds (.laz). It 
#              specifically addresses South African CRS axis-flipping conflicts 
#              by applying a brute-force coordinate inversion (X/Y) to properly 
#              align the data with the standard Lo19 (EPSG:2048) projection. 
#              The script iteratively reprojects, corrects, re-headers, and 
#              exports uncompressed .las files strictly formatted for immediate 
#              integration into the main plot-level processing workflow.
# ──────────────────────────────────────────────────────────────────────────────

library(lidR)
library(sf)

# 1. Define base directories and reference files
input_dir <- "E:/Remote Sensing Media/07. December 2025 (TLS)/03. Point Clouds"
output_dir <- "E:/Remote Sensing Media/07. December 2025 (TLS)/04. Point Clouds Clipped"
file_april <- "E:/Remote Sensing Media/23. 13 April 2026/03. Point Clouds/Top Sector Cross Hatch_13 April 2026_densified_point_cloud.las"

# Create the output directory if it doesn't already exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 2. Extract target CRS from April data
message("Extracting target CRS from April reference data...")
header_april <- readLASheader(file_april)
crs_april <- st_crs(header_april)

# 3. Setup Naming Variables (Mimicking your pipeline logic)
date_folder <- "07. December 2025 (TLS)"
file_date <- sub("^\\d+\\.\\s*", "", date_folder)
file_date_safe <- gsub(" ", "_", file_date) # Results in "December_2025_(TLS)"

# 4. Get the list of all 21 .laz files
laz_files <- list.files(input_dir, pattern = "\\.laz$", full.names = TRUE)

# 5. Loop through and process each file
for (file in laz_files) {
  
  # Extract the Plot ID (removes the folder path and the .laz extension)
  base_name <- basename(file)
  plot_id <- sub("\\.laz$", "", base_name) 
  
  # Construct the pipeline-compliant output filename
  output_filename <- paste0("Plot_", plot_id, "_", file_date_safe, ".las")
  output_path <- file.path(output_dir, output_filename)
  
  # Skip processing if the file already exists (useful if the script gets interrupted)
  if (file.exists(output_path)) {
    message("Skipping (already exists): ", output_filename)
    next
  }
  
  message("\n========================================================")
  message("Processing Plot ID: ", plot_id)
  message("========================================================")
  
  # Load the compressed .laz file
  cloud <- readLAS(file)
  
  # Reproject to standard Lo19 (EPSG:2048)
  message("  -> Reprojecting to standard Lo19...")
  cloud_reprojected <- st_transform(cloud, 2048)
  
  # BRUTE FORCE: Extract raw data to bypass bounding box offset limits
  message("  -> Extracting raw data and flipping coordinates...")
  las_data <- cloud_reprojected@data
  
  # Manually flip X and Y axes to match the negative Easting/Northing orientation
  las_data$X <- las_data$X * -1
  las_data$Y <- las_data$Y * -1
  
  # Rebuild the point cloud to recalculate valid offsets
  message("  -> Rebuilding LAS object...")
  cloud_fixed <- LAS(las_data)
  
  # Overwrite the CRS metadata
  st_crs(cloud_fixed) <- crs_april
  
  # Save the uncompressed, renamed file
  message("  -> Saving as ", output_filename, "...")
  writeLAS(cloud_fixed, output_path)
  
  # CRITICAL: Clean up RAM to prevent the loop from crashing
  rm(cloud, cloud_reprojected, cloud_fixed, las_data)
  gc()
}

message("\n========================================================")
message("BATCH CONVERSION COMPLETE!")
message("========================================================")