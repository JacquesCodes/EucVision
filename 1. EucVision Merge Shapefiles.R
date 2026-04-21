# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: SPATIAL DATA MERGING & CSV INTEGRATION PIPELINE ####
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
# Load required libraries for spatial data manipulation and data wrangling
library(sf)
library(dplyr)

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration & Path Management ####
# ──────────────────────────────────────────────────────────────────────────────
# === CONFIGURE PATHS ===
# Change this single variable for each new batch to process the correct folder!
date_folder <- "21. 31 March 2026"

# Extract the date part by removing the leading folder number, dot, and space
# (e.g., This turns "17. 02 March 2026" into "02 March 2026")
file_date <- sub("^\\d+\\.\\s*", "", date_folder)

# Replace spaces with underscores for safer file naming conventions
# (e.g., This turns "02 March 2026" into "02_March_2026")
file_date_safe <- gsub(" ", "_", file_date)

# Define the base directories (These stay constant across batches)
base_input <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/03. QGIS Extracted data"
base_output <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/04. QGIS Combined Output"
csv_path <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data analysis/00. Dataset template.csv"

# Dynamically construct the full input and output paths for the current batch
input_folder <- file.path(base_input, date_folder)
output_folder <- file.path(base_output, date_folder)

# Define final output paths for both internal backups and external drives
output_shp <- file.path(output_folder, paste0("Crown_Polygons_", file_date_safe, ".shp")) # Teams backup
base_output_file <- paste0("E:/Remote Sensing Media/", date_folder, "/08. Crown Polygons/Crown_Polygons_", file_date_safe, ".shp") # SSD drive

# ──────────────────────────────────────────────────────────────────────────────
# 3. Directory Preparation ####
# ──────────────────────────────────────────────────────────────────────────────
# Create the Teams backup output directory if it doesn't already exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Create the directory for the SSD drive if it doesn't exist
ssd_dir <- dirname(base_output_file)
if (!dir.exists(ssd_dir)) {
  dir.create(ssd_dir, recursive = TRUE)
}

# ──────────────────────────────────────────────────────────────────────────────
# 4. Shapefile Loading & Sorting ####
# ──────────────────────────────────────────────────────────────────────────────
# Fetch a list of all shapefiles in the target input folder
shp_files <- list.files(input_folder, pattern = "\\.shp$", full.names = TRUE)

# Safety check to prevent the script from running on empty directories
if (length(shp_files) == 0) {
  stop("No shapefiles found in the input directory. Check your path:\n", input_folder)
}

# --- SORT FILES NUMERICALLY ---
# Crucial Step: standard list.files sorts lexically (Plot 1, Plot 10, Plot 2). 
# We must sort numerically (Plot 1, Plot 2 ... Plot 10) to align properly with the CSV structure.
plot_numbers <- as.numeric(gsub("\\D", "", basename(shp_files)))
shp_files <- shp_files[order(plot_numbers)]

cat("Found", length(shp_files), "shapefiles in", date_folder, "Starting merge...\n")

# ──────────────────────────────────────────────────────────────────────────────
# 5. Spatial Data Merging ####
# ──────────────────────────────────────────────────────────────────────────────
# Read and combine all individual plot shapefiles into one continuous spatial dataframe
combined_sf <- lapply(shp_files, function(file) {
  
  # Read the shapefile silently to keep console output clean
  temp_sf <- st_read(file, quiet = TRUE)
  
  # Extract the plot name from the filename (e.g., "Plot 1" from "Plot 1.shp")
  plot_name <- tools::file_path_sans_ext(basename(file))
  
  # Add the plot name as a new column to trace the origin of the spatial features
  # Renamed to 'Plot_shp' to avoid naming conflicts with the incoming CSV 'Plot' column
  temp_sf <- temp_sf %>% mutate(Plot_shp = plot_name) 
  
  return(temp_sf)
  
}) %>% 
  # bind_rows safely aggregates spatial data even if column structures vary slightly
  bind_rows() 

# ──────────────────────────────────────────────────────────────────────────────
# 6. CSV Template Integration ####
# ──────────────────────────────────────────────────────────────────────────────
# Load the master dataset template
csv_data <- read.csv(csv_path)

# Data validation check to ensure shapefiles and CSV align 1:1
if(nrow(csv_data) != nrow(combined_sf)) {
  warning("The number of rows in the CSV (", nrow(csv_data), 
          ") does not match the merged shapefile (", nrow(combined_sf), ").")
}

# Bind the standard CSV columns to the left side of the newly merged shapefile data
# We recast the final object back to 'sf' to ensure spatial properties are preserved
combined_sf <- bind_cols(csv_data, combined_sf) %>% st_as_sf()

# ──────────────────────────────────────────────────────────────────────────────
# 7. Export Processed Data ####
# ──────────────────────────────────────────────────────────────────────────────
# Write the finalized shapefile to both the Teams folder and the SSD
# append = FALSE ensures clean overwriting of existing files without throwing errors
st_write(combined_sf, output_shp, append = FALSE)
st_write(combined_sf, base_output_file, append = FALSE)

cat("✅ All done. CSV successfully attached. Merged shapefile saved to:\n", output_shp, "\n")