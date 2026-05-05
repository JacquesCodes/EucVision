# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: SPATIAL DATA MERGING & CSV INTEGRATION PIPELINE
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────
# Description: Aggregates individual spatial plot shapefiles into a single 
#              continuous spatial dataframe.
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
library(sf)      # Core package for spatial vector operations
library(dplyr)   # Core package for data wrangling and piping
library(tictoc)  # For tracking script execution time

tic()
print("Initiating Spatial Data Merging & CSV Integration...")

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration & Path Management ####
# ──────────────────────────────────────────────────────────────────────────────
# === CONFIGURE PATHS ===
# Change this single variable for each new batch to process the correct folder!
date_folder <- "12. 22 January 2026"

# String manipulation for file naming
# Removes the leading folder number (e.g., "17. 02 March 2026" -> "02 March 2026")
file_date <- sub("^\\d+\\.\\s*", "", date_folder)
# Replaces spaces with underscores (e.g., "02 March 2026" -> "02_March_2026")
file_date_safe <- gsub(" ", "_", file_date)

# Define the base directories (Constant across batches)
base_input <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/03. QGIS Extracted data"
base_output <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/04. QGIS Combined Output"
csv_path <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data analysis/00. Dataset template.csv"

# Dynamically construct the full input and output paths for the current batch
input_folder <- file.path(base_input, date_folder)
output_folder <- file.path(base_output, date_folder)

# Define final output paths for both internal Teams backups and the external SSD
output_shp <- file.path(output_folder, paste0("Crown_Polygons_", file_date_safe, ".shp")) 
base_output_file <- paste0("E:/Remote Sensing Media/", date_folder, "/08. Crown Polygons/Crown_Polygons_", file_date_safe, ".shp") 

# ──────────────────────────────────────────────────────────────────────────────
# 3. Directory Preparation ####
# ──────────────────────────────────────────────────────────────────────────────
print("Verifying output directories...")

# Ensure the Teams backup output directory exists
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# Ensure the external SSD output directory exists
ssd_dir <- dirname(base_output_file)
if (!dir.exists(ssd_dir)) dir.create(ssd_dir, recursive = TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# 4. Shapefile Loading & Sorting ####
# ──────────────────────────────────────────────────────────────────────────────
print("Scanning for shapefiles and enforcing numeric sorting...")

# Fetch a list of all shapefiles in the target input folder
shp_files <- list.files(input_folder, pattern = "\\.shp$", full.names = TRUE)

# Safety check to prevent the script from running on empty directories
if (length(shp_files) == 0) {
  stop("CRITICAL ERROR: No shapefiles found in the input directory. Check path:\n", input_folder)
}

# --- CRITICAL: SORT FILES NUMERICALLY ---
# Standard list.files sorts lexically (Plot 1, Plot 10, Plot 2). 
# We must sort numerically (Plot 1, Plot 2 ... Plot 10) to align 1:1 with the CSV.
plot_numbers <- as.numeric(gsub("\\D", "", basename(shp_files)))
shp_files <- shp_files[order(plot_numbers)]

print(paste("Found", length(shp_files), "shapefiles in", date_folder, "- Starting merge..."))

# ──────────────────────────────────────────────────────────────────────────────
# 5. Spatial Data Merging ####
# ──────────────────────────────────────────────────────────────────────────────
# Read and combine all individual plot shapefiles into one continuous sf dataframe
combined_sf <- lapply(shp_files, function(file) {
  
  # Read the shapefile silently to keep the console output clean
  temp_sf <- st_read(file, quiet = TRUE)
  
  # Extract the raw plot name from the filename
  plot_name <- tools::file_path_sans_ext(basename(file))
  
  # Add the plot name as a new column to trace feature origins
  # Tagged as 'Plot_shp' to avoid namespace conflicts with the incoming CSV
  temp_sf <- temp_sf %>% mutate(Plot_shp = plot_name) 
  
  return(temp_sf)
  
}) %>% 
  # safely aggregates spatial data even if column structures vary slightly
  bind_rows() 

# ──────────────────────────────────────────────────────────────────────────────
# 6. CSV Template Integration ####
# ──────────────────────────────────────────────────────────────────────────────
print("Attaching master CSV template to spatial geometries...")

# Load the master dataset template
csv_data <- read.csv(csv_path)

# Data validation check to ensure shapefiles and CSV align 1:1 before binding
if(nrow(csv_data) != nrow(combined_sf)) {
  warning("MISMATCH DETECTED: The number of rows in the CSV (", nrow(csv_data), 
          ") does not match the merged shapefile (", nrow(combined_sf), ").")
}

# Bind the standard CSV columns to the left side of the newly merged shapefile data
# We recast the final object back to 'sf' to guarantee spatial properties survive
combined_sf <- bind_cols(csv_data, combined_sf) %>% st_as_sf()

# ──────────────────────────────────────────────────────────────────────────────
# 7. Export Processed Data & Enforce CRS ####
# ──────────────────────────────────────────────────────────────────────────────
print("Writing finalized shapefiles to disk...")

# Write the finalized shapefile to both the Teams folder and the SSD
# append = FALSE ensures clean overwriting of existing files without throwing errors
st_write(combined_sf, output_shp, append = FALSE, quiet = TRUE)
st_write(combined_sf, base_output_file, append = FALSE, quiet = TRUE)

# --- INJECT PURE EPSG:2048 WKT INTO .PRJ FILES ---
print("Enforcing strict EPSG:2048 CRS on output .prj files...")

# The strict, exact OGC Well-Known Text (WKT) for EPSG:2048
pure_epsg_2048_wkt <- 'PROJCS["Hartebeesthoek94 / Lo19",GEOGCS["Hartebeesthoek94",DATUM["Hartebeesthoek94",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6148"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4148"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",19],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Southing",SOUTH],AXIS["Westing",WEST],AUTHORITY["EPSG","2048"]]'

# Dynamically generate the .prj file paths by replacing the .shp extension
prj_teams <- sub("\\.shp$", ".prj", output_shp, ignore.case = TRUE)
prj_ssd <- sub("\\.shp$", ".prj", base_output_file, ignore.case = TRUE)

# Overwrite the newly created ESRI-style .prj files with the strict EPSG string
writeLines(pure_epsg_2048_wkt, prj_teams)
writeLines(pure_epsg_2048_wkt, prj_ssd)

print(paste("- Teams Backup successfully saved and labeled to:", output_shp))
print(paste("- External SSD successfully saved and labeled to:", base_output_file))
print("Pipeline Complete!")
toc()