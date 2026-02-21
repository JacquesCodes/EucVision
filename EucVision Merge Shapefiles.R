# Load required libraries
library(sf)
library(dplyr)

# === CONFIGURE PATHS ===
# Change this single variable for each new batch!
date_folder <- "13. 29 January 2026"

# Define the base directories (these stay constant)
base_input <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/03. QGIS Extracted data"
base_output <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/04. QGIS Combined Output"


# Dynamically construct the full paths
input_folder <- file.path(base_input, date_folder)
output_folder <- file.path(base_output, date_folder)
output_shp <- file.path(output_folder, "All_Plots.shp")
base_output_file <- paste0("E:/Remote Sensing Media/",date_folder,"/08. Crown shape file/All_Plots.shp")


# Create the output directory if it doesn't already exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Get a list of all shapefiles in the input folder
shp_files <- list.files(input_folder, pattern = "\\.shp$", full.names = TRUE)

if (length(shp_files) == 0) {
  stop("No shapefiles found in the input directory. Check your path:\n", input_folder)
}

cat("Found", length(shp_files), "shapefiles in", date_folder, "Starting merge...\n")

# Read and combine all shapefiles
combined_sf <- lapply(shp_files, function(file) {
  
  # Read the shapefile silently
  temp_sf <- st_read(file, quiet = TRUE)
  
  # Extract the plot name from the filename (e.g., "Plot 1" from "Plot 1.shp")
  plot_name <- tools::file_path_sans_ext(basename(file))
  
  # Add the plot name as a new column to keep track of the source
  temp_sf <- temp_sf %>% mutate(Plot = plot_name)
  
  return(temp_sf)
  
}) %>% 
  # bind_rows handles cases where some shapefiles might have slightly different columns
  bind_rows() 

# Write the merged shapefile to the output directory
# delete_layer = TRUE allows it to overwrite if the file already exists
st_write(combined_sf, output_shp, delete_layer = TRUE)
st_write(combined_sf, base_output_file, delete_layer = TRUE)

cat("✅ All done. Merged shapefile successfully saved to:\n", output_shp, "\n")