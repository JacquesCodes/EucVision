# Load required libraries
library(sf)
library(dplyr)

# === CONFIGURE PATHS ===
# Change this single variable for each new batch!
date_folder <- "17. 02 March 2026"

# Define the base directories (these stay constant)
base_input <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/03. QGIS Extracted data"
base_output <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/04. QGIS Combined Output"
csv_path <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data analysis/00. Dataset template.csv"

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

# --- SORT FILES NUMERICALLY ---
# Crucial: standard list.files sorts Plot 1, Plot 10, Plot 2. 
# We need Plot 1, Plot 2 ... Plot 10 so it aligns with the CSV properly!
plot_numbers <- as.numeric(gsub("\\D", "", basename(shp_files)))
shp_files <- shp_files[order(plot_numbers)]

cat("Found", length(shp_files), "shapefiles in", date_folder, "Starting merge...\n")

# Read and combine all shapefiles
combined_sf <- lapply(shp_files, function(file) {
  
  # Read the shapefile silently
  temp_sf <- st_read(file, quiet = TRUE)
  
  # Extract the plot name from the filename (e.g., "Plot 1" from "Plot 1.shp")
  plot_name <- tools::file_path_sans_ext(basename(file))
  
  # Add the plot name as a new column to keep track of the source
  temp_sf <- temp_sf %>% mutate(Plot_shp = plot_name) # Renamed to Plot_shp to avoid overlap with CSV 'Plot' column
  
  return(temp_sf)
  
}) %>% 
  # bind_rows handles cases where some shapefiles might have slightly different columns
  bind_rows() 

# --- LOAD CSV AND BIND TO LEFT ---
# Read the CSV file
csv_data <- read.csv(csv_path)

# Check if the number of rows matches
if(nrow(csv_data) != nrow(combined_sf)) {
  warning("The number of rows in the CSV (", nrow(csv_data), 
          ") does not match the merged shapefile (", nrow(combined_sf), ").")
}

# Bind the CSV columns to the left side of the shapefile data
# We re-cast it back to an 'sf' object to keep the geospatial attributes functional
combined_sf <- bind_cols(csv_data, combined_sf) %>% st_as_sf()

# Write the merged shapefile to the output directory
# delete_layer = TRUE allows it to overwrite if the file already exists
st_write(combined_sf, output_shp, delete_layer = TRUE)
st_write(combined_sf, base_output_file, delete_layer = TRUE)

cat("✅ All done. CSV attached. Merged shapefile successfully saved to:\n", output_shp, "\n")