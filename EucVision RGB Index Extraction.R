# Load necessary libraries
library(terra)
library(sf)
library(exactextractr)
library(dplyr)
library(stringr)
library(tidyr)

# --- CONFIGURATION ---
base_dir <- "E:/Remote Sensing Media"

# Get a list of all subdirectories in the base folder
folders <- list.dirs(base_dir, recursive = FALSE)

# 1. Explicitly exclude the non-dataset folders
exclude_list <- c("00. Dataset template", "000. Projects")
folders <- folders[!basename(folders) %in% exclude_list]

# 2. Filter for only the numbered dataset folders (e.g., "01. 25 February 2025")
dataset_folders <- folders[grepl("^\\d{2}\\.", basename(folders))]

all_extractions <- list()

print("Starting batch extraction...")

# --- BATCH PROCESSING LOOP ---
for (folder_path in dataset_folders) {
  
  folder_name <- basename(folder_path)
  date_str <- str_replace(folder_name, "^\\d{2}\\.\\s*", "")
  
  print("========================================")
  print(paste("Processing Date:", date_str))
  
  # Construct paths
  ortho_dir <- file.path(folder_path, "01. Orthomosaics")
  shp_dir <- file.path(folder_path, "08. Crown shape file", folder_name)
  
  # Find the .tif and .shp files
  all_tifs <- list.files(ortho_dir, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
  shp_files <- list.files(shp_dir, pattern = "\\.shp$", full.names = TRUE, ignore.case = TRUE)
  
  # Filter out the "Cross" flights so we only use the primary radiometric data
  valid_tifs <- all_tifs[!grepl("Cross", basename(all_tifs), ignore.case = TRUE)]
  
  if (length(valid_tifs) == 0 || length(shp_files) == 0) {
    print(paste("  -> Skipping: Missing TIF or SHP in", folder_name))
    next
  }
  
  # Process each valid TIF
  for (tif_file in valid_tifs) {
    print(paste("  -> Loading Raster:", basename(tif_file)))
    
    # Check if this specific TIF is a Top or Bottom split
    is_top <- grepl("Top", basename(tif_file), ignore.case = TRUE)
    is_bottom <- grepl("Bottom", basename(tif_file), ignore.case = TRUE)
    
    ortho <- rast(tif_file)
    
    # Standard RGB Band mapping (Band 1=Red, Band 2=Green, Band 3=Blue)
    R <- ortho[[1]]
    G <- ortho[[2]]
    B <- ortho[[3]]
    
    print("  -> Calculating VARI, GLI, and TGI...")
    VARI <- (G - R) / (G + R - B)
    GLI  <- ((2 * G) - R - B) / ((2 * G) + R + B)
    TGI  <- G - 0.39 * R - 0.61 * B
    
    index_stack <- c(VARI, GLI, TGI)
    names(index_stack) <- c("VARI", "GLI", "TGI")
    
    # Process Shapefiles
    for (shp_file in shp_files) {
      
      # Extract the plot number from the shapefile name (e.g., "Plot_12.shp" -> 12)
      plot_num <- as.numeric(str_extract(basename(shp_file), "\\d+"))
      
      # Routing Logic to prevent double-extraction in overlapping areas
      should_process <- TRUE
      
      if (!is.na(plot_num)) {
        if (is_top && plot_num > 21) {
          should_process <- FALSE
        } else if (is_bottom && plot_num < 22) {
          should_process <- FALSE
        }
      }
      
      # Skip to the next shapefile if it doesn't belong to this TIF
      if (!should_process) {
        next
      }
      
      print(paste("    -> Extracting polygons from:", basename(shp_file)))
      
      crowns <- st_read(shp_file, quiet = TRUE)
      
      # Match CRS
      if (st_crs(crowns) != st_crs(ortho)) {
        crowns <- st_transform(crowns, st_crs(ortho))
      }
      
      # Extract Mean values
      ext_vals <- exact_extract(index_stack, crowns, fun = 'mean', progress = FALSE)
      colnames(ext_vals) <- c("Mean_VARI", "Mean_GLI", "Mean_TGI")
      
      # Bind data, add metadata, and drop empty rows
      crowns_data <- crowns %>%
        st_drop_geometry() %>%
        bind_cols(ext_vals) %>%
        mutate(
          Date = date_str,
          Source_Raster = basename(tif_file),
          Plot_Number = plot_num  # Saving the parsed plot number for easy merging later
        ) %>%
        drop_na(Mean_VARI) 
      
      all_extractions[[length(all_extractions) + 1]] <- crowns_data
    }
  }
}

# --- FINALIZE AND EXPORT ---
print("========================================")
print("All folders processed. Merging datasets...")

final_master_df <- bind_rows(all_extractions)

output_path <- file.path(base_dir, "Master_RGB_Indices_All_Dates.csv")
write.csv(final_master_df, output_path, row.names = FALSE)

print(paste("Success! File saved to:", output_path))