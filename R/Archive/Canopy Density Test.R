# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: CANOPY VS GROUND POINT DENSITY EXTRACTION
# ──────────────────────────────────────────────────────────────────────────────
# Description: This standalone script reads the normalized plot-level point clouds
#              from step 6 of the main pipeline. It splits the data at a 20cm 
#              Z-threshold to calculate 2D ground density (points/m²) and 
#              3D voxel canopy density (points/m³).
# ──────────────────────────────────────────────────────────────────────────────

library(lidR)
library(dplyr)
library(stringr)

# ──────────────────────────────────────────────────────────────────────────────
# 1. Configuration
# ──────────────────────────────────────────────────────────────────────────────
base_dir <- "E:/Remote Sensing Media"
voxel_resolution <- 0.1 # 10cm x 10cm x 10cm voxels for the canopy

# Set to a specific folder name to test one dataset (e.g., "20. 23 March 2026 0.6cm")
# Set to NULL to run the full batch process.
target_date_override <- "30. 30 June 2026 (ALS)"

exclude_list <- c("000. Projects",
                  "00. Baseline DTM",
                  "00. Dataset Template", 
                  "01. 25 February 2025",
                  "07. December 2025 (TLS)",
                  "17. 03 March 2026 (Multispectral)",
                  "20. 24 March 2026 (Multispectral)")

folders <- list.dirs(base_dir, recursive = FALSE)
dataset_folders <- folders[grepl("^\\d{2}\\.", basename(folders)) & !basename(folders) %in% exclude_list]

if (!is.null(target_date_override)) {
  dataset_folders <- dataset_folders[basename(dataset_folders) == target_date_override]
}

# ──────────────────────────────────────────────────────────────────────────────
# 2. Batch Processing Loop
# ──────────────────────────────────────────────────────────────────────────────
for (folder_path in dataset_folders) {
  date_folder <- basename(folder_path)
  file_date <- sub("^\\d+\\.\\s*", "", date_folder)
  file_date_safe <- gsub(" ", "_", file_date)
  
  print(paste("================================================================"))
  print(paste("EXTRACTING DENSITY METRICS:", date_folder))
  print(paste("================================================================"))
  
  # Pointing to the outputs from your main pipeline
  normalised_dir <- file.path(folder_path, "06. Point Clouds Normalised")
  out_dir <- file.path(folder_path, "10. Density Metrics")
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  out_csv_path <- file.path(out_dir, paste0("Plot_Densities_", file_date_safe, ".csv"))
  
  # Skip if already processed
  if (file.exists(out_csv_path)) {
    print("-> SKIPPED: Density metrics already exist for this date.")
    next
  }
  
  # Grab all the individual plot files
  las_files <- list.files(normalised_dir, pattern = "\\.(las|laz)$", full.names = TRUE)
  
  if (length(las_files) == 0) {
    print(paste("-> SKIPPED: No normalized point clouds found in", normalised_dir))
    next
  }
  
  plot_results <- list()
  
  for (file in las_files) {
    # Extract plot ID from your naming convention (e.g., "Plot_12_23_March.las")
    plot_id_match <- str_extract(basename(file), "Plot_\\d+")
    plot_id <- ifelse(is.na(plot_id_match), basename(file), plot_id_match)
    
    # Load the specific plot
    las <- readLAS(file)
    
    # If the file is empty or corrupted, skip
    if (is.empty(las)) next
    
    # ----------------------------------------------------------------------
    # A. Ground Density (2D: points / m²)
    # ----------------------------------------------------------------------
    # Filter for points at or below 30cm
    ground_las <- filter_poi(las, Z <= 0.3)
    
    if (!is.empty(ground_las)) {
      # Calculate the 2D area of the plot footprint
      plot_area_m2 <- area(ground_las)
      ground_pts <- npoints(ground_las)
      ground_density_m2 <- ground_pts / plot_area_m2
    } else {
      ground_density_m2 <- NA
    }
    
    # ----------------------------------------------------------------------
    # B. Canopy Density (3D: points / m³)
    # ----------------------------------------------------------------------
    # Filter for points strictly above 30cm
    canopy_las <- filter_poi(las, Z > 0.3)
    
    if (!is.empty(canopy_las)) {
      # Voxelize the canopy space. This returns a data frame where 'V1' 
      # is the number of points inside each 10cm x 10cm x 10cm cube.
      voxels <- voxel_metrics(canopy_las, ~length(Z), res = voxel_resolution)
      
      # Calculate volume of a single voxel (0.1 * 0.1 * 0.1 = 0.001 m³)
      single_voxel_vol <- voxel_resolution^3
      
      # Total points in the canopy divided by the total physical volume 
      # those canopy points occupy.
      total_canopy_pts <- sum(voxels$V1)
      total_canopy_volume_m3 <- nrow(voxels) * single_voxel_vol
      
      canopy_density_m3 <- total_canopy_pts / total_canopy_volume_m3
    } else {
      canopy_density_m3 <- NA
    }
    
    # Store the iteration results
    plot_results[[length(plot_results) + 1]] <- data.frame(
      Flight_Date = file_date,
      Plot_ID = plot_id,
      Ground_Density_pts_m2 = round(ground_density_m2, 2),
      Canopy_Density_pts_m3 = round(canopy_density_m3, 2)
    )
  }
  
  # Combine all plots for this flight date into one dataframe
  if (length(plot_results) > 0) {
    flight_df <- bind_rows(plot_results)
    write.csv(flight_df, out_csv_path, row.names = FALSE)
    print(paste("-> Successfully saved:", basename(out_csv_path)))
  }
  
  # Memory management
  rm(las, ground_las, canopy_las, plot_results, flight_df)
  gc()
}

print("================================================================")
print("DENSITY EXTRACTION COMPLETE!")
print("================================================================")