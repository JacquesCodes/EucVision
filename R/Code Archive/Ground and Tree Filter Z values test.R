# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: Z-THRESHOLD SENSITIVITY ANALYSIS
# ──────────────────────────────────────────────────────────────────────────────
library(lidR)
library(dplyr)
library(stringr)

# 1. Configuration
base_dir <- "E:/Remote Sensing Media"
voxel_resolution <- 0.1 

# Pick ONE representative dataset to test your thresholds on
test_folder_name <- "20. 23 March 2026"
folder_path <- file.path(base_dir, test_folder_name)
normalised_dir <- file.path(folder_path, "06. Point Clouds Normalised")

# The thresholds you want to test
thresholds_to_test <- c(0.2, 0.3, 0.4, 0.5)

las_files <- list.files(normalised_dir, pattern = "\\.(las|laz)$", full.names = TRUE)

if (length(las_files) == 0) {
  stop("No normalized point clouds found in the target directory.")
}

plot_results <- list()

print("Starting Sensitivity Analysis...")

for (file in las_files) {
  plot_id_match <- str_extract(basename(file), "Plot_\\d+")
  plot_id <- ifelse(is.na(plot_id_match), basename(file), plot_id_match)
  
  las <- readLAS(file)
  if (is.empty(las)) next
  
  # Loop through each threshold for this specific plot
  for (t_val in thresholds_to_test) {
    
    # --- Ground Density ---
    ground_las <- filter_poi(las, Z <= t_val)
    if (!is.empty(ground_las)) {
      ground_density_m2 <- npoints(ground_las) / area(ground_las)
    } else {
      ground_density_m2 <- NA
    }
    
    # --- Canopy Density ---
    canopy_las <- filter_poi(las, Z > t_val)
    if (!is.empty(canopy_las)) {
      voxels <- voxel_metrics(canopy_las, ~length(Z), res = voxel_resolution)
      single_voxel_vol <- voxel_resolution^3
      canopy_density_m3 <- sum(voxels$V1) / (nrow(voxels) * single_voxel_vol)
    } else {
      canopy_density_m3 <- NA
    }
    
    # Store results
    plot_results[[length(plot_results) + 1]] <- data.frame(
      Plot_ID = plot_id,
      Threshold_m = t_val,
      Ground_Density_pts_m2 = round(ground_density_m2, 2),
      Canopy_Density_pts_m3 = round(canopy_density_m3, 2)
    )
  }
}

# Combine and save
sensitivity_df <- bind_rows(plot_results)

# Sort it nicely so you can read it easily: Group by Plot, then by Threshold
sensitivity_df <- sensitivity_df %>% arrange(Plot_ID, Threshold_m)

out_csv <- file.path(folder_path, "Threshold_Sensitivity_Test.csv")
write.csv(sensitivity_df, out_csv, row.names = FALSE)

print(paste("Done! Review the results at:", out_csv))