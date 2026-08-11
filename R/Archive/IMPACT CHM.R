library(terra)
library(ggplot2)
library(ggspatial)
library(tidyterra)

# 1. Define the exact path to the IMPACT Site CHM
path_chm <- "E:/Remote Sensing Media/23. 13 April 2026/07. Canopy Height Models IMPACT/IMPACT_Site_CHM_Single_13_April_2026.tif"

# --- Safely Read in the Data ---
if (file.exists(path_chm)) {
  las_chm <- rast(path_chm)
  
  # Filter out the rogue photogrammetry noise
  las_chm[las_chm > 7] <- NA  # Replace 15 with your realistic upper threshold
  las_chm[las_chm < 0] <- 0    # Cleans up any sub-surface negative noise
  
  # --- Smooth the CHM ---
  # Define the moving window 'w' (5x5 matrix based on your original code)
  w <- matrix(1, nrow = 5, ncol = 5) 
  
  # Apply the focal smoothing function to the filtered CHM
  smoothed_chm <- terra::focal(las_chm, w = w, fun = max, na.rm = TRUE)
  
  message("Loaded and smoothed: IMPACT_Site_CHM_Single_13_April_2026.tif")
  
} else {
  message("Error: Could not find the specified CHM file. Please check the file path.")
  smoothed_chm <- NULL 
}

# --- Plotting Section ---

# Only generate the map if the smoothed CHM was successfully created
if (!is.null(smoothed_chm)) {
  chm_map <- ggplot() +
    # 1. Add the Smoothed Canopy Height Model raster layer
    geom_spatraster(data = smoothed_chm) +
    
    # 2. Define the color palette for the raster
    scale_fill_viridis_c(
      name = "Canopy\nHeight (m)", 
      option = "viridis", 
      na.value = "transparent"
    ) +
    
    # 3. Add the Scale Bar
    annotation_scale(
      location = "br", 
      width_hint = 0.4, 
      text_cex = 0.8
    ) +
    
    # 4. Add the North Arrow 
    annotation_north_arrow(
      location = "tr", 
      which_north = "false",
      style = north_arrow_fancy_orienteering()
    ) +
    
    # 5. Formatting and Titles (Updated for the IMPACT site)
    labs(
      title = "Smoothed Canopy Height Model",
      subtitle = "IMPACT Site - 13 April 2026",
      x = "Longitude",
      y = "Latitude"
    ) +
    
    # 6. Coordinate system
    coord_sf() +
    
    # 7. Theme adjustments
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "in")
    )
  
  # View the map
  print(chm_map)
}