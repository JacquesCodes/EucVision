library(lidR)
library(terra)
library(ggplot2)
library(ggspatial)
library(tidyterra)

# 1. Define variables based on your specific file path
date_folder <- "17. 02 March 2026"
myPath <- paste0("E:/Remote Sensing Media/", date_folder, "/")
Number <- 1

# Define the exact file path
las_file <- paste0(myPath, "06. Point clouds normalised/Plot_", Number, "_classified_normalised.las")

# 2. Read the LAS file 
# We include the filter to drop ground points, noise, and keep only first returns to speed up CHM generation
las_normalised <- readLAS(las_file, filter = "-drop_class 2 -drop_z_below 0 -keep_first")

# 3. Define the resolutions to test (in meters)
resolutions <- c(0.01, 0.02, 0.05, 0.10, 0.20)

# Create an empty list to store the plot objects if you want to access them later
chm_plots <- list()

# 4. Loop through each resolution
for (res in resolutions) {
  
  # Print progress to the console
  message(paste0("Generating CHM at ", res * 100, "cm resolution..."))
  
  # Generate the CHM
  # Using your previously selected point-to-raster algorithm with TIN interpolation
  chm <- rasterize_canopy(las_normalised, 
                          res = res, 
                          algorithm = p2r(subcircle = 0.015, na.fill = tin()))
  
  # Create the map using your ggplot template
  chm_map <- ggplot() +
    # Add the Canopy Height Model raster layer
    geom_spatraster(data = chm) +
    
    # Define the color palette for the raster
    scale_fill_viridis_c(
      name = "Canopy Height (m)", 
      option = "viridis", 
      na.value = "transparent"
    ) +
    
    # Add the Scale Bar
    annotation_scale(
      location = "br", 
      width_hint = 0.4, 
      text_cex = 0.8
    ) +
    
    # Add the North Arrow
    annotation_north_arrow(
      location = "tr", 
      which_north = "true", 
      style = north_arrow_fancy_orienteering()
    ) +
    
    # Formatting and Titles (Dynamically updating the subtitle with the resolution)
    labs(
      title = paste("Canopy Height Model - Plot", Number),
      subtitle = paste0("Extracted from DJI Matrice 3D SFM | Resolution: ", res * 100, "cm"), 
      x = "Longitude",
      y = "Latitude"
    ) +
    
    # Coordinate system
    coord_sf() +
    
    # Theme adjustments
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "in")
    )
  
  # Print the map to the viewer
  print(chm_map)
  
  # Store the plot in the list using the resolution as the name (e.g., "res_1cm")
  list_name <- paste0("res_", res * 100, "cm")
  chm_plots[[list_name]] <- chm_map
}

message("All resolutions generated and plotted successfully!")