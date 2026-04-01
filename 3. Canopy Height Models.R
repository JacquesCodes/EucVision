library(lidR)
library(RCSF)
library(RMCC)
library(gstat)
library(sf)
library(tictoc)
library(geometry)
library(dplyr)
library(future)
library(terra)
library(rgl)
library(ggplot2)
library(ggspatial)
library(tidyterra)

# Change this single variable for each new batch!
date_folder <- "19. 16 March 2026"

# Extract the date part by removing the leading folder number, dot, and space
file_date <- sub("^\\d+\\.\\s*", "", date_folder)

# Replace spaces with underscores for safer file naming conventions
file_date_safe <- gsub(" ", "_", file_date)

# My path to the remote sensing dataset
myPath <- paste0("E:/Remote Sensing Media/",date_folder,"/")

# Plot number
Number <- 38

# Dynamically construct the single-date file names
name_clipped <- paste0("Plot_", Number, "_", file_date_safe)
name_classified <- paste0(name_clipped, "_classified")
name_normalised <- paste0(name_classified, "_normalised")
name_chm <- paste0(name_normalised, "_chm")

# Build the full file paths with Title Case directories
path_normalised <- paste0(myPath, "06. Point Clouds Normalised/", name_normalised, ".las")
path_chm <- paste0(myPath, "07. Canopy Height Models/", name_chm, ".tif")

# --- Safely Read in the Data ---
if (file.exists(path_normalised)) {
  las_normalised <- readLAS(path_normalised)
  message("Loaded: ", name_normalised, ".las")
} else {
  message("Skipped: Could not find folder or file for 06. Point Clouds Normalised")
  las_normalised <- NULL
}

if (file.exists(path_chm)) {
  las_chm <- rast(path_chm)
  
  # Filter out the rogue photogrammetry noise
  las_chm[las_chm > 15] <- NA  # Replace 8 with your realistic upper threshold
  las_chm[las_chm < 0] <- 0   # Optional: Cleans up any sub-surface negative noise
  
  message("Loaded: ", name_chm, ".tif")
} else {
  message("Skipped: Could not find folder or file for 07. Canopy Height Models")
  las_chm <- NULL
}

# --- Plotting Section ---

# Plot normalised las (Removed "axis = TRUE" to prevent the warning message!)
# if (!is.null(las_normalised)) {
#   plot(las_normalised, color = "RGB", size = 2, bg = "white")
#   plot(las_normalised, size = 2, bg = "white")
# }

# Only generate the map if the CHM loaded successfully
if (!is.null(las_chm)) {
  chm_map <- ggplot() +
    # 1. Add the Canopy Height Model raster layer
    geom_spatraster(data = las_chm) +
    
    # 2. Add the tree locations (uncomment if PlotTrees is loaded)
    # geom_sf(data = PlotTrees, color = "red", size = 1.5, shape = 16) +
    
    # 3. Define the color palette for the raster
    scale_fill_viridis_c(
      name = "Canopy Height (m)", 
      option = "viridis", 
      na.value = "transparent"
    ) +
    
    # 4. Add the Scale Bar
    annotation_scale(
      location = "br", 
      width_hint = 0.4, 
      text_cex = 0.8
    ) +
    
    # 5. Add the North Arrow (Moved INSIDE to prevent title collisions)
    annotation_north_arrow(
      location = "tr", 
      which_north = "false",
      style = north_arrow_fancy_orienteering()
    ) +
    
    # 6. Formatting and Titles
    labs(
      title = "Canopy Height Model",
      subtitle = paste0("Plot ",Number," - ",file_date),
      x = "Longitude",
      y = "Latitude"
    ) +
    
    # 7. Coordinate system
    coord_sf() +
    
    # 8. Theme adjustments
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
      # FIX: Angle the longitude text so it never overlaps
      axis.text.x = element_text(angle = 45, hjust = 1),
      # Normalized margins since the compass is no longer pushed outside
      plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "in")
    )
  
  # View the map
  print(chm_map)
}