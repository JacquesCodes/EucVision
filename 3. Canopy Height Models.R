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
date_folder <- "20. 16 March 2026"

# My path to the remote sensing dataset
myPath <- paste0("E:/Remote Sensing Media/",date_folder,"/")

#Plot number
Number <- 5

las_normalised <- readLAS(paste0(myPath,"06. Point clouds normalised/Plot_",Number, "_classified_normalised.las"))
las_chm <- rast(paste0(myPath,"07. Canopy Height Models/Plot_",Number, "_classified_normalised_chm.tif"))

# Plot normalised las with axes turned on for scale reference
plot(las_normalised, color = "RGB", size = 2, bg = "white", axis = TRUE)
plot(las_normalised, size = 2, bg = "white", axis = TRUE)

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
    which_north = "false", # Changed to "true" so it points correctly
    style = north_arrow_fancy_orienteering()
  ) +
  
  # 6. Formatting and Titles
  labs(
    title = paste("Canopy Height Model - Plot", Number),
    subtitle = "Extracted from DJI Matrice 3D SFM", # Removed trailing comma
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