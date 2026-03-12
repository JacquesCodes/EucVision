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
date_folder <- "17. 02 March 2026"

# My path to the remote sensing dataset
myPath <- paste0("E:/Remote Sensing Media/",date_folder,"/")

#Plot number
Number <- 38

las_normalised <- readLAS(paste0(myPath,"06. Point clouds normalised/Plot_",Number, "_classified_normalised.las"))
las_chm <- rast(paste0(myPath,"07. Canopy Height Models/Plot_",Number, "_classified_normalised_chm.tif"))

# Plot normalised las with axes turned on for scale reference
plot(las_normalised, color = "RGB", size = 2, bg = "white", axis = TRUE)

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
  
  # 4. Add the Scale Bar from ggspatial (bottom right)
  annotation_scale(
    location = "br", 
    width_hint = 0.4, 
    text_cex = 0.8
  ) +
  
  # 5. Add the North Arrow from ggspatial (Moved outside, top right)
  annotation_north_arrow(
    location = "tr", 
    which_north = "false",
    # Negative padding pushes the arrow OUTSIDE the main plot panel
    pad_x = unit(-0.6, "in"), 
    pad_y = unit(-0.6, "in"), 
    style = north_arrow_fancy_orienteering()
  ) +
  
  # 6. Formatting and Titles
  labs(
    title = paste("Canopy Height Model - Plot", Number),
    subtitle = "Extracted from DJI Matrice 3D SFM,",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  # 7. CRITICAL: Turn off clipping so the compass isn't hidden by the plot edges
  coord_sf(clip = "off") +
  
  # 8. Use a clean theme and expand margins for the compass
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
    # Expand top (t) and right (r) margins to create canvas space for the compass
    plot.margin = margin(t = 1, r = 1, b = 0.2, l = 0.2, unit = "in")
  )

# View the map
print(chm_map)