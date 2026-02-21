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

# Change this single variable for each new batch!
date_folder <- "13. 29 January 2026"

# My path to the remote sensing dataset
myPath <- paste0("E:/Remote Sensing Media/",date_folder,"/")

#Plot number
Number <- 17

las <- readLAS(paste0(myPath,"04. Point clouds clipped/Plot_",Number,".las"))
las_classified <- readLAS(paste0(myPath,"05. Point clouds ground classified/Plot_",Number, "_classified.las"))
las_normalised <- readLAS(paste0(myPath,"06. Point clouds normalised/Plot_",Number, "_classified_normalised.las"))
las_chm <- rast(paste0(myPath,"07. Canopy Height Models/Plot_",Number, "_classified_normalised_chm.tif"))

trees <- st_read(paste0(myPath,"08. Crown shape file/All_Plots.shp"))

PlotTrees <- trees[trees$Plot == paste0("Plot_",Number),]

# # Ensure both have an ID column
# PlotTrees$ID <- 1:nrow(PlotTrees)
# # Calculate metrics
# tree_heights <- terra::extract(las_chm, PlotTrees, fun = max, na.rm = TRUE)
# # Join results back using the ID
# trees_with_heights <- left_join(PlotTrees, st_drop_geometry(tree_heights), by = "ID")

# Plot cropped las
# plot(las, size = 4, bg = "white")
plot(las, size = 2, color = "RGB", bg = "white")

# Plot classified las
gnd <- filter_ground(las_classified)
plot(gnd, size = 4, bg = "white")

# Plot DTM
dtm_tin_0 <- rasterize_terrain(las_classified, res = 1, algorithm = tin())
plot_dtm3d(dtm_tin_0, bg = "white") 

# Plot normalised las
plot(las_normalised,bg = "white")

# Plot canopy height model
plot(las_chm, col = height.colors(50))
plot(PlotTrees, add = TRUE, col = "red")


# 8.1 Individual Tree Detection (ITD)

MinimumTreeHeight <- 0.5

# Local Maximum Filter with variable windows size
# f <- function(z) {1 * z + 0.5}
# lmf_algorithm <- lmf(ws = f, hmin = MinimumTreeHeight, shape = "circular")

# create Local Maximum Filter (lmf) function for the "ws" search
lmf_algorithm <- lmf(ws = 3, hmin = MinimumTreeHeight, shape = "circular")

# Locate trees in a circle with a diameter of "ws" in meters
ttops <- locate_trees(las = las_chm, algorithm = lmf_algorithm)

# Tree detection results in 2D
plot(las_chm, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)

# Tree detection results can also be visualized in 3D!
x <- plot(las_normalised, bg = "white", size = 2)
add_treetops3d(x, ttops, radius = 0.15, fastTransparency = TRUE, alpha = 0.8)




# Video

# --- Define spin motion ---
spin <- spin3d(axis = c(0, 0, 1), rpm = 6)

# Doesn't work if you make the plot fullscreen.
# --- Save the animation ---
movie3d(
  movie = "Tree_Tops_animation",   # base name for output
  dir = getwd(),                   # output directory
  spin,                            # the animation function
  duration = 10,                   # 10 seconds
  fps = 25,                        # frames per second
  clean = TRUE,                    # remove frames after combining
  type = "gif"                     # try making a .gif directly
)


# Testing csf cloth resolution = 1
tic()
las_classified_2 <- classify_ground(las, csf(sloop_smooth = TRUE, 
                                             class_threshold = 0.05, 
                                             cloth_resolution = 0.5,
                                             rigidness = 1,
                                             time_step = 0.65))

plot(las_classified_2, color = "Classification", size = 3, bg = "white")
print("CSF Ground filtering time:")
toc()

# Classified las
gnd <- filter_ground(las_classified_2)
plot(gnd, size = 3, bg = "white")

# Classified ground
dtm_tin_0 <- rasterize_terrain(las_classified_2, res = 0.01, algorithm = tin())
plot_dtm3d(dtm_tin_0, bg = "white") 

# Normalised 
norm <- normalize_height(las = las_classified_2, algorithm = tin())
plot(norm, size = 1, color = "RGB", bg = "white")