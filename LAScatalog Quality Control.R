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


#Plot number
Number <- 17

las <- readLAS(paste0("E:/Remote Sensing Media/0. R Projects/Point Cloud/1. Clipped/Plot ",Number,".las"))
las_classified <- readLAS(paste0("E:/Remote Sensing Media/0. R Projects/Point Cloud/2. Ground Classified/Plot ",Number, "_classified.las"))
las_normalised <- readLAS(paste0("E:/Remote Sensing Media/0. R Projects/Point Cloud/3. Normalised/Plot ",Number, "_classified_normalised.las"))
las_chm <- rast(paste0("E:/Remote Sensing Media/0. R Projects/Point Cloud/4. Canopy Height Model/Plot ",Number, "_classified_normalised_chm.tif"))

trees <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/QGIS Combined Output/All_Plots.shp")

PlotTrees <- trees[trees$Plot == paste0("Plot ",Number),]

# Ensure both have an ID column
PlotTrees$ID <- 1:nrow(PlotTrees)
# Calculate metrics
tree_heights <- terra::extract(las_chm, PlotTrees, fun = max, na.rm = TRUE)
# Join results back using the ID
trees_with_heights <- left_join(PlotTrees, st_drop_geometry(tree_heights), by = "ID")

# Cropped las
plot(las, size = 4, bg = "#F1F8F8")

# Classified las
las_check(las_classified)
gnd <- filter_ground(las_classified)
plot(gnd, size = 3, bg = "#F1F8F8")

dtm_tin_0 <- rasterize_terrain(las_classified, res = 1, algorithm = tin())
plot_dtm3d(dtm_tin_0, bg = "#F1F8F8") 

plot(las_normalised,bg = "white")

plot(las_chm, col = height.colors(50))
plot(PlotTrees, add = TRUE, col = "red")



# 8.1 Individual Tree Detection (ITD)

MinimumTreeHeight <- 0.6

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
x <- plot(las_normalised, bg = "white", size = 4)
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

