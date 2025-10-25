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

#Plot number
Number <- 2

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
plot(las)

# Classified las
las_check(las_classified)
gnd <- filter_ground(las_classified)
plot(gnd, size = 3, bg = "white")

plot(las_normalised)

plot(las_chm)
plot(PlotTrees, add = TRUE, col = "red")














