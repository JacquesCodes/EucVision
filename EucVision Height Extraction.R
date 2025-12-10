library(lidR)
library(RCSF)
library(RMCC)
library(gstat)
library(sf)
library(tictoc)
library(geometry)
library(dplyr)
library(future)
library(sp)
library(terra)


# Read in already processed Canopy height models ####

# My path to the remote sensing dataset
myPath <- "E:/Remote Sensing Media/09. 14 November 2025/"

# Read in tree location for height extraction
# Make sure its CRS are 2048 please. Otherwise convert in QGIS with layer CRS function and then save as function.
trees <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/QGIS Combined Output/30 October 2025/All_Plots_Reprojected.shp")
# trees <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/QGIS Combined Output/1 September 2025/All_Plots.shp")

# Read in Canopy height models
ctg_chm <-rast(paste0(myPath,"07. Canopy Height Models/rasterize_canopy.vrt"))

# Plot Canopy Height Models and tree shape files
plot(ctg_chm,range=c(-0.5,4))
plot(trees$geometry, add = TRUE, col = "red")

# Extract tree heights ####

tic()
# Ensure both have an ID column
trees$ID <- 1:nrow(trees)

# Calculate metrics
tree_heights <- terra::extract(ctg_chm, trees, fun = max, na.rm = TRUE)

# Join results back using the ID
trees_with_heights <- left_join(trees, st_drop_geometry(tree_heights), by = "ID")

# Save to shape and excel file
st_write(trees_with_heights, paste0(myPath,"08. Tree heights/All Plots.shp"))
st_write(trees_with_heights, paste0(myPath,"08. Tree heights/All Plots.csv"))
toc()