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
library(exactextractr) # Added exactextractr

# Read in already processed Canopy height models ####

# Change this single variable for each new batch!
date_folder <- "03. 30 October 2025"

# Read in tree location for height extraction
# Make sure its CRS are 2048 please. Otherwise convert in QGIS with CRS layer function and then save as function.
trees <- st_read(paste0("E:/Remote Sensing Media/",date_folder,"/08. Crown shape file/All_Plots.shp"))

# Safely remove old columns if they exist to prevent duplicate name/ESRI errors
trees <- trees %>%
  select(-any_of(c("grop_ld","group_ulid", "N_GM", "id", "N_FG", "N_BG", "BBox")))

# Read in Canopy height models
ctg_chm <- rast(paste0("E:/Remote Sensing Media/",date_folder,"/07. Canopy Height Models/rasterize_canopy.vrt"))

# Plot Canopy Height Models and tree shape files
plot(ctg_chm, range=c(-0.5,4))
plot(trees$geometry, add = TRUE, col = "red")

# Extract tree heights ####
# Test

tic()
# Calculate metrics using exact_extract (Outputs directly as a vector)
trees$Tree_Height <- exact_extract(ctg_chm, trees, 'max')

# Save to shapefile (completely overwriting old files to prevent schema errors)
st_write(trees, paste0("E:/Remote Sensing Media/",date_folder,"/09. Tree heights/All Plots.shp"), delete_dsn = TRUE)

# Save lightweight CSV without the messy spatial geometry text
write.csv(st_drop_geometry(trees), paste0("E:/Remote Sensing Media/",date_folder,"/09. Tree heights/All Plots.csv"), row.names = FALSE)
toc()