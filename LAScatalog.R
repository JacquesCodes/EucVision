library(lidR)
library(RCSF)
library(RMCC)
library(gstat)
library(sf)
library(tictoc)
library(geometry)

################################################################################
# DO NOT USE LAZ.!! ONLY USE LAS. IT IMPROVES PERFORMANCE X10!
################################################################################

# Read in headers of Las. files in a folder
ctg <- readLAScatalog("E:/Remote Sensing Media/0. R Projects/Point Cloud/Clipped/")
plot(ctg)

las_check(ctg)

trees <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/QGIS Combined Output/All_Plots.shp")

# View trees on plots:
plot(ctg)
plot(trees, add = TRUE, col = "red")

plan(multisession)

opt_independent_files(ctg) <- TRUE
################################################################################
# No ground points found. Impossible to compute a DTM.
# Write code to batch process DTM
################################################################################

# Write to disk rather than memory:
opt_output_files(ctg) <- paste0("E:/Remote Sensing Media/0. R Projects/Point Cloud/LAScatalog", "{*}_normalised")
# A point cloud-based normalization without a raster is also possible:
ctg_norm <- normalize_height(ctg, tin())

