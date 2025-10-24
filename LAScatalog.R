library(lidR)
library(RCSF)
library(RMCC)
library(gstat)
library(sf)
library(tictoc)
library(geometry)
library(dplyr)
library(future)

################################################################################
# DO NOT USE LAZ.!! ONLY USE LAS. IT IMPROVES PERFORMANCE X10!
################################################################################

# Read in headers of Las. files in a folder called a catalog (ctg)
ctg <- readLAScatalog("E:/Remote Sensing Media/0. R Projects/Point Cloud/1. Clipped/")
plot(ctg)

las_check(ctg)

trees <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/QGIS Combined Output/All_Plots.shp")

# View trees on plots:
plot(ctg)
plot(trees, add = TRUE, col = "red")

# Multiple threads mode
plan(multisession)

# Each plot is independent and buffers are not needed
opt_independent_files(ctg) <- TRUE

tic()
# Write to disk rather than memory:
opt_output_files(ctg) <- paste0("E:/Remote Sensing Media/0. R Projects/Point Cloud/2. Ground Classified/", "{*}_classified")
# Ground classifications:
mycsf <- csf(sloop_smooth = TRUE, class_threshold = 1, cloth_resolution = 1, time_step = 1)
ctg_classified <- classify_ground(ctg, mycsf)
toc()

tic()
# Write to disk rather than memory:
opt_output_files(ctg_classified) <- paste0("E:/Remote Sensing Media/0. R Projects/Point Cloud/3. Normalised/", "{*}_normalised")
# A point cloud-based normalization without a raster:
ctg_normalised <- normalize_height(ctg_classified, tin())
toc()

tic()
# Write to disk rather than memory:
opt_output_files(ctg_normalised) <- paste0("E:/Remote Sensing Media/0. R Projects/Point Cloud/4. Canopy Height Model/", "{*}_chm")
# Rasterize canopy with interpolation:
ctg_chm <- rasterize_canopy(ctg_normalised, res = 0.01, algorithm = p2r(na.fill = tin()))
print("Rasterize canopy time:")
toc()

tic()
# Ensure both have an ID column
trees$ID <- 1:nrow(trees)

# Calculate metrics
tree_heights <- terra::extract(ctg_chm, trees, fun = max, na.rm = TRUE)

# Join results back using the ID
trees_with_heights <- left_join(trees, st_drop_geometry(tree_heights), by = "ID")

# Save to file
st_write(trees_with_heights, paste0("E:/Remote Sensing Media/0. R Projects/Point Cloud/Heights/All Plots.shp"), delete_dsn = TRUE)

print("Finish Time:")
toc()








