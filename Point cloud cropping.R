library(lidR)
library(RCSF)
library(RMCC)
library(gstat)
library(sf)
library(tictoc)
library(geometry)
library(dplyr)
library(future)
library(rgl)
library(magick)

################################################################################
# DO NOT USE LAZ.!! ONLY USE LAS. IT IMPROVES PERFORMANCE X10!
################################################################################



# 1. Link to .las file ####

# Link October .las file
# Top
Link <- "E:/Remote Sensing Media/7. 30 October 2025/Point cloud/Lourensford_Top_30 October_group1_densified_point_cloud.las"
# Bottom
Link <- "E:/Remote Sensing Media/7. 30 October 2025/Point cloud/Lourensford_30 Ocotober_Bottom_group1_densified_point_cloud.las"

# Link to bigger .las file
Link <- "E:/Remote Sensing Media/6. September 2025/Point Cloud/SU Lourensford September 2025_point cloud-001.las"

# Link to smaller 50cm/pixel .las file
Link <- "E:/Remote Sensing Media/6. September 2025/Point Cloud/SU Lourensford September 2025_point cloud_50cm.las"

# Link to the March IAS .las
Link <- "E:/Remote Sensing Media/4. March 2025/DJI Matrice 3TD/RGB/15m/M3E 15mAGL_pointcloud.laz"

Link <- "E:/Remote Sensing Media/4. March 2025/DJI Matrice 3TD/Thermal/39m/M3TD 39mAGL_pointcloud.las"


Link <- "E:/Remote Sensing Media/0. R Projects/2. 30 October 2025/1. Clipped/Plot 18.las"
  


# Link to Clipped .las file
Plot <- 18
Folder <- "1. Clipped"
Link <- paste0("E:/Remote Sensing Media/0. R Projects/Point Cloud/",Folder,"/Plot ",Plot,".las")


Link <- "E:/Remote Sensing Media/0. R Projects/Point Cloud/3. Normalised/Plot 37_classified_normalised.las"

# You can filter attributes out if needed to decrease RAM usage by over 50%!
tic()
las <- readLAS(Link, select = "xyzRGB")
print("Las. read in time:")
toc()

plot(las, bg = "white", color = "RGB", size = 3)

# Create a GIF in the current working directory
# Spin around the Z-axis for 10 seconds
spin <- spin3d(axis = c(0, 0, 1), rpm = 6)
movie3d(spin,dir = getwd(),duration = 10, movie = "my_rgl_animation.gif", type = "gif", clean = TRUE)


# 2. Read and check point cloud ####

# Read .las file into memory
las <- readLAS(Link)

# Check .las details
print(las)
las_check(las)

# 3. How to clip a point cloud ####

# Region of interest (ROI)
# las <- clip_circle(las, x = -6985, y = -3766340, radius = 5)
# plot(las, bg = "white", size = 4)

# Find extent through plot boundaries

skip_nums <- c(24, 58, 75, 36, 49, 67)

# NOTE: Check where i starts to count! I changed it!
for (i in setdiff(15:75, skip_nums)) {
tic()  
  PlotNumber <- i
  shape_file <- st_read(paste0("E:/Remote Sensing Media/1. QGIS Projects/Michelle/Michelle QGIS Project/1 September 2025/Plot ",PlotNumber,".shp"))
  extent <- ext(shape_file)
  las_clipped <- clip_rectangle(las, xleft = extent[1], xright = extent[2], ybottom = extent[3], ytop = extent[4])
  writeLAS(las_clipped, file.path("E:/Remote Sensing Media/0. R Projects/Point Cloud/Clipped", paste0("Plot ",PlotNumber,".las")), index = FALSE)
  
  print(paste0("Cropping done for plot ",PlotNumber))
toc()
}

tic()  
# October Crop
for (i in 22:75) {
# for (i in 1:21) {

  PlotNumber <- i
  shape_file <- st_read(paste0("E:/Remote Sensing Media/1. QGIS Projects/1. Plot boundaries for cropping/Plots shape files/id_",PlotNumber,".shp"))
  extent <- ext(shape_file)
  las_clipped <- clip_rectangle(las, xleft = extent[1], xright = extent[2], ybottom = extent[3], ytop = extent[4])
  writeLAS(las_clipped, file.path("E:/Remote Sensing Media/0. R Projects/2. 30 October 2025/1. Clipped/", paste0("Plot ",PlotNumber,".las")), index = FALSE)
  
  print(paste0("Cropping done for plot ",PlotNumber))

}
toc()




las <- las_clipped

# 4. Ground classification ####


# Progressive Morphological Filter (PMF) 
tic()
ws <- seq(3, 12, 3)
th <- seq(0.1, 1.5, length.out = length(ws))
las1 <- classify_ground(las, algorithm = pmf(ws = ws, th = th))
plot(las1, color = "Classification", size = 3, bg = "lightblue")
print("PMF Ground filtering time:")
toc()

# Always keep cloth_resolution 0.5
# Keep class threshold below 0.01 otherwise it starts to exclude short trees

tic()
####### Winner Winner Chicken Dinner! ####
las4 <- classify_ground(las, csf(sloop_smooth = TRUE, class_threshold = 0.01, cloth_resolution = 0.5, time_step = 1))
plot(las4, color = "Classification", size = 3, bg = "lightblue")
print("CSF Ground filtering time:")
toc()




# Multiscale Curvature Classification (MCC)
# Preferred but takes longer
# las1 <- classify_ground(las, mcc(1.5,0.3))
# plot(las1, color = "Classification", size = 3, bg = "lightblue") 


# Display ground points
# gnd <- filter_ground(las1)
# plot(gnd, size = 3, bg = "white")


# 5. Digital terrain model ####

# # Preferred DTM model with the Triangular irregular network (TIN) algorithm
# dtm_tin <- rasterize_terrain(las1, res = 0.1, algorithm = tin())
# plot_dtm3d(dtm_tin, bg = "white") 
# 
# # DTM model with the Kriging algorithm
# dtm_kriging_1 <- rasterize_terrain(las1, algorithm = kriging(k = 40))
# plot_dtm3d(dtm_kriging_1, bg = "white") 


# 6. Height normalization ####

# Normalised 
tic()
nlas <- normalize_height(las4, tin())
plot(nlas, size = 4, bg = "white")
# writeLAS(nlas, file.path("E:/Remote Sensing Media/0. R Projects/Point Cloud/Normalised Plots/", paste0("Plot ",PlotNumber,".las")), index = FALSE)
print("Height normalisation time:")
toc()

plot(nlas, size= 2, color = "RGB", bg = "white")

# 7. Digital Surface Model (DSM) and Canopy Height model (CHM) ####

# Digital Surface Model (DSM) and Canopy Height model (CHM)
# chm <- rasterize_canopy(nlas, algorithm = p2r())
# col <- height.colors(25)
# plot(chm, col = col)

tic()
# Rasterize canopy with interpolation
chm <- rasterize_canopy(nlas, res = 0.01, algorithm = p2r(na.fill = tin()))
# writeRaster(chm,"E:/Remote Sensing Media/0. R Projects/Point Cloud/Canopy Height Model/", paste0("Plot ",PlotNumber,".las"), index = FALSE, overwrite = TRUE)
print("Rastierize canopy time:")
toc()


# smoothed <- terra::focal(chm, w, fun = mean, na.rm = TRUE)
# plot(smoothed, col = col)

# 7.4 Post-processing a CHM ####

# subcircle ONLY NEEDED FOR LOW POINT DENSITY!

# fill.na <- function(x, i=5) { if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }}
# w <- matrix(1, 3, 3)
# 
# chm <- rasterize_canopy(nlas, res = 0.01, algorithm = p2r(subcircle = 0.15), pkg = "terra")
# filled <- terra::focal(chm, w, fun = fill.na)
# smoothed <- terra::focal(chm, w, fun = mean, na.rm = TRUE)
# 
# chms <- c(chm, filled, smoothed)
# names(chms) <- c("Base", "Filled", "Smoothed")
# plot(chms, col = col)


# 8. Individual tree detection and segmentation ####
# 
# # 8.1 Individual Tree Detection (ITD)
# 
# MinimumTreeHeight <- 0.5
# 
# # Local Maximum Filter with variable windows size
# # f <- function(z) {1 * z + 0.5}
# # lmf_algorithm <- lmf(ws = f, hmin = MinimumTreeHeight, shape = "circular")
# 
# # create Local Maximum Filter (lmf) function for the "ws" search
# lmf_algorithm <- lmf(ws = 1, hmin = MinimumTreeHeight, shape = "circular")
# 
# # Locate trees in a circle with a diameter of "ws" in meters
# ttops <- locate_trees(las = nlas[nlas$Z>= 0], algorithm = lmf_algorithm)
# 
# # Tree detection results in 2D
# plot(chm, col = height.colors(50))
# plot(sf::st_geometry(ttops), add = TRUE, pch = 3)
# 
# # Tree detection results can also be visualized in 3D!
# x <- plot(nlas, bg = "white", size = 4)
# add_treetops3d(x, ttops, radius = 0.1, fastTransparency = TRUE, alpha = 1)
# 
# # Individual Tree Segmentation (ITS)
# 
# algo <- dalponte2016(chm, ttops)
# las_segmented <- segment_trees(nlas, algo) # segment point cloud
# plot(las_segmented, bg = "white", size = 4, color = "treeID") # visualize trees

# 8.1 Individual Tree Detection (ITD)

MinimumTreeHeight <- 0.8

# Local Maximum Filter with variable windows size
# f <- function(z) {1 * z + 0.5}
# lmf_algorithm <- lmf(ws = f, hmin = MinimumTreeHeight, shape = "circular")

# create Local Maximum Filter (lmf) function for the "ws" search
lmf_algorithm <- lmf(ws = 2, hmin = MinimumTreeHeight, shape = "circular")

# Locate trees in a circle with a diameter of "ws" in meters
ttops <- locate_trees(las = chm, algorithm = lmf_algorithm)

# Tree detection results in 2D
plot(chm, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)




# Tree detection results can also be visualized in 3D!
x <- plot(nlas, bg = "white", color = "RGB", size = 2)
add_treetops3d(x, ttops, radius = 0.15, fastTransparency = TRUE, alpha = 0.8)

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

# 9. Tree metrics ####

# trees <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/QGIS Extracted data/1 September 2025/Plot 37.shp")

# Read in shape file
trees <- st_read(paste0("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/QGIS Extracted data/1 September 2025/Plot ",Plot,".shp"))

# Ensure both have an ID column
trees$ID <- 1:nrow(trees)

# Calculate metrics
tree_heights <- terra::extract(chm, trees, fun = max, na.rm = TRUE)

# Join results back using the ID
trees_with_heights <- left_join(trees, st_drop_geometry(tree_heights), by = "ID")

# Save to file
st_write(trees_with_heights, paste0("E:/Remote Sensing Media/0. R Projects/Point Cloud/Heights/Plot ",PlotNumber,".shp"), delete_dsn = TRUE)

print("Finish Time:")
toc()











# Test zone ####
# Testing different DSM/CHM algorithms and settings

# 1.162803 cm / point on average
dsm <- rasterize_canopy(las, res = 0.01, algorithm = p2r(na.fill = tin()))
smoothed <- terra::focal(dsm, w, fun = mean, na.rm = TRUE)
plot(smoothed, col = col)


# Read in shape file
trees <- st_read(paste0("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/QGIS Extracted data/1 September 2025/Plot ",Plot,".shp"))

# Ensure both have an ID column
trees$ID <- 1:nrow(trees)

# Calculate metrics
tree_heights <- terra::extract(dsm, trees, fun = max, na.rm = TRUE)

# Join results back using the ID
trees_with_heights <- left_join(trees, st_drop_geometry(tree_heights), by = "ID")

mean(trees_with_heights$Z)
