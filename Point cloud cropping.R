library(lidR)
library(RCSF)
library(RMCC)
library(gstat)
library(sf)

# Unified Parallel and Distributed Processing in R for Everyone 
# – a library that facilitates parallel processing of point cloud data.
library(future)

# 1. Link to .las file ####

# Link to bigger .las file
#Link <- "E:/Remote Sensing Media/6. September 2025/Point Cloud/SU Lourensford September 2025_point cloud-001.copc.laz"

# Link to smaller 50cm/pixel .las file
Link <- "E:/Remote Sensing Media/6. September 2025/Point Cloud/SU Lourensford September 2025_point cloud_50cm.las"

# Link to Clipped .las file
Plot <- 37
Folder <- "Clipped"
Folder <- "Clipped_50cm"
Link <- paste0("E:/Remote Sensing Media/0. R Projects/Point Cloud/",Folder,"/clipped_Plot_",Plot,".las")

# You can filter attributes out if needed
las <- readLAS(Link, select = "xyzi")

# 2. Read and check point cloud ####

# Read .las file into memory
# las <- readLAS(Link)

# Check .las details
print(las)

# 3. How to clip a point cloud ####

# Region of interest (ROI)
las <- clip_circle(las, x = -6985, y = -3766340, radius = 5)
plot(las, bg = "white", size = 4)

# Find extent through plot boundaries

PlotNumber <- 37
shape_file <- st_read(paste0("E:/Remote Sensing Media/1. QGIS Projects/Michelle/Michelle QGIS Project/1 September 2025/Plot ",PlotNumber,".shp"))
extent <- ext(shape_file)
las <- clip_rectangle(las, xleft = extent[1], xright = extent[2], ybottom = extent[3], ytop = extent[4])
writeLAS(las, file.path("E:/Remote Sensing Media/0. R Projects/Point Cloud/Clipped_50cm", paste0("clipped_Plot_",PlotNumber,".las")), index = FALSE)

# 4. Ground classification ####

# Progressive Morphological Filter (PMF) 
ws <- seq(3, 12, 3)
th <- seq(0.1, 1.5, length.out = length(ws))
las1 <- classify_ground(las, algorithm = pmf(ws = ws, th = th))
plot(las1, color = "Classification", size = 3, bg = "lightblue")

# Multiscale Curvature Classification (MCC)
# Preferred but takes longer
las1 <- classify_ground(las, mcc(1.5,0.3))
plot(las1, color = "Classification", size = 3, bg = "lightblue") 


# Display ground points
gnd <- filter_ground(las1)
plot(gnd, size = 3, bg = "white")


# 5. Digital terrain model ####

# Preferred DTM model with the Triangular irregular network (TIN) algorithm
dtm_tin <- rasterize_terrain(las1, res = 0.1, algorithm = tin())
plot_dtm3d(dtm_tin, bg = "white") 

# DTM model with the Kriging algorithm
dtm_kriging_1 <- rasterize_terrain(las1, algorithm = kriging(k = 40))
plot_dtm3d(dtm_kriging_1, bg = "white") 


# 6. Height normalization ####

# Normalised 
nlas <- normalize_height(las1, tin())
plot(nlas, size = 4, bg = "white")


# 7. Digital Surface Model (DSM) and Canopy Height model (CHM) ####

# Digital Surface Model (DSM) and Canopy Height model (CHM)
chm <- rasterize_canopy(nlas, algorithm = p2r())
col <- height.colors(25)
plot(chm, col = col)

# Rasterize canopy with interpolation
chm <- rasterize_canopy(nlas, res = 0.5, p2r(0.2, na.fill = tin()))
plot(chm, col = col)

# 7.4 Post-processing a CHM

fill.na <- function(x, i=5) { if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }}
w <- matrix(1, 3, 3)

chm <- rasterize_canopy(nlas, res = 0.01, algorithm = p2r(subcircle = 0.15), pkg = "terra")
filled <- terra::focal(chm, w, fun = fill.na)
smoothed <- terra::focal(chm, w, fun = mean, na.rm = TRUE)

chms <- c(chm, filled, smoothed)
names(chms) <- c("Base", "Filled", "Smoothed")
plot(chms, col = col)


# 8. Individual tree detection and segmentation ####

# 8.1 Individual Tree Detection (ITD)

MinimumTreeHeight <- 0.5

# Local Maximum Filter with variable windows size
# f <- function(z) {1 * z + 0.5}
# lmf_algorithm <- lmf(ws = f, hmin = MinimumTreeHeight, shape = "circular")

# create Local Maximum Filter (lmf) function for the "ws" search
lmf_algorithm <- lmf(ws = 1, hmin = MinimumTreeHeight, shape = "circular")

# Locate trees in a circle with a diameter of "ws" in meters
ttops <- locate_trees(las = nlas[nlas$Z>= 0], algorithm = lmf_algorithm)

# Tree detection results in 2D
plot(chm, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)

# Tree detection results can also be visualized in 3D!
x <- plot(nlas, bg = "white", size = 4)
add_treetops3d(x, ttops, radius = 0.1, fastTransparency = TRUE, alpha = 1)

# Individual Tree Segmentation (ITS)

algo <- dalponte2016(chm, ttops)
las_segmented <- segment_trees(nlas, algo) # segment point cloud
plot(las_segmented, bg = "white", size = 4, color = "treeID") # visualize trees

# 9. Tree metrics ####

trees <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/QGIS Extracted data/1 September 2025/Plot 37.shp")
tree_heights <- terra::extract(chm, trees, fun = max, na.rm = TRUE)






