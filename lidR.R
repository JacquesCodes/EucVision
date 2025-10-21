library(lidR)
library(RCSF)
library(RMCC)
library(gstat)
library(sf)

# Website: https://r-lidar.github.io/lidRbook/

# 2. Initiate .las ####

# Link to smaller 50cm/pixel .las file
#Link <- "E:/Remote Sensing Media/6. September 2025/Point Cloud/SU Lourensford September 2025_point cloud_50cm.las"

# Link to bigger .las file
Link <- "E:/Remote Sensing Media/6. September 2025/Point Cloud/SU Lourensford September 2025_point cloud-001.copc.laz"

# 2.1 Read .las file
#las <- readLAS(Link)

# You can filter attributes out if needed
las <- readLAS(Link, select = "xyzi")

# Check .las details
print(las)

# 2.2 Check if Point cloud is valid
las_check(las)

# How to clip a point cloud ####

# Region of interest (ROI)
roi <- clip_circle(las, x = -6990, y = -3766340, radius = 10)
plot(roi, bg = "white", size = 4)


shape_file <- st_read("E:/Remote Sensing Media/1. QGIS Projects/Michelle/Michelle QGIS Project/1 September 2025/Plot 37.shp")

clipped_las <- clip_rectangle(las, xleft = -6988.282, xright = -6984.501, ybottom=-3766353, ytop=-3766326)

writeLAS(clipped_las, file.path("E:/Remote Sensing Media/6. September 2025/Point Cloud", 'clipped_file.las'), index = FALSE)

# 3. Render .las ####

# 3.1 Basic 3D rendering with rgl
plot(las, bg = "lightblue", color = "RGB", size = 3)

# Size = Size of each point

# 4. Ground Classification ####

## 4.1 Progressive Morphological Filter ####

ws <- seq(3, 12, 3)
th <- seq(0.1, 1.5, length.out = length(ws))
las1 <- classify_ground(las, algorithm = pmf(ws = ws, th = th))
plot(las1, color = "Classification", size = 3, bg = "lightblue") 

## Display ground points

gnd <- filter_ground(las1)
plot(gnd, size = 3, bg = "white")


## 4.2 Cloth Simulation Function ####

mycsf <- csf(sloop_smooth = TRUE, class_threshold = 1, cloth_resolution = 1, time_step = 1)
las2 <- classify_ground(las, mycsf)
plot(las2, color = "Classification", size = 3, bg = "lightblue") 


## 4.3 Multiscale Curvature Classification (MCC) ####
# Preferred

las3 <- classify_ground(roi, mcc(1.5,0.3))
plot(las3, color = "Classification", size = 3, bg = "lightblue") 


# 5. Digital terrain model ####

# Topography .las example
LASfile <- system.file("extdata", "Topography.laz", package="lidR")
las0 <- readLAS(LASfile, select = "xyzc")
plot(las, size = 3, bg = "white")

## 5.1 Triangular irregular network (TIN) ####
# Preferred

# DTM model with the TIN algorithm for the example
dtm_tin_0 <- rasterize_terrain(las0, res = 1, algorithm = tin())
plot_dtm3d(dtm_tin_0, bg = "white") 

# DTM model with the TIN algorithm for my .las with classified ground points,
# Ground points classified with 1. Progressive Morphological Filter
dtm_tin_1 <- rasterize_terrain(las1, res = 1, algorithm = tin())
plot_dtm3d(dtm_tin_1, bg = "white") 

## 5.2 Invert distance weighting (IDW) ####

dtm_idw_1 <- rasterize_terrain(las1, algorithm = knnidw(k = 10L, p = 2))
plot_dtm3d(dtm_idw_1, bg = "white") 

## 5.3 Kriging ####

dtm_kriging_1 <- rasterize_terrain(las1, algorithm = kriging(k = 40))
plot_dtm3d(dtm_kriging_1, bg = "white") 

# 6. Height normalization ####

# # normalized las (nlas) = las - dtm
# nlas <- las - dtm_tin_1
# 
# # Note that not all flat surface's z value is 0
# hist(filter_ground(nlas)$Z, breaks = seq(-1000, 1000, 1), main = "", xlab = "Elevation")

# How to fix Z values
# PREFFERED way to normalize heights. Don't las - dtm!
nlas <- normalize_height(las3, tin())
hist(filter_ground(nlas)$Z, breaks = seq(-0.6, 0.6, 0.01), main = "", xlab = "Elevation")

# Normalised Heat map
plot(nlas, size = 4, bg = "white")

# Normalised RGB map
plot(nlas, color = "RGB", size = 4, bg = "white")


# How to clip a point cloud ####

# Region of interest (ROI)
roi <- clip_circle(nlas, x = -6990, y = -3766340, radius = 30)
plot(roi, bg = "white", size = 4)

# 7. Digital Surface Model (DSDM) and Canopy Height model (CHM) ####

chm <- rasterize_canopy(roi, res = 1, algorithm = p2r())
col <- height.colors(25)
plot(chm, col = col)

# Rasterize canopy with interpolation
chm <- rasterize_canopy(roi, res = 0.5, p2r(0.2, na.fill = tin()))
plot(chm, col = col)

## 7.4 Post-processing a CHM ####

fill.na <- function(x, i=5) { if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }}
w <- matrix(1, 3, 3)

chm <- rasterize_canopy(roi, res = 0.5, algorithm = p2r(subcircle = 0.15), pkg = "terra")
filled <- terra::focal(chm, w, fun = fill.na)
smoothed <- terra::focal(chm, w, fun = mean, na.rm = TRUE)

chms <- c(chm, filled, smoothed)
names(chms) <- c("Base", "Filled", "Smoothed")
plot(chms, col = col)

# 8 Individual tree detection and segmentation ####

## 8.1 Individual Tree Detection (ITD) ####

# Locate trees in a circle with a diameter of "ws" in meters
ttops <- locate_trees(roi, lmf(ws = 1))

chm_roi <- rasterize_canopy(roi, res = 0.5, p2r(0.2, na.fill = tin()))
plot(chm_roi, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)

# Tree detection results can also be visualized in 3D!
x <- plot(roi, bg = "white", size = 4)
add_treetops3d(x, ttops)

# Local Maximum Filter with variable windows size

f <- function(x) {x * 0.2 + 0.2}
heights <- seq(0,4,1)
ws <- f(heights)
plot(heights, ws, type = "l", ylim = c(0,6))

ttops <- locate_trees(las, lmf(f))

plot(chm, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)

# Segment point cloud

algo1 <- dalponte2016(chm, ttops)
algo2 <- li2012()
las_algo1 <- segment_trees(nlas, algo1, attribute = "IDdalponte")
las_algo2 <- segment_trees(nlas, algo2, attribute = "IDli")

x <- plot(las_algo1, bg = "white", size = 4, color = "IDdalponte", colorPalette = pastel.colors(200))
#> The argument 'coloPalette' is deprecated. Use 'pal' instead
plot(las, add = x + c(100,0), bg = "white", size = 4, color = "IDli", colorPalette = pastel.colors(200))
#> The argument 'coloPalette' is deprecated. Use 'pal' instead


x <- plot(las_algo2, bg = "white", size = 4, color = "IDdalponte", colorPalette = pastel.colors(200))
#> The argument 'coloPalette' is deprecated. Use 'pal' instead
plot(las, add = x + c(100,0), bg = "white", size = 4, color = "IDli", colorPalette = pastel.colors(200))
#> The argument 'coloPalette' is deprecated. Use 'pal' instead