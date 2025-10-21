library(lidR)
library(RCSF)
library(RMCC)
library(gstat)
library(sf)

# Link to bigger .las file
#Link <- "E:/Remote Sensing Media/6. September 2025/Point Cloud/SU Lourensford September 2025_point cloud-001.copc.laz"

# Clipped .las
# Link <- "E:/Remote Sensing Media/6. September 2025/Point Cloud/clipped_Plot_37.las"

# Link to smaller 50cm/pixel .las file
Link <- "E:/Remote Sensing Media/6. September 2025/Point Cloud/SU Lourensford September 2025_point cloud_50cm.las"

# You can filter attributes out if needed
las <- readLAS(Link, select = "xyzi")

# Check .las details
print(las)

# Region of interest (ROI)
las <- clip_circle(las, x = -6990, y = -3766340, radius = 10)
plot(las, bg = "white", size = 4)

# How to clip a point cloud ####

# # Find extent through plot boundaries
# shape_file <- st_read("E:/Remote Sensing Media/1. QGIS Projects/Michelle/Michelle QGIS Project/1 September 2025/Plot 37.shp")
# clipped_las <- clip_rectangle(las, xleft = -6988.282, xright = -6984.501, ybottom=-3766353, ytop=-3766326)
# writeLAS(clipped_las, file.path("E:/Remote Sensing Media/6. September 2025/Point Cloud", 'clipped_Plot_37.las'), index = FALSE)

## 4.1 Progressive Morphological Filter (PMF) for ground classification ####
ws <- seq(3, 12, 3)
th <- seq(0.1, 1.5, length.out = length(ws))
las1 <- classify_ground(las, algorithm = pmf(ws = ws, th = th))
plot(las1, color = "Classification", size = 3, bg = "lightblue") 

# Display ground points
gnd <- filter_ground(las1)
plot(gnd, size = 3, bg = "white")

# DTM model with the TIN algorithm
dtm_tin <- rasterize_terrain(las1, res = 0.1, algorithm = tin())
plot_dtm3d(dtm_tin, bg = "white") 

# DTM model with the Kriging algorithm
dtm_kriging_1 <- rasterize_terrain(las1, algorithm = kriging(k = 40))
plot_dtm3d(dtm_kriging_1, bg = "white") 

# Normalised 
nlas <- normalize_height(las1, tin())
plot(nlas, size = 4, bg = "white")

# Digital Surface Model (DSDM) and Canopy Height model (CHM)
chm <- rasterize_canopy(nlas, algorithm = p2r())
col <- height.colors(25)
plot(chm, col = col)

# Rasterize canopy with interpolation
chm <- rasterize_canopy(nlas, res = 0.5, p2r(0.2, na.fill = tin()))
plot(chm, col = col)

## 7.4 Post-processing a CHM ####

fill.na <- function(x, i=5) { if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }}
w <- matrix(1, 3, 3)

chm <- rasterize_canopy(nlas, res = 0.5, algorithm = p2r(subcircle = 0.15), pkg = "terra")
filled <- terra::focal(chm, w, fun = fill.na)
smoothed <- terra::focal(chm, w, fun = mean, na.rm = TRUE)

chms <- c(chm, filled, smoothed)
names(chms) <- c("Base", "Filled", "Smoothed")
plot(chms, col = col)

# 8 Individual tree detection and segmentation ####

## 8.1 Individual Tree Detection (ITD) ####

# Locate trees in a circle with a diameter of "ws" in meters
ttops <- locate_trees(nlas, lmf(ws = 1))

chm <- rasterize_canopy(nlas, res = 1, p2r(0.2, na.fill = tin()))
plot(chm, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)

# Tree detection results can also be visualized in 3D!
x <- plot(nlas, bg = "white", size = 4)
add_treetops3d(x, ttops, radius = 0.1, fastTransparency = TRUE, alpha = 0.5)

# Local Maximum Filter with variable windows size

f <- function(x) {x * 0.2 + 0.2}
heights <- seq(0,4,1)
ws <- f(heights)
plot(heights, ws, type = "l", ylim = c(0,6))

ttops <- locate_trees(nlas, lmf(f))

plot(chm, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)

# Tree detection results can also be visualized in 3D!
x <- plot(nlas, bg = "white", size = 4)
add_treetops3d(x, ttops, radius = 0.1, fastTransparency = TRUE, alpha = 0.5)

# Individual Tree Segmentation (ITS)

algo1 <- dalponte2016(chm, ttops)
algo2 <- li2012()
las <- segment_trees(las, algo1, attribute = "IDdalponte")
las <- segment_trees(las, algo2, attribute = "IDli")

x <- plot(las, bg = "white", size = 4, color = "IDdalponte", colorPalette = pastel.colors(200))
#> The argument 'coloPalette' is deprecated. Use 'pal' instead
plot(las, add = x + c(100,0), bg = "white", size = 4, color = "IDli", colorPalette = pastel.colors(200))
#> The argument 'coloPalette' is deprecated. Use 'pal' instead