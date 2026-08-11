library(terra)

# 1. Load the skewed February raster
skewed_ortho <- rast("E:/Remote Sensing Media/01. 25 February 2025/02. Digital Surface Models/LFDSMFeb25-05cm.tif")

# 2. Create the data frame of your stable tie points (in EPSG:4326)
pts_df <- data.frame(
  id = 1:3,
  ref_x    = c(19.075345557, 19.075358271, 19.075358381),
  ref_y    = c(-34.023412141, -34.023413146, -34.023421714), 
  skew_x   = c(19.075356327, 19.075368787, 19.075369085), 
  skew_y   = c(-34.023419667, -34.023420726, -34.023429354)  
)

# 3. Convert coordinates into spatial objects (SpatVectors) in EPSG:4326
ref_pts <- vect(cbind(pts_df$ref_x, pts_df$ref_y), crs="epsg:4326")
skew_pts <- vect(cbind(pts_df$skew_x, pts_df$skew_y), crs="epsg:4326")

# 4. Project the points to match the raster's exact metric coordinate system
ref_pts_proj <- project(ref_pts, crs(skewed_ortho))
skew_pts_proj <- project(skew_pts, crs(skewed_ortho))

# 5. Extract the new metric coordinates
coords_ref <- crds(ref_pts_proj)
coords_skew <- crds(skew_pts_proj)

# 6. Calculate the true metric offset with the SA Axis Fix
# By reversing the subtraction on the X-axis, we flip the Westing difference
# into a standard shift (Negative value corrects the horizontal drift).
dx <- mean(coords_skew[,1] - coords_ref[,1]) 
dy <- mean(coords_ref[,2] - coords_skew[,2]) # Keep as is (Negative moves it correctly North)

print(paste("Fixed X offset in meters:", round(dx, 3)))
print(paste("Fixed Y offset in meters:", round(dy, 3)))

# 7. Apply the corrected metric shift to the raster
aligned_ortho <- shift(skewed_ortho, dx = dx, dy = dy)

# 8. Export the correctly aligned orthomosaic
writeRaster(aligned_ortho, "E:/Remote Sensing Media/01. 25 February 2025/02. Digital Surface Models/DSM_25_February_2025.tif", overwrite=TRUE)