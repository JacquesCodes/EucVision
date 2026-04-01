library(lidR)
library(terra)
library(future)

# Enable parallel processing to speed up the point clouds processing
# (Keep workers at 4 or 6 so you don't run out of memory per thread)
plan(multisession, workers = 8)

# 1. Load your baseline point cloud catalog
baseline_ctg <- readLAScatalog("E:/Remote Sensing Media/03. 30 October 2025/03. Point clouds")

# 2. Memory and File Management
opt_chunk_size(baseline_ctg) <- 50
opt_chunk_buffer(baseline_ctg) <- 10
opt_output_files(baseline_ctg) <- "E:/Remote Sensing Media/00. Baseline DTM and plot cropping/Classified_Chunks_2/Tile_{XLEFT}_{YBOTTOM}_ground"

# --- THE SPEED OPTIMIZATION ---
# Only load spatial coordinates into RAM, dropping heavy RGB data
opt_select(baseline_ctg) <- "xyz"
# ------------------------------

# 3. Classify the ground
baseline_classified_ctg <- classify_ground(baseline_ctg, csf(sloop_smooth = TRUE,
                                                             class_threshold = 0.01,
                                                             cloth_resolution = 2,
                                                             time_step = 1))


# Resolution 20 recommended for forestry, 10 for mountains and 50 for cities
# baseline_classified_ctg <- classify_ground(baseline_ctg, ptd(20))


# 4. Rasterize the terrain to create your master DTM

baseline_classified_ctg <- readLAScatalog("E:/Remote Sensing Media/00. Baseline DTM and plot cropping/Classified_Chunks_30_October_2025")

# rasterize_terrain handles catalogs automatically. 
# A 10ha raster at 10cm resolution is very small (~40MB), so it is safe to load into RAM.
master_dtm <- rasterize_terrain(baseline_classified_ctg,
                                res = 0.05,
                                algorithm = tin())

# 4. Save the final DTM to disk
writeRaster(master_dtm, "E:/Remote Sensing Media/00. Baseline DTM and plot cropping/Master_Baseline_DTM_30_October_2025_2.tif", overwrite = TRUE)




library(terra)

# 1. Define the file path
path <- "E:/Remote Sensing Media/00. Baseline DTM and Plot Cropping/Master_Baseline_DTM_30_October_2025.tif"

# 2. Load the raster
# Note: Use forward slashes (/) or double backslashes (\\) in R paths
dtm <- rast(path)

# 3. Basic inspection
print(dtm)

# 4. Simple 2D Visualization
plot(dtm, main = "Baseline DTM - October 2025", col = terrain.colors(100))

# Calculate slope and aspect
slp <- terrain(dtm, "slope", unit = "radians")
asp <- terrain(dtm, "aspect", unit = "radians")

# Calculate hillshade
hill <- shade(slp, asp, angle = 45, direction = 315)

# Plot hillshade with the DTM overlaid (transparent)
plot(hill, col = grey(0:100/100), legend = FALSE, main = "DTM with Hillshade")
plot(dtm, col = terrain.colors(100, alpha = 0.5), add = TRUE)


