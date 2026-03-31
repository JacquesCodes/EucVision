library(lidR)
library(terra)
library(future)

# Enable parallel processing to speed up the 500 million points
# (Keep workers at 4 or 6 so you don't run out of memory per thread)
plan(multisession, workers = 8)

# 1. Load your baseline point cloud catalog
baseline_ctg <- readLAScatalog("E:/Remote Sensing Media/02. 01 September 2025 (DJI M300)/03. Point clouds")

# 2. Memory and File Management
opt_chunk_size(baseline_ctg) <- 25
opt_chunk_buffer(baseline_ctg) <- 5
opt_output_files(baseline_ctg) <- "E:/Remote Sensing Media/00. Baseline DTM and plot cropping/Classified_Chunks/Tile_{XLEFT}_{YBOTTOM}_ground"

# --- THE SPEED OPTIMIZATION ---
# Only load spatial coordinates into RAM, dropping heavy RGB data
opt_select(baseline_ctg) <- "xyz"
# ------------------------------

# 3. Classify the ground
baseline_classified_ctg <- classify_ground(baseline_ctg, csf(sloop_smooth = TRUE, 
                                                             class_threshold = 0.01, 
                                                             cloth_resolution = 0.5, 
                                                             time_step = 1))

# 4. Rasterize the terrain to create your master DTM
# rasterize_terrain handles catalogs automatically. 
# A 10ha raster at 10cm resolution is very small (~40MB), so it is safe to load into RAM.
master_dtm <- rasterize_terrain(baseline_classified_ctg, res = 0.1, algorithm = tin())

# 4. Save the final DTM to disk
writeRaster(master_dtm, "E:/Remote Sensing Media/00. Baseline DTM and plot cropping/Master_Baseline_DTM.tif", overwrite = TRUE)