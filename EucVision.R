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
# DO NOT USE .laz ONLY USE .las IT IMPROVES PERFORMANCE x10!
################################################################################


# Top Catalog(ctg)
ctg <- readLAScatalog("E:/Remote Sensing Media/7. 30 October 2025/Point cloud/")

las_check(ctg)

plot(ctg)

opt_output_files(ctg) <- "E:/Remote Sensing Media/0. R Projects/2. 30 October 2025/Test/{XLEFT}_{YBOTTOM}" # label outputs based on coordinates
opt_chunk_buffer(ctg) <- 0
opt_chunk_size(ctg) <- 50 # retile to 250 m
small <- catalog_retile(ctg) # apply retile
plot(small) # some plotting

# Multiple threads mode
plan(multisession)

# Each plot is independent and buffers are not needed
opt_independent_files(small) <- TRUE

tic()
# Write to disk rather than memory:
opt_output_files(small) <- paste0("E:/Remote Sensing Media/0. R Projects/2. 30 October 2025/Test/", "{*}_classified")
# Ground classifications :
small_classified <- classify_ground(small, csf(sloop_smooth = TRUE, class_threshold = 0.01, cloth_resolution = 0.5, time_step = 1))
toc()

tic()
# Write to disk rather than memory:
opt_output_files(small_classified) <- paste0("E:/Remote Sensing Media/0. R Projects/2. 30 October 2025/Test/", "{*}_normalised")
# A point cloud-based normalization without a raster:
small_normalised <- normalize_height(las = small_classified, algorithm = tin())
toc()

tic()
# Write to disk rather than memory:
opt_output_files(small_normalised) <- paste0("E:/Remote Sensing Media/0. R Projects/2. 30 October 2025/Test/", "{*}_chm")
# Rasterize canopy with interpolation:
small_chm <- rasterize_canopy(small_normalised, res = 0.01, algorithm = p2r(na.fill = tin()))
print("Rasterize canopy time:")
toc()
