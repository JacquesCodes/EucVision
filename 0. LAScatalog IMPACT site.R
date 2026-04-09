library(lidR)
library(RCSF)
library(RMCC)
library(gstat)
library(sf)
library(tictoc)
library(geometry)
library(dplyr)
library(future)
library(sp)
library(terra)
library(exactextractr)

# Read in all point clouds and shape files ####

# Change this single variable for each new batch!
date_folder <- "21. 31 March 2026"

# Extract the date part 
file_date <- sub("^\\d+\\.\\s*", "", date_folder)
file_date_safe <- gsub(" ", "_", file_date)

# Read in point clouds into a catalog (ctg)
ctg <- readLAScatalog(paste0("E:/Remote Sensing Media/",date_folder,"/03. Point Clouds"))

# Read in shape files for individual plot boundaries
plots_buffered_unsorted <- st_read(paste0("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LidR Boundaries/EucVision LidR Boundaries.shp"))
plots <- plots_buffered_unsorted[order(plots_buffered_unsorted$id), ]

IMPACT_Boundaries <- st_read("E:/Remote Sensing Media/00. Baseline DTM and Plot Cropping/Impact Boundaries.shp")

# Read in tree shape files for height extraction
trees <- st_read(paste0("E:/Remote Sensing Media/", date_folder, "/08. Crown Polygons/Crown_Polygons_", file_date_safe, ".shp"))

# Automatically check and transform to EPSG: 2048 if it doesn't match
if (is.na(st_crs(trees)$epsg) || st_crs(trees)$epsg != 2048) {
  trees <- st_transform(trees, 2048)
  print("Transformed CRS to 2048 successfully.")
}

trees <- trees %>%
  select(-any_of(c("group_ulid", "N_GM", "id", "N_FG", "N_BG", "BBox")))

# Crop IMPACT site ####

plan(multisession)
opt_independent_files(ctg) <- FALSE
opt_select(ctg) <- "xyz"

tic()
# Changed output naming since it's likely a single boundary polygon now, not 75 plots with IDs
opt_output_files(ctg) <- paste0("E:/Remote Sensing Media/",date_folder,"/04. Point Clouds Clipped IMPACT/IMPACT_Site_", file_date_safe)
ctg_clipped <- clip_roi(ctg, IMPACT_Boundaries)
toc()

# Classify plots ####

plan(multisession, workers = 6)
opt_independent_files(ctg_clipped) <- TRUE
opt_select(ctg_clipped) <- "xyz"

# CRITICAL FOR WHOLE SITE PROCESSING: Define chunk size and buffer
# This splits the single large file across your 6 CPU workers and prevents edge artifacts
opt_chunk_size(ctg_clipped) <- 200 # Process in 200m x 200m chunks
opt_chunk_buffer(ctg_clipped) <- 10  # 10m buffer around chunks to prevent edge artifacts

tic()
opt_output_files(ctg_clipped) <- paste0("E:/Remote Sensing Media/",date_folder,"/05. Point Clouds Ground Classified IMPACT/", "IMPACT_Tile_{XLEFT}_{YBOTTOM}_classified", file_date_safe)
ctg_classified <- classify_ground(ctg_clipped, csf(sloop_smooth = TRUE, 
                                                   class_threshold = 0.15, 
                                                   cloth_resolution = 1.5, 
                                                   rigidness = 3,
                                                   time_step = 1))
toc()

# # Generate Digital Terrain Model (DTM) ####
# 
# # Ensure you create this folder in your directory first!
# # e.g., "E:/Remote Sensing Media/19. 16 March 2026/10. Digital Terrain Models/"
# 
# ctg_classified <- readLAScatalog(paste0("E:/Remote Sensing Media/",date_folder,"/05. Point Clouds Ground Classified"))
# 
# plan(multisession, workers = 6)
# opt_independent_files(ctg_classified) <- TRUE
# opt_select(ctg_classified) <- "xyzc"
# 
# # Maintain the exact same chunking and buffers used for classification
# opt_chunk_size(ctg_classified) <- 200 
# opt_chunk_buffer(ctg_classified) <- 10  
# 
# tic()
# opt_output_files(ctg_classified) <- paste0("E:/Remote Sensing Media/",date_folder,"/10. Digital Terrain Models/", "IMPACT_Tile_{XLEFT}_{YBOTTOM}_dtm_", file_date_safe)
# 
# # Generate DTM. rasterize_terrain automatically targets class 2 (ground) points.
# # A resolution of 0.5m (50cm) is usually perfect for DTMs to smooth out micro-noise,
# # but you can change 'res' to 0.05 if you want it to perfectly match your CHM pixel-for-pixel.
# ctg_dtm <- rasterize_terrain(ctg_classified, 
#                              res = 0.5, 
#                              algorithm = tin())
# print("Rasterize DTM time:")
# toc()

# Normalize plots ####

plan(multisession, workers = 6)
opt_independent_files(ctg_classified) <- TRUE
opt_select(ctg_classified) <- "xyzc"

# Maintain chunking and buffers for normalization
opt_chunk_size(ctg_classified) <- 200 
opt_chunk_buffer(ctg_classified) <- 10 

tic()
opt_output_files(ctg_classified) <- paste0("E:/Remote Sensing Media/",date_folder,"/06. Point Clouds Normalised IMPACT/", "IMPACT_Tile_{XLEFT}_{YBOTTOM}_normalised", file_date_safe)
ctg_normalised <- normalize_height(las = ctg_classified, algorithm = tin())
toc()

# Rasterize plots ####

ctg_normalised <- readLAScatalog(paste0("E:/Remote Sensing Media/",date_folder,"/06. Point Clouds Normalised IMPACT"))

plan(multisession, workers = 4)
opt_independent_files(ctg_normalised) <- TRUE
opt_select(ctg_normalised) <- "xyz"
opt_filter(ctg_normalised) <- "-drop_z_below 0 -drop_z_above 30"

# Maintain chunking and buffers for rasterization
opt_chunk_size(ctg_normalised) <- 200 
opt_chunk_buffer(ctg_normalised) <- 10 

tic()
opt_output_files(ctg_normalised) <- paste0("E:/Remote Sensing Media/",date_folder,"/07. Canopy Height Models IMPACT/", "IMPACT_Tile_{XLEFT}_{YBOTTOM}_chm", file_date_safe)
ctg_chm <- rasterize_canopy(ctg_normalised,
                            res = 0.05,
                            algorithm = p2r(na.fill = tin()))
print("Rasterize canopy time:")
toc()

# Extract tree heights ####

tic()

# NEW: Since rasterize_canopy output multiple tif tiles to disk based on chunks,
# we need to combine them virtually for exact_extract to read them as one continuous site.
chm_files <- list.files(paste0("E:/Remote Sensing Media/",date_folder,"/07. Canopy Height Models IMPACT/"), 
                        pattern = "\\.tif$", full.names = TRUE)
site_chm_vrt <- terra::vrt(chm_files)

# Calculate metrics using exact_extract against the virtual raster
trees$Tree_Height <- exact_extract(site_chm_vrt, trees, 'max')

# Save to shapefile
st_write(trees, paste0("E:/Remote Sensing Media/",date_folder,"/09. Crown Metrics IMPACT/Crown_Metrics_", file_date_safe, ".shp"), delete_dsn = TRUE)

# Save lightweight CSV
write.csv(st_drop_geometry(trees), paste0("E:/Remote Sensing Media/",date_folder,"/09. Crown Metrics IMPACT/Crown_Metrics_", file_date_safe, ".csv"), row.names = FALSE)
toc()