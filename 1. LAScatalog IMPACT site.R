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

# Change this single variable for each new batch!
date_folder <- "22. 08 April 2026"

# Read in all point clouds and shape files ####

# Extract the date part 
file_date <- sub("^\\d+\\.\\s*", "", date_folder)
file_date_safe <- gsub(" ", "_", file_date)

# --- CREATE MISSING IMPACT DIRECTORIES ---
# List all the required output directories for the IMPACT workflow
impact_dirs <- c(
  paste0("E:/Remote Sensing Media/", date_folder, "/04. Point Clouds Clipped IMPACT"),
  paste0("E:/Remote Sensing Media/", date_folder, "/05. Point Clouds Ground Classified IMPACT"),
  paste0("E:/Remote Sensing Media/", date_folder, "/06. Point Clouds Normalised IMPACT"),
  paste0("E:/Remote Sensing Media/", date_folder, "/07. Canopy Height Models IMPACT"),
  paste0("E:/Remote Sensing Media/", date_folder, "/09. Crown Metrics IMPACT")
)

# Loop through the list and create any folders that do not exist
for (dir in impact_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("Created missing directory:", dir, "\n")
  }
}

# --- BATCH PROCESSING SETUP ---
las_folder <- paste0("E:/Remote Sensing Media/", date_folder, "/03. Point Clouds")
all_las_files <- list.files(las_folder, pattern = "\\.(las|laz)$", full.names = TRUE, ignore.case = TRUE)

# Identify Top and Bottom files based on filenames
top_files <- all_las_files[grepl("Top", basename(all_las_files), ignore.case = TRUE)]
bot_files <- all_las_files[grepl("Bottom", basename(all_las_files), ignore.case = TRUE)]

# Read in shape files for Top and Bottom cropping
IMPACT_Top <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LAScatalog Boundaries/Top IMPACT Boundaries.shp")
IMPACT_Bottom <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LAScatalog Boundaries/Bottom IMPACT Boundaries.shp")

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

tic()

if (length(top_files) > 0 && length(bot_files) > 0) {
  print("Top and Bottom point clouds detected. Cropping separately...")
  
  # Read into separate catalogs
  ctg_top <- readLAScatalog(top_files)
  ctg_bot <- readLAScatalog(bot_files)
  
  # Configure engine settings for TOP
  opt_independent_files(ctg_top) <- FALSE
  opt_select(ctg_top) <- "xyz"
  opt_output_files(ctg_top) <- paste0("E:/Remote Sensing Media/",date_folder,"/04. Point Clouds Clipped IMPACT/IMPACT_Site_Top_", file_date_safe)
  
  # Configure engine settings for BOTTOM
  opt_independent_files(ctg_bot) <- FALSE
  opt_select(ctg_bot) <- "xyz"
  opt_output_files(ctg_bot) <- paste0("E:/Remote Sensing Media/",date_folder,"/04. Point Clouds Clipped IMPACT/IMPACT_Site_Bottom_", file_date_safe)
  
  # Crop catalogs separately
  print("Cropping Top IMPACT Boundary...")
  ctg_clipped_top <- clip_roi(ctg_top, IMPACT_Top)
  
  print("Cropping Bottom IMPACT Boundary...")
  ctg_clipped_bot <- clip_roi(ctg_bot, IMPACT_Bottom)
  
  # Re-read the entire clipped directory as a single catalog for downstream steps
  ctg_clipped <- readLAScatalog(paste0("E:/Remote Sensing Media/", date_folder, "/04. Point Clouds Clipped IMPACT"))
  
} else {
  print("Single point cloud or no Top/Bottom distinction detected. Make sure your files contain 'Top' and 'Bottom' in the names.")
}

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

# Combine Canopy Height Models Tiffs and Extract tree heights ####

tic()

# Since rasterize_canopy output multiple tif tiles to disk based on chunks,
# we need to combine them virtually for exact_extract to read them as one continuous site.
chm_files <- list.files(paste0("E:/Remote Sensing Media/",date_folder,"/07. Canopy Height Models IMPACT/"), 
                        pattern = "\\.tif$", full.names = TRUE)
site_chm_vrt <- terra::vrt(chm_files)

# Write the virtual raster out to a single physical .tif file
single_chm_path <- paste0("E:/Remote Sensing Media/",date_folder,"/07. Canopy Height Models IMPACT/IMPACT_Site_CHM_Single_", file_date_safe, ".tif")
terra::writeRaster(site_chm_vrt, filename = single_chm_path, overwrite = TRUE)

# Optional: Delete the individual chunk tiles to save disk space
file.remove(chm_files)

# Re-read the single physical file for your extraction (replaces the VRT in memory)
site_chm_single <- terra::rast(single_chm_path)

# Calculate metrics using exact_extract against the NEW single continuous raster
trees$Tree_Height <- exact_extract(site_chm_single, trees, 'max')

# Save to shapefile
st_write(trees, paste0("E:/Remote Sensing Media/",date_folder,"/09. Crown Metrics IMPACT/Crown_Metrics_", file_date_safe, ".shp"), delete_dsn = TRUE)

# Save lightweight CSV
write.csv(st_drop_geometry(trees), paste0("E:/Remote Sensing Media/",date_folder,"/09. Crown Metrics IMPACT/Crown_Metrics_", file_date_safe, ".csv"), row.names = FALSE)
toc()