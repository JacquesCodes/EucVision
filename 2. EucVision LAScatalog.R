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
date_folder <- "19. 16 March 2026"

# Read in all point clouds and shape files ####

# Extract the date part by removing the leading folder number, dot, and space
# This turns "17. 02 March 2026" into "02 March 2026"
file_date <- sub("^\\d+\\.\\s*", "", date_folder)

# Replace spaces with underscores for safer file naming conventions
# This turns "02 March 2026" into "02_March_2026"
file_date_safe <- gsub(" ", "_", file_date)

# Read in shape files for individual plot boundaries
plots_buffered_unsorted <- st_read(paste0("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LidR Boundaries/EucVision LidR Boundaries.shp"))
plots <- plots_buffered_unsorted[order(plots_buffered_unsorted$id), ]

# --- BATCH PROCESSING & CATALOG SETUP ---
las_folder <- paste0("E:/Remote Sensing Media/", date_folder, "/03. Point Clouds")
all_las_files <- list.files(las_folder, pattern = "\\.(las|laz)$", full.names = TRUE, ignore.case = TRUE)

# Identify Top and Bottom files based on filenames
top_files <- all_las_files[grepl("Top", basename(all_las_files), ignore.case = TRUE)]
bot_files <- all_las_files[grepl("Bottom", basename(all_las_files), ignore.case = TRUE)]

# Read in tree shape files for height extraction ####
trees <- st_read(paste0("E:/Remote Sensing Media/", date_folder, "/08. Crown Polygons/Crown_Polygons_", file_date_safe, ".shp"))

# Automatically check and transform to EPSG: 2048 if it doesn't match
if (is.na(st_crs(trees)$epsg) || st_crs(trees)$epsg != 2048) {
  trees <- st_transform(trees, 2048)
  print("Transformed CRS to 2048 successfully.")
}

trees <- trees %>%
  select(-any_of(c("group_ulid", "N_GM", "id", "N_FG", "N_BG", "BBox")))

# # View catalog, plots and trees
# plot(ctg)
# plot(plots$geometry,add = TRUE)
# plot(trees$geometry, add = TRUE, col = "red")

# Crop plots ####

# Multiple CPU threads mode
plan(multisession)

# tic() & toc() is to check code running time
tic()

# Check if we need to split processing based on Top/Bottom files
if (length(top_files) > 0 && length(bot_files) > 0) {
  print("Top and Bottom point clouds detected. Splitting processing by Plot ID...")
  
  # Split plot boundaries by ID
  plots_top <- plots %>% filter(id <= 21)
  plots_bot <- plots %>% filter(id >= 22)
  
  # Read into separate catalogs
  ctg_top <- readLAScatalog(top_files)
  ctg_bot <- readLAScatalog(bot_files)
  
  # Configure engine settings for TOP
  opt_independent_files(ctg_top) <- FALSE
  opt_select(ctg_top) <- "xyz"
  # Note: using {id} instead of {ID} to match your dataframe's column name casing
  opt_output_files(ctg_top) <- paste0("E:/Remote Sensing Media/", date_folder, "/04. Point Clouds Clipped/", "Plot_{id}_", file_date_safe)
  
  # Configure engine settings for BOTTOM
  opt_independent_files(ctg_bot) <- FALSE
  opt_select(ctg_bot) <- "xyz"
  opt_output_files(ctg_bot) <- paste0("E:/Remote Sensing Media/", date_folder, "/04. Point Clouds Clipped/", "Plot_{id}_", file_date_safe)
  
  # Crop catalogs separately
  print("Cropping Top plots (1-21)...")
  ctg_clipped_top <- clip_roi(ctg_top, plots_top)
  
  print("Cropping Bottom plots (22-75)...")
  ctg_clipped_bot <- clip_roi(ctg_bot, plots_bot)
  
  # Re-read the entire clipped directory as a single catalog for downstream steps
  ctg_clipped <- readLAScatalog(paste0("E:/Remote Sensing Media/", date_folder, "/04. Point Clouds Clipped"))
  
} else {
  print("Single point cloud or no Top/Bottom distinction detected. Processing entirely...")
  
  # Normal processing for a single or unstructured catalog
  ctg <- readLAScatalog(las_folder)
  opt_independent_files(ctg) <- FALSE
  opt_select(ctg) <- "xyz"
  opt_output_files(ctg) <- paste0("E:/Remote Sensing Media/", date_folder, "/04. Point Clouds Clipped/", "Plot_{id}_", file_date_safe)
  
  # Crop plots
  ctg_clipped <- clip_roi(ctg, plots)
}

toc()

# Classify plots ####

plan(multisession)
opt_independent_files(ctg_clipped) <- TRUE
opt_select(ctg_clipped) <- "xyz"

tic()
# Write to disk rather than memory:
opt_output_files(ctg_clipped) <- paste0("E:/Remote Sensing Media/",date_folder,"/05. Point Clouds Ground Classified/", "{*}_classified")
# Ground classifications :
ctg_classified <- classify_ground(ctg_clipped, csf(sloop_smooth = TRUE, 
                                                   class_threshold = 0.01, 
                                                   cloth_resolution = 1, 
                                                   time_step = 1))
toc()

# Class_threshold 
# => The distance to the simulated cloth to classify a point cloud into ground and non-ground. 
# The default is 0.5. 
# Need to be set no larger than the smallest tree. 0.01 preferred for best height estimations.
# The higher the value the higher the ground classifications become. 

# Cloth_resolution 
# => The distance between particles in the cloth. 
# This is usually set to the average distance of the points in the point cloud. 
# The default value is 0.5.
# PREFFERED = 1
# DO NOT MAKE LOWER THAN 1 as it classify trees as ground points because of canopy closure
# Needed to make value to 1 for 1m x 1m plots 
# Otherwise the cloth falls between the points and classify trees.

# Normalize plots ####

# If in the future the user are unable to classify ground. Use a baseline raster from previous datasets:
# ctg_normalised <- normalize_height(las = ctg_classified, algorithm = tin(), dtm = )

plan(multisession)
opt_independent_files(ctg_classified) <- TRUE
# Only load in x-,y-,z- coordinates and "c" the classification values into RAM
opt_select(ctg_classified) <- "xyzc"

tic()
# Write to disk rather than memory:
opt_output_files(ctg_classified) <- paste0("E:/Remote Sensing Media/",date_folder,"/06. Point Clouds Normalised/", "{*}_normalised")
# A point cloud-based normalization without a raster:
ctg_normalised <- normalize_height(las = ctg_classified, algorithm = tin())
toc()

# Rasterize plots ####

# You can optimize processing by utilizing RAM better in 4 ways:
# 1. Use smaller chunk/plot sizes
# 2. Use opt_select() to load only needed fields into memory
# 3. Decrease amount of active workers (threads) as each workers uses own RAM
# 4. Exclude ground points and sub-surface noise

ctg_normalised <- readLAScatalog(paste0("E:/Remote Sensing Media/",date_folder,"/06. Point Clouds Normalised"))

# Limit the amount of workers (threads) if you don't have enough RAM. Each worker uses own RAM.
# plan(multisession)
plan(multisession, workers = 6)
opt_independent_files(ctg_normalised) <- TRUE
opt_select(ctg_normalised) <- "xyz"

# Drop ground points and sub-surface noise
opt_filter(ctg_normalised) <- "-drop_z_below 0 -drop_z_above 30"

tic()
# Write to disk rather than memory:
opt_output_files(ctg_normalised) <- paste0("E:/Remote Sensing Media/",date_folder,"/07. Canopy Height Models/", "{*}_chm")
# Rasterize canopy with interpolation:

ctg_chm <- rasterize_canopy(ctg_normalised,
                            res = 0.05,
                            algorithm = p2r(na.fill = tin()))
print("Rasterize canopy time:")
toc()

# Extract tree heights ####

tic()
# Calculate metrics using exact_extract (Outputs directly as a vector)
trees$Tree_Height <- exact_extract(ctg_chm, trees, 'max')

# Save to shapefile (completely overwriting old files to prevent schema errors)
st_write(trees, paste0("E:/Remote Sensing Media/",date_folder,"/09. Crown Metrics/Crown_Metrics_", file_date_safe, ".shp"), delete_dsn = TRUE)

# Save lightweight CSV without the messy spatial geometry text
write.csv(st_drop_geometry(trees), paste0("E:/Remote Sensing Media/",date_folder,"/09. Crown Metrics/Crown_Metrics_", file_date_safe, ".csv"), row.names = FALSE)
toc()