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

# Read in all point clouds and shape files ####

# My path to the remote sensing dataset
myPath <- "E:/Remote Sensing Media/03. Feb 2025/"

# Read in point clouds into a catalog (ctg)
ctg <- readLAScatalog(paste0(myPath,"03. Point clouds"))

# Read in shape files for points of interest boundaries like my plots
plots_buffered_unsorted <- st_read(paste0("E:/Remote Sensing Media/00. QGIS Projects/1. Plot boundaries for cropping/R Plots.shp"))
plots <- plots_buffered_unsorted[order(plots_buffered_unsorted$id), ]

# Read in tree location for height extraction at the end
trees <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/QGIS Combined Output/All_Plots.shp")

# View catalog, plots and trees
plot(ctg)
plot(plots$geometry,add = TRUE)
plot(trees$geometry, add = TRUE, col = "red")

# Crop plots ####

# Multiple CPU threads mode
plan(multisession)
# .las files overlaps and are not independent
opt_independent_files(ctg) <- FALSE
# Only load needed variables
opt_select(ctg) <- "xyz"

# tic() & toc() is to check code running time
tic()
# Write to disk rather than memory:
opt_output_files(ctg) <- paste0(myPath,"04. Clipped/", "Plot_{ID}")
# Crop plots
ctg_clipped <- clip_roi(ctg, plots)
toc()

# Classify plots ####

plan(multisession)
opt_independent_files(ctg_clipped) <- TRUE
opt_select(ctg_clipped) <- "xyz"

tic()
# Write to disk rather than memory:
opt_output_files(ctg_clipped) <- paste0(myPath,"05. Ground classified/", "{*}_classified")
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
# The cloth falls between the points and classify trees.

# Normalize plots ####

# If in the future I am unable to classify ground well. I should use a baseline raster instead:
# ctg_normalised <- normalize_height(las = ctg_classified, algorithm = tin(), dtm = )

plan(multisession)
opt_independent_files(ctg_classified) <- TRUE
opt_select(ctg_classified) <- "xyz"

tic()
# Write to disk rather than memory:
opt_output_files(ctg_classified) <- paste0(myPath,"06. Normalised/", "{*}_normalised")
# A point cloud-based normalization without a raster:
ctg_normalised <- normalize_height(las = ctg_classified, algorithm = tin())
toc()

# Rasterize plots ####

# You can optimize processing by utilizing RAM better in 3 ways:
# 1. Use smaller chunk/plot sizes
# 2. Use opt_select() to load only needed fields into memory
# 3. Decrease amount of active workers (threads) as each workers uses own RAM

# Limit the amount of workers (threads) if you don't have enough RAM. Each worker uses own RAM.
# plan(multisession)
plan(multisession, workers = 6)
opt_independent_files(ctg_normalised) <- TRUE
opt_select(ctg_normalised) <- "xyz"

tic()
# Write to disk rather than memory:
opt_output_files(ctg_normalised) <- paste0(myPath,"07. Canopy Height Model/", "{*}_chm")
# Rasterize canopy with interpolation:
ctg_chm <- rasterize_canopy(ctg_normalised, res = 0.01, algorithm = p2r(na.fill = tin()))
print("Rasterize canopy time:")
toc()

# Extract trees' heights ####

tic()
# Ensure both have an ID column
trees$ID <- 1:nrow(trees)

# Calculate metrics
tree_heights <- terra::extract(ctg_chm, trees, fun = max, na.rm = TRUE)

# Join results back using the ID
trees_with_heights <- left_join(trees, st_drop_geometry(tree_heights), by = "ID")

# Save to file
st_write(trees_with_heights, paste0(myPath,"08. Heights/All Plots.shp"))
st_write(trees_with_heights, paste0(myPath,"08. Heights/All Plots.csv"))
toc()









