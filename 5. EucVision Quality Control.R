library(lidR)
library(RCSF)
library(RMCC)
library(gstat)
library(sf)
library(tictoc)
library(geometry)
library(dplyr)
library(future)
library(terra)
library(rgl)

# Change this single variable for each new batch!
date_folder <- "19. 16 March 2026"

# Plot number
Number <- 37

# Extract the date part by removing the leading folder number, dot, and space
file_date <- sub("^\\d+\\.\\s*", "", date_folder)

# Replace spaces with underscores for safer file naming conventions
file_date_safe <- gsub(" ", "_", file_date)

# My path to the remote sensing dataset
myPath <- paste0("E:/Remote Sensing Media/",date_folder,"/")

# Dynamically construct the single-date file names
name_clipped <- paste0("Plot_", Number, "_", file_date_safe)
name_classified <- paste0(name_clipped, "_classified")
name_normalised <- paste0(name_classified, "_normalised")
name_chm <- paste0(name_normalised, "_chm")

# Build the full file paths
path_clipped <- paste0(myPath, "04. Point Clouds Clipped/", name_clipped, ".las")
path_classified <- paste0(myPath, "05. Point Clouds Ground Classified/", name_classified, ".las")
path_normalised <- paste0(myPath, "06. Point Clouds Normalised/", name_normalised, ".las")
path_chm <- paste0(myPath, "07. Canopy Height Models/", name_chm, ".tif")

# --- Safely Read in the Data ---

if (file.exists(path_clipped)) {
  las <- readLAS(path_clipped)
  message("Loaded: ", name_clipped, ".las")
} else {
  message("Skipped: Could not find folder or file for 04. Point Clouds Clipped")
  las <- NULL
}

if (file.exists(path_classified)) {
  las_classified <- readLAS(path_classified)
  message("Loaded: ", name_classified, ".las")
} else {
  message("Skipped: Could not find folder or file for 05. Point Clouds Ground Classified")
  las_classified <- NULL
}

if (file.exists(path_normalised)) {
  las_normalised <- readLAS(path_normalised)
  message("Loaded: ", name_normalised, ".las")
} else {
  message("Skipped: Could not find folder or file for 06. Point Clouds Normalised")
  las_normalised <- NULL
}

if (file.exists(path_chm)) {
  las_chm <- rast(path_chm)
  
  # Filter out the rogue photogrammetry noise
  las_chm[las_chm > 15] <- NA  # Replace 8 with your realistic upper threshold
  las_chm[las_chm < 0] <- 0   # Optional: Cleans up any sub-surface negative noise
  
  message("Loaded: ", name_chm, ".tif")
} else {
  message("Skipped: Could not find folder or file for 07. Canopy Height Models")
  las_chm <- NULL
}

# Safely load and transform the Crown Polygons shapefile
path_trees <- paste0(myPath, "08. Crown Polygons/Crown_Polygons_", file_date_safe, ".shp")

if (file.exists(path_trees)) {
  trees <- st_read(path_trees)
  
  # Crucial CRS fix to prevent "spatial index out of range" errors!
  if (is.na(st_crs(trees)$epsg) || st_crs(trees)$epsg != 2048) {
    trees <- st_transform(trees, 2048)
    message("  -> Transformed Crown Polygons CRS to EPSG:2048")
  }
  
  # Filter for just the plot we are currently visualizing
  PlotTrees <- trees[trees$Plot == paste0("Plot_", Number),]
  message("Loaded and filtered Crown Polygons for Plot ", Number)
  
} else {
  message("Skipped: Could not find folder or file for 08. Crown Polygons")
  PlotTrees <- NULL
}

# --- Plotting Section ---

# Plot cropped las
if (!is.null(las)) {
  plot(las, size = 4, bg = "white")
}

# Plot classified las and DTM
if (!is.null(las_classified)) {
  gnd <- filter_ground(las_classified)
  plot(gnd, size = 4, bg = "white")
  
  dtm_tin_0 <- rasterize_terrain(las_classified, res = 1, algorithm = tin())
  plot_dtm3d(dtm_tin_0, bg = "white")
}

# Plot normalised las
if (!is.null(las_normalised)) {
  plot(las_normalised, bg = "white")
}

# Plot canopy height model
if (!is.null(las_chm)) {
  plot(las_chm, col = height.colors(50))
}


# --- 8.1 Individual Tree Detection (ITD) ---

MinimumTreeHeight <- 0.5

# create Local Maximum Filter (lmf) function for the "ws" search
lmf_algorithm <- lmf(ws = 3, hmin = MinimumTreeHeight, shape = "circular")

# Only run tree detection if the CHM exists
if (!is.null(las_chm)) {
  # Locate trees in a circle with a diameter of "ws" in meters
  ttops <- locate_trees(las = las_chm, algorithm = lmf_algorithm)
  
  # Tree detection results in 2D
  plot(las_chm, col = height.colors(50))
  plot(sf::st_geometry(ttops), add = TRUE, pch = 3)
  
  # Tree detection results can also be visualized in 3D!
  # This requires BOTH the normalised point cloud and the CHM to exist
  if (!is.null(las_normalised)) {
    x <- plot(las_normalised, bg = "white", size = 2)
    add_treetops3d(x, ttops, radius = 0.15, fastTransparency = TRUE, alpha = 0.8)
  } else {
    message("Cannot plot 3D treetops: Normalised point cloud is missing.")
  }
}

# # Video
# 
# # --- Define spin motion ---
# spin <- spin3d(axis = c(0, 0, 1), rpm = 6)
# 
# # Doesn't work if you make the plot fullscreen.
# # --- Save the animation ---
# movie3d(
#   movie = "Tree_Tops_animation",   # base name for output
#   dir = getwd(),                   # output directory
#   spin,                            # the animation function
#   duration = 10,                   # 10 seconds
#   fps = 25,                        # frames per second
#   clean = TRUE,                    # remove frames after combining
#   type = "gif"                     # try making a .gif directly
# )