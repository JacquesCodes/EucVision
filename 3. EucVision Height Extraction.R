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

# Read in already processed Canopy height models ####

# Change this single variable for each new batch!
date_folder <- "23. 13 April 2026"

# Extract the date part for file naming
file_date <- sub("^\\d+\\.\\s*", "", date_folder)
file_date_safe <- gsub(" ", "_", file_date)

# 1. Read in tree locations (Updated folder and file name)
path_trees <- paste0("E:/Remote Sensing Media/", date_folder, "/08. Crown Polygons/Crown_Polygons_", file_date_safe, ".shp")
trees <- st_read(path_trees)

# Automatically check and transform to EPSG: 2048 if it doesn't match
if (is.na(st_crs(trees)$epsg) || st_crs(trees)$epsg != 2048) {
  trees <- st_transform(trees, 2048)
  message("Transformed CRS to 2048 successfully.")
}

# Safely remove old columns if they exist to prevent duplicate name/ESRI errors
trees <- trees %>%
  select(-any_of(c("grop_ld", "group_ulid", "N_GM", "id", "N_FG", "N_BG", "BBox")))

# 2. Read in Canopy height models (Updated folder name to Title Case)
# Keeping the .vrt reference as requested
ctg_chm <- rast(paste0("E:/Remote Sensing Media/", date_folder, "/07. Canopy Height Models/rasterize_canopy.vrt"))

# 3. Define the color palette based on your species
species_colors <- c(
  "Cladocalyx"    = "#336998",
  "Grandis"       = "#97dde3",
  "Cloeziana"     = "#ffffff",
  "Urophylla"     = "#e3acff",
  "Grandis clone" = "#ff7da0"
)

# 4. Plot the Canopy Height Model
plot(ctg_chm, range = c(0, 5), main = "Canopy Height Model with Tree Species")

# 5. Plot the tree shapefiles overlaid, mapping the 'Species' column to the colors
# Note: 'border = "black"' adds a black outline so the white polygons don't blend in
plot(trees$geometry, 
     add = TRUE, 
     col = species_colors[trees$Species])

# Extract tree heights ####

tic()
# 6. Calculate metrics using exact_extract (Outputs directly as a vector)
trees$Tree_Height <- exact_extract(ctg_chm, trees, 'max')

# 7.Save to Crown Metrics (Updated folder and file name)
output_path_shp <- paste0("E:/Remote Sensing Media/", date_folder, "/09. Crown Metrics/Crown_Metrics_", file_date_safe, ".shp")
output_path_csv <- paste0("E:/Remote Sensing Media/", date_folder, "/09. Crown Metrics/Crown_Metrics_", file_date_safe, ".csv")

# Save to shapefile (completely overwriting old files to prevent schema errors)
st_write(trees, output_path_shp, delete_dsn = TRUE)

# Save lightweight CSV without the messy spatial geometry text
write.csv(st_drop_geometry(trees), output_path_csv, row.names = FALSE)

message("✅ Extraction complete. Metrics saved to: ", output_path_csv)
toc()