library(lidR)
library(sf)
library(tictoc)
library(dplyr)
library(terra)
library(exactextractr) 

# 1. Setup and Data Loading ####

# Set batch date folder 
date_folder <- "Benchmark"

# Read in tree location for height extraction
trees <- st_read(paste0("E:/Remote Sensing Media/", date_folder, "/08. Crown shape file/All_Plots.shp"))

# Clean up the shapefile by removing unnecessary columns 
# (This also removes the old 'id' column, fixing the duplicate name crash!)
trees <- trees %>% 
  select(-group_ulid, -N_GM, -id, -N_FG, -N_BG, -BBox)

# Read in Canopy height models (VRT)
ctg_chm <- rast(paste0("E:/Remote Sensing Media/", date_folder, "/07. Canopy Height Models/rasterize_canopy.vrt"))

# Create a fresh, unique ID column for joining
trees$ID <- 1:nrow(trees)

# ---------------------------------------------------------
# 2. Benchmark: terra::extract ####
# ---------------------------------------------------------
print("Starting Benchmark: terra...")
tic("terra_timer")

# Calculate metrics using terra
tree_heights_terra <- terra::extract(ctg_chm, trees, fun = max, na.rm = TRUE)

# Rename the extracted column to strictly 7 characters to satisfy ESRI Shapefile limits
colnames(tree_heights_terra)[2] <- "MaxHt_T"

# Join results back
trees_terra_final <- left_join(trees, tree_heights_terra, by = "ID")

# Save outputs as Shapefile (.shp), completely wiping the old files to prevent DBF schema errors
st_write(trees_terra_final, paste0("E:/Remote Sensing Media/", date_folder, "/09. Tree heights/All_Plots_Terra.shp"), delete_dsn = TRUE)
write.csv(st_drop_geometry(trees_terra_final), paste0("E:/Remote Sensing Media/", date_folder, "/09. Tree heights/All_Plots_Terra.csv"), row.names = FALSE)

terra_time <- toc(log = TRUE)

# ---------------------------------------------------------
# 3. Benchmark: exactextractr::exact_extract ####
# ---------------------------------------------------------
print("Starting Benchmark: exactextractr...")
tic("exact_timer")

# Calculate metrics using exact_extract and assign a short 7-character column name
trees$MaxHt_E <- exact_extract(ctg_chm, trees, 'max')

# Save outputs as Shapefile (.shp)
st_write(trees, paste0("E:/Remote Sensing Media/", date_folder, "/09. Tree heights/All_Plots_Exact.shp"), delete_dsn = TRUE)
write.csv(st_drop_geometry(trees), paste0("E:/Remote Sensing Media/", date_folder, "/09. Tree heights/All_Plots_Exact.csv"), row.names = FALSE)

exact_time <- toc(log = TRUE)

# ---------------------------------------------------------
# 4. Final Comparison Results ####
# ---------------------------------------------------------
cat("\n--- BENCHMARK RESULTS ---\n")
print(tic.log(format = TRUE))