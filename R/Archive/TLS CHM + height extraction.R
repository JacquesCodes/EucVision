# ==============================================================================
# DYNAMIC TLS HEIGHT EXTRACTION SCRIPT (CORRECTED)
# ==============================================================================

# 1. Setup and Imports
library(lidR)            
library(sf)              
library(terra)           
library(dplyr)           
library(exactextractr)   
library(tictoc)          

# 2. Configuration & Paths
dataset_dir <- "E:/Remote Sensing Media/07. December 2025 (TLS)"
file_date_safe <- "December_2025_(TLS)"

normalised_dir <- file.path(dataset_dir, "06. Point Clouds Normalised")
chm_dir <- file.path(dataset_dir, "07. Canopy Height Models")
polygons_dir <- file.path(dataset_dir, "08. Crown Polygons")
metrics_dir <- file.path(dataset_dir, "09. Crown Metrics")

single_chm_path <- file.path(chm_dir, paste0("Master_Site_CHM_Single_", file_date_safe, ".tif"))
crown_shp_path <- file.path(polygons_dir, paste0("Crown_Polygons_", file_date_safe, ".shp"))

master_csv_path <- "C:/Users/jakev/Downloads/02. Processed Master Dataset.csv"

# 3. CHM Generation from Normalised Clouds
tic("CHM Rasterization complete")
print("Loading normalised catalog and applying general atmospheric filter...")

ctg_normalised <- readLAScatalog(normalised_dir)
opt_independent_files(ctg_normalised) <- TRUE
opt_select(ctg_normalised) <- "xyz"

# Removes extreme atmospheric outliers (birds/clouds high up) 
# Leaves the actual canopy fine-tuning to the dynamic polygon script below.
opt_filter(ctg_normalised) <- "-drop_z_below 0 -drop_z_above 15"
opt_output_files(ctg_normalised) <- file.path(chm_dir, "{*}_chm")

print("Generating clean, individual Canopy Height Models...")
ctg_chm <- rasterize_canopy(ctg_normalised, res = 0.05, algorithm = p2r(na.fill = tin()))
toc()

# 4. Spatial Data & Reference Data Loading
print("Loading polygons and dynamic height reference data...")

trees <- st_read(crown_shp_path, quiet = TRUE)
if (is.na(st_crs(trees)$epsg) || st_crs(trees)$epsg != 2048) {
  trees <- st_transform(trees, 2048)
}
trees <- trees %>% select(-any_of(c("group_ulid", "N_GM", "id", "N_FG", "N_BG", "BBox")))

# FIX 1: Correct the 10-character shapefile truncation
if ("Compartmen" %in% names(trees)) {
  names(trees)[names(trees) == "Compartmen"] <- "Compartment"
}

# FIX 2: Safely construct Tree_ID
if(!"Tree_ID" %in% names(trees)) {
  trees$Tree_ID <- paste(trees$Compartment, trees$Line, trees$Plot, trees$Tree, sep="_")
}

# Load the Processed Master Dataset for dynamic capping
ref_data <- read.csv(master_csv_path)

# Extract the maximum known Adjusted_Height for each tree to use as our baseline expected height
expected_heights <- ref_data %>%
  group_by(Tree_ID) %>%
  summarise(Expected_Height = max(Adjusted_Height, na.rm = TRUE)) %>%
  ungroup()

# FIX 3: Trim any invisible whitespaces from both datasets before joining
trees$Tree_ID <- trimws(as.character(trees$Tree_ID))
expected_heights$Tree_ID <- trimws(as.character(expected_heights$Tree_ID))

# Merge the expected heights into the spatial trees object
trees <- left_join(trees, expected_heights, by = "Tree_ID")

# FIX 4: Add a debug check so you know if the join fails again!
na_count <- sum(is.na(trees$Expected_Height))
print(paste("Number of trees missing historical data (defaulting to 3.0m fallback):", na_count, "out of", nrow(trees)))

# If there are any truly new/missing trees without historical data, give them the safe fallback
trees$Expected_Height[is.na(trees$Expected_Height)] <- 3.0 

# 5. Metric Extraction & Consolidation (DYNAMIC FILTERING)
tic("Metric extraction complete")
print("Extracting metrics using dynamic per-tree thresholds...")

# Gather the newly generated individual CHMs
chm_files <- list.files(chm_dir, pattern = "\\.tif$", full.names = TRUE)
# Exclude the old Master CHM if it still exists in the folder so it doesn't stitch itself
chm_files <- chm_files[!grepl("Master_Site_CHM_Single", chm_files)]

site_chm_vrt <- terra::vrt(chm_files)

# STEP A: Extract ALL pixels inside each polygon (returns a list of dataframes)
extracted_pixels <- exact_extract(site_chm_vrt, trees)

# STEP B: Dynamically filter pixels and find the true max
# We allow a 0.5m buffer above the expected historical height for recent growth
growth_buffer <- 0.5 

trees$Tree_Height <- sapply(seq_along(extracted_pixels), function(i) {
  # Get the raw pixels for this specific tree polygon
  pixels <- extracted_pixels[[i]]$value
  
  # Calculate the dynamic threshold for this specific tree
  max_allowed_height <- trees$Expected_Height[i] + growth_buffer
  
  # Drop NA values and drop noise pixels that exceed the threshold
  valid_pixels <- pixels[!is.na(pixels) & pixels <= max_allowed_height]
  
  # Return the highest valid pixel, or NA if no valid pixels exist
  if(length(valid_pixels) > 0) {
    return(max(valid_pixels))
  } else {
    return(NA) 
  }
})

# Write the final physical TIFF file (No more terra::clamp!)
terra::writeRaster(site_chm_vrt, filename = single_chm_path, overwrite = TRUE)

# Export Outputs
print("Exporting Shapefile and CSV...")
st_write(trees, file.path(metrics_dir, paste0("Crown_Metrics_", file_date_safe, ".shp")), delete_dsn = TRUE, quiet = TRUE)
write.csv(st_drop_geometry(trees), file.path(metrics_dir, paste0("Crown_Metrics_", file_date_safe, ".csv")), row.names = FALSE)

toc()
print("Process finished successfully.")