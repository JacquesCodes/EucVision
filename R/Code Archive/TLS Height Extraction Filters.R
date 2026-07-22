# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: DYNAMIC TLS HEIGHT EXTRACTION & FILTERING PIPELINE ####
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────
# Description: Automates the extraction of maximum tree heights from plot-level 
#              Canopy Height Models (CHMs) using dynamically calculated, per-tree 
#              height thresholds. It cross-references historical baseline data 
#              to establish expected tree heights, applies a strict physiological 
#              growth buffer constraint to eliminate canopy noise or anomalies, 
#              and extracts the true maximum pixel value within each mapped 
#              crown polygon. Output metrics are rigorously formatted and 
#              exported as both spatial (Shapefile) and tabular (CSV) datasets.
# ──────────────────────────────────────────────────────────────────────────────

# 1. Setup and Imports
library(sf)              
library(terra)           
library(dplyr)           
library(exactextractr)   

# 2. Configuration & Paths
dataset_dir <- "E:/Remote Sensing Media/07. December 2025 (TLS)"

date_folder <- basename(dataset_dir)

# Extract the date part and create a safe filename format
# (e.g., "17. 02 March 2026" -> "02_March_2026")
file_date_safe <- gsub(" ", "_", sub("^\\d+\\.\\s*", "", date_folder))

chm_dir <- file.path(dataset_dir, "07. Canopy Height Models")
polygons_dir <- file.path(dataset_dir, "08. Crown Polygons")
metrics_dir <- file.path(dataset_dir, "09. Crown Metrics")

single_chm_path <- file.path(chm_dir, paste0("Master_Site_CHM_Single_", file_date_safe, ".tif"))
crown_shp_path <- file.path(polygons_dir, paste0("Crown_Polygons_", file_date_safe, ".shp"))

master_csv_path <- "C:/Users/jakev/Downloads/02. Processed Master Dataset.csv"

# 3. Spatial Data & Reference Data Loading
print("Loading polygons and dynamic height reference data...")

trees <- st_read(crown_shp_path, quiet = TRUE)
if (is.na(st_crs(trees)$epsg) || st_crs(trees)$epsg != 2048) {
  trees <- st_transform(trees, 2048)
}
trees <- trees %>% select(-any_of(c("group_ulid", "N_GM", "id", "N_FG", "N_BG", "BBox")))

# FIX 1: Account for the exact shapefile abbreviation "Cmprtmn" or "Compartmen"
if ("Cmprtmn" %in% names(trees)) {
  names(trees)[names(trees) == "Cmprtmn"] <- "Compartment"
} else if ("Compartmen" %in% names(trees)) {
  names(trees)[names(trees) == "Compartmen"] <- "Compartment"
}

# FIX 2: Safely construct Tree_ID matching the strict Python 2-decimal format
if(!"Tree_ID" %in% names(trees)) {
  # Forces tree numbers like 3.1 to become exactly "3.10"
  formatted_tree <- sprintf("%.2f", as.numeric(trees$Tree))
  trees$Tree_ID <- paste(trees$Compartment, trees$Line, trees$Plot, formatted_tree, sep="_")
}

# Load the Processed Master Dataset
ref_data <- read.csv(master_csv_path)

# FIX 3: STRICT BASELINE FILTER - Only use adjusted heights from November 28, 2025
ref_data_filtered <- ref_data %>%
  filter(Date == "2025-11-28")

# Extract the known Adjusted_Height for each tree from that specific date
expected_heights <- ref_data_filtered %>%
  group_by(Tree_ID) %>%
  summarise(
    # Safely handle completely blank trees to avoid -Inf errors
    Expected_Height = suppressWarnings(
      if(all(is.na(Adjusted_Height))) NA_real_ else max(Adjusted_Height, na.rm = TRUE)
    )
  ) %>%
  ungroup()

# Trim any invisible whitespaces from both datasets before joining
trees$Tree_ID <- trimws(as.character(trees$Tree_ID))
expected_heights$Tree_ID <- trimws(as.character(expected_heights$Tree_ID))

# Merge the expected heights into the spatial trees object
trees <- left_join(trees, expected_heights, by = "Tree_ID")

# Debug check: Validate how many of the 3,144 trees joined successfully
na_count <- sum(is.na(trees$Expected_Height))
print(paste("Number of trees missing historical Nov 28 data (defaulting to 3.0m fallback):", na_count, "out of", nrow(trees)))

# Apply the safe fallback only to trees truly missing from the Nov 28 dataset
trees$Expected_Height[is.na(trees$Expected_Height)] <- 3.0 

# 4. Metric Extraction & Consolidation (DYNAMIC FILTERING)
print("Extracting metrics using dynamic per-tree thresholds...")

# Gather the previously generated individual CHMs
chm_files <- list.files(chm_dir, pattern = "\\.tif$", full.names = TRUE)
chm_files <- chm_files[!grepl("Master_Site_CHM_Single", chm_files)]

# Create a virtual raster from the individual CHMs
site_chm_vrt <- terra::vrt(chm_files)

# STEP A: Extract ALL pixels inside each polygon
extracted_pixels <- exact_extract(site_chm_vrt, trees)

# STEP B: Dynamically filter pixels and find the true max
# Set the rigid 0.5m growth buffer above the Nov 28 baseline
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

# Optional: Write the final physical TIFF file for the entire site
terra::writeRaster(site_chm_vrt, filename = single_chm_path, overwrite = TRUE)

# Export Outputs
print("Exporting Shapefile and CSV...")
st_write(trees, file.path(metrics_dir, paste0("Crown_Metrics_", file_date_safe, ".shp")), delete_dsn = TRUE, quiet = TRUE)
write.csv(st_drop_geometry(trees), file.path(metrics_dir, paste0("Crown_Metrics_", file_date_safe, ".csv")), row.names = FALSE)

print("Extraction finished successfully.")