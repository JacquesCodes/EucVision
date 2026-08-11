# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: 5cm POINT DENSITY RASTER GENERATION (HIT MAPS)
# ──────────────────────────────────────────────────────────────────────────────
# Description: Generates high-resolution (5cm) 2D density rasters for visual 
#              inspection of SfM reconstruction gaps. Separates Canopy (>0.3m) 
#              from Ground (<=0.3m) so soil density doesn't mask canopy holes.
# ──────────────────────────────────────────────────────────────────────────────

library(lidR)
library(terra)
library(stringr)
library(sf) # Added for shapefile handling

# ──────────────────────────────────────────────────────────────────────────────
# 1. Configuration
# ──────────────────────────────────────────────────────────────────────────────
base_dir <- "E:/Remote Sensing Media"
raster_res <- 1    # 5cm pixel resolution
z_threshold <- 0.3    # Empty trunk space cut-off

# Added the DTM and Plot boundaries required for the stitching mask
baseline_dtm_path <- "E:/Remote Sensing Media/00. Baseline DTM/IMPACT_OAL_Baseline_DTM_ESRI_102562.tif"
plots_shp_path <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. QGIS Shapefiles/1. LAScatalog Plot Boundaries/LAScatalog_Plot_Boundaries_ESRI_102562.shp"

print("Loading static spatial data for masking...")
baseline_dtm <- rast(baseline_dtm_path)
plots <- st_read(plots_shp_path, quiet = TRUE)
st_crs(plots) <- st_crs(baseline_dtm) # Force CRS match

# Test on a single date first to check the outputs in QGIS
target_date_override <- "20. 23 March 2026"

exclude_list <- c("000. Projects", "00. Baseline DTM", "00. Dataset Template",
                  "01. 25 February 2025",
                  "02. 01 September 2025",
                  "07. December 2025 (TLS)",
                  "17. 02 March 2026",
                  "17. 02 March 2026 2.4",
                  "17. 02 March 2026 19.2",
                  "20. 23 March 2026 0.6cm")

folders <- list.dirs(base_dir, recursive = FALSE)
dataset_folders <- folders[grepl("^\\d{2}\\.", basename(folders)) & !basename(folders) %in% exclude_list]

if (!is.null(target_date_override)) {
  dataset_folders <- dataset_folders[basename(dataset_folders) == target_date_override]
}

# ──────────────────────────────────────────────────────────────────────────────
# 2. Batch Processing Loop
# ──────────────────────────────────────────────────────────────────────────────
for (folder_path in dataset_folders) {
  
  date_folder <- basename(folder_path)
  file_date <- sub("^\\d+\\.\\s*", "", date_folder)
  file_date_safe <- gsub(" ", "_", file_date)
  
  print(paste("================================================================"))
  print(paste("GENERATING DENSITY RASTERS:", date_folder))
  print(paste("================================================================"))
  
  # Pointing to the outputs from Step 6 of your master pipeline
  normalised_dir <- file.path(folder_path, "06. Point Clouds Normalised")
  out_dir <- file.path(folder_path, "11. Density Rasters")
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  las_files <- list.files(normalised_dir, pattern = "\\.(las|laz)$", full.names = TRUE)
  
  if (length(las_files) == 0) {
    print(paste("-> SKIPPED: No normalized point clouds found in", normalised_dir))
    next
  }
  
  for (file in las_files) {
    plot_id_match <- str_extract(basename(file), "Plot_\\d+")
    plot_id <- ifelse(is.na(plot_id_match), basename(file), plot_id_match)
    
    # Check if this specific plot raster already exists to save time
    canopy_raster_path <- file.path(out_dir, paste0(plot_id, "_Canopy_Density_5cm_", file_date_safe, ".tif"))
    if (file.exists(canopy_raster_path)) next
    
    las <- readLAS(file)
    if (is.empty(las)) next
    
    # ----------------------------------------------------------------------
    # Canopy Hit Map (> 0.3m)
    # ----------------------------------------------------------------------
    canopy_las <- filter_poi(las, Z > z_threshold)
    
    if (!is.empty(canopy_las)) {
      # rasterize_density calculates points per square meter by default.
      # We use it to create a 2D map of where the canopy points actually exist.
      canopy_density_rast <- rasterize_density(canopy_las, res = raster_res)
      
      # REMOVED the NA to 0 assignment here to fix the bounding box overlap
      terra::writeRaster(canopy_density_rast, filename = canopy_raster_path, overwrite = TRUE)
    }
    
    # ----------------------------------------------------------------------
    # Ground Hit Map (<= 0.3m) - Optional, but great for side-by-side contrast
    # ----------------------------------------------------------------------
    ground_las <- filter_poi(las, Z <= z_threshold)
    
    if (!is.empty(ground_las)) {
      ground_raster_path <- file.path(out_dir, paste0(plot_id, "_Ground_Density_5cm_", file_date_safe, ".tif"))
      ground_density_rast <- rasterize_density(ground_las, res = raster_res)
      
      # REMOVED the NA to 0 assignment here to fix the bounding box overlap
      terra::writeRaster(ground_density_rast, filename = ground_raster_path, overwrite = TRUE)
    }
    
    print(paste("   -> Saved 5cm Rasters for:", plot_id))
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # 3. Stitch Master Hit Maps (Site-Level Consolidation)
  # ────────────────────────────────────────────────────────────────────────────
  print("Stitching individual plots into master site-level Hit Maps...")
  
  # Fixed file names (removed undefined 'flight_type')
  master_canopy_path <- file.path(out_dir, paste0("Master_Site_Canopy_Density_5cm_", file_date_safe, ".tif"))
  master_ground_path <- file.path(out_dir, paste0("Master_Site_Ground_Density_5cm_", file_date_safe, ".tif"))
  
  # Fixed pattern (removed undefined 'flight_type')
  canopy_files <- list.files(out_dir, pattern = "^Plot_.*Canopy_Density_5cm_.*\\.tif$", full.names = TRUE)
  
  if (length(canopy_files) > 0) {
    # 1. Stitch transparently
    site_canopy_vrt <- terra::vrt(canopy_files)
    
    # 2. Create a precise 1/NA mask from your QGIS polygons
    plot_mask <- terra::rasterize(plots, site_canopy_vrt, field = 1)
    
    # 3. If inside a plot (mask == 1) AND the density is NA (a hole), make it 0. Otherwise keep original.
    # We write directly to disk during the ifel() command to save RAM.
    terra::ifel(plot_mask == 1 & is.na(site_canopy_vrt), 0, site_canopy_vrt, 
                filename = master_canopy_path, overwrite = TRUE)
    
    print("   -> Successfully saved: Master Canopy Hit Map (Masked)")
    
    # Clean up individual plot rasters to save space (optional)
    file.remove(canopy_files)
  }
  
  ground_files <- list.files(out_dir, pattern = "^Plot_.*Ground_Density_5cm_.*\\.tif$", full.names = TRUE)
  
  if (length(ground_files) > 0) {
    site_ground_vrt <- terra::vrt(ground_files)
    plot_mask_grnd <- terra::rasterize(plots, site_ground_vrt, field = 1)
    
    terra::ifel(plot_mask_grnd == 1 & is.na(site_ground_vrt), 0, site_ground_vrt, 
                filename = master_ground_path, overwrite = TRUE)
    
    print("   -> Successfully saved: Master Ground Hit Map (Masked)")
    
    # Clean up individual plot rasters to save space (optional)
    file.remove(ground_files)
  }
  
  # Clean up memory
  rm(las, canopy_las, ground_las, canopy_density_rast, ground_density_rast)
  gc()
}

print("================================================================")
print("DENSITY RASTER GENERATION COMPLETE!")
print("================================================================")