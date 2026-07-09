# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: SITE-WIDE COVERAGE & DENSITY METRICS (0.6cm vs 3cm)
# ──────────────────────────────────────────────────────────────────────────────
# Description: Generates permanent 0.25m Hit Count rasters for Ground and Canopy.
#              Extracts total Plot Area, Point Counts, and true Planar Densities.
#              Applies a strict Boolean overlay to calculate sub-1cm % Coverage.
#
# UPDATED (v3):
#  1. A date is processed as long as EITHER 3cm or 0.6cm point clouds exist.
#  2. Fallback logic correctly identifies and routes historical layout files vs 
#     modern layout files based on date.
#  3. single-catalog fallback added for combined Top+Bottom files.
#  4. Has_3cm_Flight / Has_06cm_Flight columns track data availability in the CSV.
#  5. Exclude list added to skip specific folders.
#  6. Target Date Override added for rapid single-folder testing.
# ──────────────────────────────────────────────────────────────────────────────

library(lidR)
library(sf)
library(terra)
library(exactextractr)
library(dplyr)
library(stringr)
library(future)

# ──────────────────────────────────────────────────────────────────────────────
# 1. Configuration & Spatial Data
# ──────────────────────────────────────────────────────────────────────────────
base_dir <- "E:/Remote Sensing Media"
baseline_dtm_path <- "E:/Remote Sensing Media/00. Baseline DTM/Compartments_DTM_from_UAV_ESRI_102562.tif"
plots_shp_path <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. QGIS Shapefiles/1. LAScatalog Plot Boundaries/LAScatalog_Plot_Boundaries_ESRI_102562.shp"
impact_bounds_path <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. QGIS Shapefiles/4. IMPACT Plot & Compartment Boundaries/IMPACT_Plot_&_Compartment_Boundaries_ESRI_102562.shp"

raster_res <- 0.25    # 0.25m pixels (0.0625 m^2 per pixel)
z_threshold <- 0.3    

# Exception Dates
gsd_3cm_only_exception_dates <- c("16 March 2026")
modern_layout_start_date <- as.Date("2026-03-23")

# Folder Exclusion List
exclude_list <- c("000. Projects", "00. Baseline DTM", "00. Dataset Template", 
                  "01. 25 February 2025","02. 01 September 2025", "07. December 2025 (TLS)",
                  "17. 02 March 2026","17. 02 March 2026 2.4","17. 02 March 2026 19.2",
                  "20. 23 March 2026 0.6cm","30. 30 June 2026 (ALS)")

# ADDED: Target Date Override for testing
# Set to a specific folder name (e.g., "16. 16 March 2026") to test, or NULL for full batch
target_date_override <- NULL

print("Loading static spatial data...")
baseline_dtm <- rast(baseline_dtm_path)

plots <- st_read(plots_shp_path, quiet = TRUE)
st_crs(plots) <- st_crs(baseline_dtm)
plots <- plots[order(plots$id), ]

impact_bounds <- st_read(impact_bounds_path, quiet = TRUE)
st_crs(impact_bounds) <- st_crs(baseline_dtm)

top_bound <- impact_bounds[grepl("TOP", impact_bounds$Name, ignore.case = TRUE), ]
bot_bound <- impact_bounds[grepl("BOTTOM", impact_bounds$Name, ignore.case = TRUE), ]

# Combined boundary used as the clip ROI when a single point cloud covers both
combined_bound <- st_union(st_geometry(top_bound), st_geometry(bot_bound))
combined_bound <- st_sf(geometry = combined_bound)
st_crs(combined_bound) <- st_crs(baseline_dtm)

plan(multisession, workers = 6)

# ──────────────────────────────────────────────────────────────────────────────
# 2. Master Folder List & Master Dual-Flight Loop
# ──────────────────────────────────────────────────────────────────────────────
folders <- list.dirs(base_dir, recursive = FALSE)

# Apply regex match AND filter out folders in the exclude_list
dataset_folders <- folders[grepl("^\\d{2}\\.", basename(folders)) & !basename(folders) %in% exclude_list]

# Apply Target Date Override if not NULL
if (!is.null(target_date_override)) {
  dataset_folders <- dataset_folders[basename(dataset_folders) == target_date_override]
  print("================================================================")
  print(paste("⚠️ TARGET DATE OVERRIDE ACTIVE: Testing only ->", target_date_override))
  print("================================================================")
}

for (folder_path in dataset_folders) {
  
  las_folder_std <- file.path(folder_path, "03. Point Clouds")
  las_folder_fine <- file.path(las_folder_std, "0.6 GSD")
  
  if (!dir.exists(las_folder_std)) next
  
  date_folder <- basename(folder_path)
  file_date <- sub("^\\d+\\.\\s*", "", date_folder)
  file_date_safe <- gsub(" ", "_", file_date)
  
  print(paste("================================================================"))
  print(paste("EXTRACTING SITE-WIDE METRICS:", date_folder))
  print(paste("================================================================"))
  
  out_dir <- file.path(folder_path, "12. Coverage Analysis")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  out_csv_path <- file.path(out_dir, paste0("Planar_Coverage_and_Density_", file_date_safe, ".csv"))
  if (file.exists(out_csv_path)) {
    print("-> SKIPPED: Analysis CSV already exists.")
    next
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # 3. Helper Function: Process Point Clouds into Permanent Hit Count Rasters
  # ────────────────────────────────────────────────────────────────────────────
  generate_site_rasters <- function(las_files, name_tag) {
    
    if (length(las_files) == 0) {
      return(list(ground = NULL, canopy = NULL))
    }
    
    top_files <- las_files[grepl("Top", basename(las_files), ignore.case = TRUE)]
    bot_files <- las_files[grepl("Bottom", basename(las_files), ignore.case = TRUE)]
    
    merged_grnd_path <- file.path(out_dir, paste0("Master_Ground_HitCount_0.25m_", name_tag, "_", file_date_safe, ".tif"))
    merged_can_path  <- file.path(out_dir, paste0("Master_Canopy_HitCount_0.25m_", name_tag, "_", file_date_safe, ".tif"))
    
    grnd_to_merge <- list()
    can_to_merge <- list()
    
    if (length(top_files) > 0 && length(bot_files) > 0) {
      
      print(paste0("   [", name_tag, "] Splitting processing by Compartment (Top/Bottom)..."))
      
      # --- Process Top ---
      ctg_top <- readLAScatalog(top_files)
      st_crs(ctg_top) <- st_crs(baseline_dtm) 
      opt_select(ctg_top) <- "xyz"
      opt_chunk_buffer(ctg_top) <- 0
      
      top_las <- clip_roi(ctg_top, top_bound)
      if (!is.empty(top_las)) {
        top_las <- normalize_height(top_las, baseline_dtm)
        
        top_grnd <- filter_poi(top_las, Z <= z_threshold)
        if (!is.empty(top_grnd)) grnd_to_merge[[length(grnd_to_merge) + 1]] <- pixel_metrics(top_grnd, ~length(Z), res = raster_res)
        
        top_can <- filter_poi(top_las, Z > z_threshold)
        if (!is.empty(top_can)) can_to_merge[[length(can_to_merge) + 1]] <- pixel_metrics(top_can, ~length(Z), res = raster_res)
      }
      
      # --- Process Bottom ---
      ctg_bot <- readLAScatalog(bot_files)
      st_crs(ctg_bot) <- st_crs(baseline_dtm)
      opt_select(ctg_bot) <- "xyz"
      opt_chunk_buffer(ctg_bot) <- 0
      
      bot_las <- clip_roi(ctg_bot, bot_bound)
      if (!is.empty(bot_las)) {
        bot_las <- normalize_height(bot_las, baseline_dtm)
        
        bot_grnd <- filter_poi(bot_las, Z <= z_threshold)
        if (!is.empty(bot_grnd)) grnd_to_merge[[length(grnd_to_merge) + 1]] <- pixel_metrics(bot_grnd, ~length(Z), res = raster_res)
        
        bot_can <- filter_poi(bot_las, Z > z_threshold)
        if (!is.empty(bot_can)) can_to_merge[[length(can_to_merge) + 1]] <- pixel_metrics(bot_can, ~length(Z), res = raster_res)
      }
      
    } else {
      
      print(paste0("   [", name_tag, "] Processing single catalog entirely (Top+Bottom combined)..."))
      
      ctg_all <- readLAScatalog(las_files)
      st_crs(ctg_all) <- st_crs(baseline_dtm)
      opt_select(ctg_all) <- "xyz"
      opt_chunk_buffer(ctg_all) <- 0
      
      all_las <- clip_roi(ctg_all, combined_bound)
      if (!is.empty(all_las)) {
        all_las <- normalize_height(all_las, baseline_dtm)
        
        all_grnd <- filter_poi(all_las, Z <= z_threshold)
        if (!is.empty(all_grnd)) grnd_to_merge[[length(grnd_to_merge) + 1]] <- pixel_metrics(all_grnd, ~length(Z), res = raster_res)
        
        all_can <- filter_poi(all_las, Z > z_threshold)
        if (!is.empty(all_can)) can_to_merge[[length(can_to_merge) + 1]] <- pixel_metrics(all_can, ~length(Z), res = raster_res)
      }
    }
    
    # --- Merge and Save ---
    results <- list(ground = NULL, canopy = NULL)
    
    if (length(grnd_to_merge) > 0) {
      site_grnd <- terra::mosaic(sprc(grnd_to_merge), fun = "max")
      crs(site_grnd) <- crs(baseline_dtm)
      terra::writeRaster(site_grnd, merged_grnd_path, overwrite = TRUE)
      results$ground <- site_grnd
    }
    
    if (length(can_to_merge) > 0) {
      site_can <- terra::mosaic(sprc(can_to_merge), fun = "max")
      crs(site_can) <- crs(baseline_dtm)
      terra::writeRaster(site_can, merged_can_path, overwrite = TRUE)
      results$canopy <- site_can
    }
    
    return(results)
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # 4. Classify files into 3cm / 0.6cm buckets, then generate rasters
  # ────────────────────────────────────────────────────────────────────────────
  main_files          <- list.files(las_folder_std, pattern = "\\.(las|laz)$", full.names = TRUE, recursive = FALSE)
  fine_subfolder_files <- list.files(las_folder_fine, pattern = "\\.(las|laz)$", full.names = TRUE, recursive = FALSE)
  
  if (length(fine_subfolder_files) > 0) {
    std_files  <- main_files
    fine_files <- fine_subfolder_files
    
  } else if (length(main_files) > 0) {
    parsed_date <- suppressWarnings(as.Date(trimws(file_date), format = "%d %B %Y"))
    is_modern_date <- !is.na(parsed_date) && parsed_date >= modern_layout_start_date
    
    if (trimws(file_date) %in% gsd_3cm_only_exception_dates) {
      print(paste0("-> Main-folder-only files classified as 3cm GSD (", file_date, " exception)."))
      std_files  <- main_files
      fine_files <- character(0)
      
    } else if (is_modern_date) {
      print(paste0("-> WARNING: '", file_date, "' is on/after ", modern_layout_start_date,
                   " but has no '0.6 GSD' subfolder. Defaulting main-folder files to 3cm GSD - please verify."))
      std_files  <- main_files
      fine_files <- character(0)
      
    } else {
      print(paste0("-> Main-folder-only files classified as 0.6cm GSD (historical pre-23-March-2026 layout)."))
      std_files  <- character(0)
      fine_files <- main_files
    }
    
  } else {
    std_files  <- character(0)
    fine_files <- character(0)
  }
  
  if (length(std_files) == 0 && length(fine_files) == 0) {
    print("-> SKIPPED: No usable point cloud files found for this date.")
    next
  }
  
  rasters_3cm  <- generate_site_rasters(std_files, "3cm")
  rasters_06cm <- generate_site_rasters(fine_files, "06cm")
  
  has_3cm  <- !is.null(rasters_3cm$canopy)
  has_06cm <- !is.null(rasters_06cm$canopy)
  
  if (!has_3cm && !has_06cm) {
    print("Error generating rasters (no usable ground/canopy points found). Skipping.")
    next
  }
  
  if (has_3cm && has_06cm) {
    print("-> Both 3cm and 0.6cm data available - running full coverage comparison.")
  } else if (has_3cm) {
    print("-> Only 3cm data available for this date - Coverage_Percent will be NA.")
  } else {
    print("-> Only 0.6cm data available for this date - Coverage_Percent will be NA.")
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # 5. The Mathematical Extraction
  # ────────────────────────────────────────────────────────────────────────────
  print("-> Extracting Plot Metrics and calculating Coverage...")
  
  if (has_3cm) {
    plots$PointSum_3cm_Grnd <- exact_extract(rasters_3cm$ground, plots, 'sum')
    plots$PointSum_3cm_Can  <- exact_extract(rasters_3cm$canopy, plots, 'sum')
    
    mask_3cm_grnd <- terra::ifel(rasters_3cm$ground > 0, 1, NA)
    mask_3cm_can  <- terra::ifel(rasters_3cm$canopy > 0, 1, NA)
    
    plots$Pix_3cm_Grnd <- exact_extract(mask_3cm_grnd, plots, 'sum')
    plots$Pix_3cm_Can  <- exact_extract(mask_3cm_can, plots, 'sum')
  }
  
  if (has_06cm) {
    plots$PointSum_06cm_Grnd <- exact_extract(rasters_06cm$ground, plots, 'sum')
    
    if (has_3cm) {
      can_06cm_aligned <- terra::resample(rasters_06cm$canopy, rasters_3cm$canopy, method = "near")
    } else {
      can_06cm_aligned <- rasters_06cm$canopy
    }
    
    plots$PointSum_06cm_Can <- exact_extract(can_06cm_aligned, plots, 'sum')
    
    mask_06cm_grnd <- terra::ifel(rasters_06cm$ground > 0, 1, NA)
    mask_06cm_can  <- terra::ifel(can_06cm_aligned > 0, 1, NA)
    
    plots$Pix_06cm_Grnd <- exact_extract(mask_06cm_grnd, plots, 'sum')
    plots$Pix_06cm_Can  <- exact_extract(mask_06cm_can, plots, 'sum')
    
    if (has_3cm) {
      strict_mask_06cm_can <- terra::ifel(mask_3cm_can == 1 & mask_06cm_can == 1, 1, NA)
      plots$Pix_06cm_Can_Strict <- exact_extract(strict_mask_06cm_can, plots, 'sum')
    }
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # 6. Structuring the CSV Math
  # ────────────────────────────────────────────────────────────────────────────
  pixel_area_m2 <- raster_res^2 # 0.0625 sq meters
  
  select_cols <- "id"
  if (has_3cm)  select_cols <- c(select_cols, "PointSum_3cm_Grnd", "PointSum_3cm_Can", "Pix_3cm_Grnd", "Pix_3cm_Can")
  if (has_06cm) select_cols <- c(select_cols, "PointSum_06cm_Grnd", "PointSum_06cm_Can", "Pix_06cm_Grnd", "Pix_06cm_Can")
  if (has_3cm && has_06cm) select_cols <- c(select_cols, "Pix_06cm_Can_Strict")
  
  results_df <- st_drop_geometry(plots) %>%
    select(all_of(select_cols)) %>%
    rename(Plot_ID = id) %>%
    mutate(
      Plot_ID = paste0("Plot_", Plot_ID),
      Flight_Date = file_date,
      Has_3cm_Flight = has_3cm,
      Has_06cm_Flight = has_06cm
    )
  
  # --- 3cm Metrics ---
  if (has_3cm) {
    results_df <- results_df %>%
      mutate(
        Ground_Area_3cm_m2 = round(coalesce(Pix_3cm_Grnd, 0) * pixel_area_m2, 2),
        Ground_Points_3cm = round(coalesce(PointSum_3cm_Grnd, 0)),
        Ground_Density_3cm = ifelse(Ground_Area_3cm_m2 > 0, round(Ground_Points_3cm / Ground_Area_3cm_m2, 2), 0),
        
        Canopy_Area_3cm_m2 = round(coalesce(Pix_3cm_Can, 0) * pixel_area_m2, 2),
        Canopy_Points_3cm = round(coalesce(PointSum_3cm_Can, 0)),
        Canopy_Density_3cm = ifelse(Canopy_Area_3cm_m2 > 0, round(Canopy_Points_3cm / Canopy_Area_3cm_m2, 2), 0)
      )
  } else {
    results_df <- results_df %>%
      mutate(
        Ground_Area_3cm_m2 = NA_real_, Ground_Points_3cm = NA_real_, Ground_Density_3cm = NA_real_,
        Canopy_Area_3cm_m2 = NA_real_, Canopy_Points_3cm = NA_real_, Canopy_Density_3cm = NA_real_
      )
  }
  
  # --- 0.6cm Metrics ---
  if (has_06cm) {
    results_df <- results_df %>%
      mutate(
        Ground_Area_06cm_m2 = round(coalesce(Pix_06cm_Grnd, 0) * pixel_area_m2, 2),
        Ground_Points_06cm = round(coalesce(PointSum_06cm_Grnd, 0)),
        Ground_Density_06cm = ifelse(Ground_Area_06cm_m2 > 0, round(Ground_Points_06cm / Ground_Area_06cm_m2, 2), 0),
        
        Canopy_Area_06cm_m2 = round(coalesce(Pix_06cm_Can, 0) * pixel_area_m2, 2),
        Canopy_Points_06cm = round(coalesce(PointSum_06cm_Can, 0)),
        Canopy_Density_06cm = ifelse(Canopy_Area_06cm_m2 > 0, round(Canopy_Points_06cm / Canopy_Area_06cm_m2, 2), 0)
      )
  } else {
    results_df <- results_df %>%
      mutate(
        Ground_Area_06cm_m2 = NA_real_, Ground_Points_06cm = NA_real_, Ground_Density_06cm = NA_real_,
        Canopy_Area_06cm_m2 = NA_real_, Canopy_Points_06cm = NA_real_, Canopy_Density_06cm = NA_real_
      )
  }
  
  # --- Strict Coverage Test ---
  if (has_3cm && has_06cm) {
    results_df <- results_df %>%
      mutate(
        Canopy_Area_06cm_Strict_m2 = round(coalesce(Pix_06cm_Can_Strict, 0) * pixel_area_m2, 2),
        Coverage_Percent = ifelse(Canopy_Area_3cm_m2 > 0, round((Canopy_Area_06cm_Strict_m2 / Canopy_Area_3cm_m2) * 100, 2), 0)
      )
  } else {
    results_df <- results_df %>%
      mutate(Canopy_Area_06cm_Strict_m2 = NA_real_, Coverage_Percent = NA_real_)
  }
  
  results_df <- results_df %>%
    select(Plot_ID, Flight_Date, Has_3cm_Flight, Has_06cm_Flight,
           Ground_Area_3cm_m2, Ground_Points_3cm, Ground_Density_3cm,
           Canopy_Area_3cm_m2, Canopy_Points_3cm, Canopy_Density_3cm,
           Ground_Area_06cm_m2, Ground_Points_06cm, Ground_Density_06cm,
           Canopy_Area_06cm_m2, Canopy_Points_06cm, Canopy_Density_06cm,
           Canopy_Area_06cm_Strict_m2, Coverage_Percent)
  
  write.csv(results_df, out_csv_path, row.names = FALSE)
  print(paste("-> Successfully saved:", basename(out_csv_path)))
  
  # Cleanup Memory
  rm(rasters_3cm, rasters_06cm, results_df)
  if (exists("can_06cm_aligned")) rm(can_06cm_aligned)
  gc()
}

print("================================================================")
print("SITE-WIDE EXTRACTION COMPLETE!")
print("================================================================")


# ──────────────────────────────────────────────────────────────────────────────
# 7. Compile All Date-Specific CSVs into a Single Master File
# ──────────────────────────────────────────────────────────────────────────────
print("================================================================")
print(" COMPILING ALL SITE-WIDE METRICS INTO MASTER CSV...")
print("================================================================")

master_output_dir <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data analysis"
master_output_path <- file.path(master_output_dir, "06. Planar Coverage and Point Density.csv")

# Create master directory if it doesn't exist
if (!dir.exists(master_output_dir)) dir.create(master_output_dir, recursive = TRUE)

# Re-scan the active dataset folders to find generated CSVs
all_analysis_csvs <- list.files(
  path = dataset_folders, 
  pattern = "^Planar_Coverage_and_Density_.*\\.csv$", 
  full.names = TRUE, 
  recursive = TRUE
)

# Ensure we restrict the paths specifically to the "12. Coverage Analysis" subfolders
valid_csvs <- all_analysis_csvs[grepl("12\\. Coverage Analysis", all_analysis_csvs)]

if (length(valid_csvs) > 0) {
  print(paste("Found", length(valid_csvs), "CSV sheets. Combining now..."))
  
  # Read and bind all data frames together
  master_df <- valid_csvs %>%
    lapply(read.csv, stringsAsFactors = FALSE) %>%
    bind_rows()
  
  # Write the unified master CSV
  write.csv(master_df, master_path <- master_output_path, row.names = FALSE)
  print(paste("🎉 SUCCESS: Master summary compiled and saved to:", master_output_path))
} else {
  print("⚠️ WARNING: No valid analysis CSV files found to combine.")
}