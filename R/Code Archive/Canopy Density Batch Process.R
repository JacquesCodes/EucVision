# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: DUAL-GSD ON-THE-FLY DENSITY EXTRACTION (SPACE SAVER)
# ──────────────────────────────────────────────────────────────────────────────
# Description: Automates the extraction of Canopy (3D) and Ground (2D) point 
#              densities. Now explicitly checks for and separately processes 
#              nested "0.6 GSD" subfolders to ensure same-day multi-GSD flights 
#              are kept statistically isolated.
# ──────────────────────────────────────────────────────────────────────────────

library(lidR)
library(sf)
library(terra)
library(dplyr)
library(stringr)
library(future)

# ──────────────────────────────────────────────────────────────────────────────
# 1. Configuration & Static Data
# ──────────────────────────────────────────────────────────────────────────────
base_dir <- "E:/Remote Sensing Media"
baseline_dtm_path <- "E:/Remote Sensing Media/00. Baseline DTM/Ultimate_Ensemble_Baseline_DTM_ESRI_102562.tif"
plots_shp_path <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. QGIS Shapefiles/1. LAScatalog Plot Boundaries/LAScatalog_Plot_Boundaries_ESRI_102562.shp"

voxel_resolution <- 0.1  # 10cm voxels
z_threshold <- 0.3       # Locked-in empty trunk space threshold

# Set to a specific folder to test, or NULL for full batch
target_date_override <- NULL

print("Loading static spatial data...")
baseline_dtm <- rast(baseline_dtm_path)
plots_buffered_unsorted <- st_read(plots_shp_path, quiet = TRUE)
st_crs(plots_buffered_unsorted) <- st_crs(baseline_dtm) 
plots <- plots_buffered_unsorted[order(plots_buffered_unsorted$id), ]

# Enable parallel processing for lidR 
plan(multisession, workers = 6)

exclude_list <- c("000. Projects", "00. Baseline DTM", "00. Dataset Template", 
                  "01. 25 February 2025", "07. December 2025 (TLS)",
                  "17. 02 March 2026","17. 02 March 2026 2.4","17. 02 March 2026 19.2","20. 23 March 2026 0.6cm","30. 30 June 2026 (ALS)")

folders <- list.dirs(base_dir, recursive = FALSE)
dataset_folders <- folders[grepl("^\\d{2}\\.", basename(folders)) & !basename(folders) %in% exclude_list]

if (!is.null(target_date_override)) {
  dataset_folders <- dataset_folders[basename(dataset_folders) == target_date_override]
}

# ──────────────────────────────────────────────────────────────────────────────
# 2. Master Batch Loop
# ──────────────────────────────────────────────────────────────────────────────
for (folder_path in dataset_folders) {
  
  date_folder <- basename(folder_path)
  file_date <- sub("^\\d+\\.\\s*", "", date_folder)
  file_date_safe <- gsub(" ", "_", file_date)
  
  print(paste("================================================================"))
  print(paste("PROCESSING:", date_folder))
  print(paste("================================================================"))
  
  # Define base directories
  las_folder_main <- file.path(folder_path, "03. Point Clouds")
  las_folder_fine <- file.path(las_folder_main, "0.6 GSD")
  out_dir         <- file.path(folder_path, "10. Density Metrics")
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Map out the flight variations present for this date
  flights_to_process <- list()
  
  # 1. Grab main folder files (recursive = FALSE so it ignores the 0.6 folder)
  main_files <- list.files(las_folder_main, pattern = "\\.(las|laz)$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
  if (length(main_files) > 0) {
    flights_to_process[["Standard"]] <- main_files
  }
  
  # 2. Grab 0.6cm subfolder files if they exist
  if (dir.exists(las_folder_fine)) {
    fine_files <- list.files(las_folder_fine, pattern = "\\.(las|laz)$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
    if (length(fine_files) > 0) {
      flights_to_process[["Fine_0.6cm"]] <- fine_files
    }
  }
  
  if (length(flights_to_process) == 0) {
    print("-> SKIPPED: No raw point clouds found.")
    next
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # 3. Sub-Flight Loop (Processes Main and Fine GSDs separately)
  # ────────────────────────────────────────────────────────────────────────────
  for (flight_type in names(flights_to_process)) {
    
    print(paste("-> Initiating sub-routine for:", flight_type, "flight..."))
    current_las_files <- flights_to_process[[flight_type]]
    
    # Tag the CSV so they don't overwrite each other
    out_csv_path <- file.path(out_dir, paste0("Plot_Densities_", file_date_safe, "_", flight_type, ".csv"))
    
    if (file.exists(out_csv_path)) {
      print(paste("   -> SKIPPED:", flight_type, "metrics already exist."))
      next
    }
    
    # Temporary directories scoped to the flight type
    temp_clip_dir  <- file.path(folder_path, paste0("TEMP_Clipped_", flight_type))
    temp_norm_dir  <- file.path(folder_path, paste0("TEMP_Normalised_", flight_type))
    dir.create(temp_clip_dir, showWarnings = FALSE)
    dir.create(temp_norm_dir, showWarnings = FALSE)
    
    # --- STEP A: Cropping ---
    top_files <- current_las_files[grepl("Top", basename(current_las_files), ignore.case = TRUE)]
    bot_files <- current_las_files[grepl("Bottom", basename(current_las_files), ignore.case = TRUE)]
    
    if (length(top_files) > 0 && length(bot_files) > 0) {
      plots_top <- plots %>% filter(id <= 21)
      plots_bot <- plots %>% filter(id >= 22)
      
      ctg_top <- readLAScatalog(top_files)
      opt_independent_files(ctg_top) <- FALSE
      opt_select(ctg_top) <- "xyz"
      opt_output_files(ctg_top) <- file.path(temp_clip_dir, paste0("Plot_{id}_", file_date_safe))
      suppressMessages(clip_roi(ctg_top, plots_top))
      
      ctg_bot <- readLAScatalog(bot_files)
      opt_independent_files(ctg_bot) <- FALSE
      opt_select(ctg_bot) <- "xyz"
      opt_output_files(ctg_bot) <- file.path(temp_clip_dir, paste0("Plot_{id}_", file_date_safe))
      suppressMessages(clip_roi(ctg_bot, plots_bot))
    } else {
      ctg <- readLAScatalog(current_las_files)
      opt_independent_files(ctg) <- FALSE
      opt_select(ctg) <- "xyz"
      opt_output_files(ctg) <- file.path(temp_clip_dir, paste0("Plot_{id}_", file_date_safe))
      suppressMessages(clip_roi(ctg, plots))
    }
    
    # --- STEP B: Normalization ---
    ctg_clipped <- readLAScatalog(temp_clip_dir)
    opt_independent_files(ctg_clipped) <- TRUE
    opt_select(ctg_clipped) <- "xyz"
    opt_output_files(ctg_clipped) <- file.path(temp_norm_dir, "{ORIGINALFILENAME}_norm")
    
    suppressWarnings(normalize_height(las = ctg_clipped, algorithm = baseline_dtm))
    
    # --- STEP C: Density Extraction ---
    norm_files <- list.files(temp_norm_dir, pattern = "\\.(las|laz)$", full.names = TRUE)
    plot_results <- list()
    
    for (file in norm_files) {
      plot_id_match <- str_extract(basename(file), "Plot_\\d+")
      plot_id <- ifelse(is.na(plot_id_match), basename(file), plot_id_match)
      
      las <- readLAS(file)
      if (is.empty(las)) next
      
      # Ground Density (2D)
      ground_las <- filter_poi(las, Z <= z_threshold)
      if (!is.empty(ground_las)) {
        ground_density_m2 <- npoints(ground_las) / area(ground_las)
      } else {
        ground_density_m2 <- NA
      }
      
      # Canopy Density (3D Voxels)
      canopy_las <- filter_poi(las, Z > z_threshold)
      if (!is.empty(canopy_las)) {
        voxels <- voxel_metrics(canopy_las, ~length(Z), res = voxel_resolution)
        single_voxel_vol <- voxel_resolution^3
        canopy_density_m3 <- sum(voxels$V1) / (nrow(voxels) * single_voxel_vol)
      } else {
        canopy_density_m3 <- NA
      }
      
      plot_results[[length(plot_results) + 1]] <- data.frame(
        Flight_Date = file_date,
        Flight_Type = flight_type,
        Plot_ID = plot_id,
        Ground_Density_pts_m2 = round(ground_density_m2, 2),
        Canopy_Density_pts_m3 = round(canopy_density_m3, 2)
      )
    }
    
    if (length(plot_results) > 0) {
      flight_df <- bind_rows(plot_results)
      write.csv(flight_df, out_csv_path, row.names = FALSE)
      print(paste("   -> Successfully saved:", basename(out_csv_path)))
    }
    
    # --- STEP D: Housekeeping ---
    unlink(temp_clip_dir, recursive = TRUE)
    unlink(temp_norm_dir, recursive = TRUE)
    rm(ctg_clipped, norm_files, plot_results, flight_df)
    gc()
  }
}

print("================================================================")
print("ALL DONE! Dual-GSD logic processed and disks remain uncluttered.")
print("================================================================")