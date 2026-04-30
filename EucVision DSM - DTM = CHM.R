# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: BATCH PLOT-LEVEL DSM-DTM PROCESSING & DATA MERGING PIPELINE ####
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────
# Description: Automates a two-part raster math and data consolidation workflow. 
#              Part 1 generates site-wide Canopy Height Models (CHMs) by 
#              mosaicing masked Digital Surface Models (DSMs) and subtracting a 
#              master baseline DTM, subsequently extracting maximum tree heights. 
#              Part 2 aggregates these tabular metrics across all temporal 
#              flights, seamlessly joins them with ground-truth field data, 
#              applies statistical outlier filtering, and exports a cleaned, 
#              chronologically sorted Master Dataset for analysis.
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
# Install missing packages if necessary
required_packages <- c("sf", "tictoc", "dplyr", "future", "terra", "exactextractr", "stringr", "readr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(sf)
library(tictoc)
library(dplyr)
library(future)
library(terra)
library(exactextractr)
library(stringr)
library(readr)

# Force R to use English for date parsing to prevent locale-specific errors
Sys.setlocale("LC_TIME", "C")

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration & Path Management ####
# ──────────────────────────────────────────────────────────────────────────────
base_dir <- "E:/Remote Sensing Media"
baseline_dtm_path <- "E:/Remote Sensing Media/00. Baseline DTM/Ultimate_Ensemble_Baseline_DTM.tif"

baseline_dtm <- rast(baseline_dtm_path)

# --- RUN CONTROLS ---
# target_date_override <- NULL
target_date_override <- "07. 28 November 2025"

# --- EXCLUDE LIST ---
exclude_list <- c("000. Projects",
                  "00. Baseline DTM",
                  "00. Dataset Template", 
                  "01. 25 February 2025",
                  "17. 03 March 2026 (Multispectral)",
                  "20. 24 March 2026 (Multispectral)")

folders <- list.dirs(base_dir, recursive = FALSE)
dataset_folders <- folders[grepl("^\\d{2}\\.", basename(folders)) & !basename(folders) %in% exclude_list]

if (!is.null(target_date_override)) {
  dataset_folders <- dataset_folders[basename(dataset_folders) == target_date_override]
  if (length(dataset_folders) == 0) {
    stop("Target date folder not found! Please check the spelling and try again.")
  }
}

# --- MERGE DESTINATION PATHS ---
dest_master_csv <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data analysis/01. DSM - DTM = CHM/01. Master Dataset_RasterMath.csv"
field_measurements_csv <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data analysis/02. Field Measurements.csv"
dest_backup_dir <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/05. Crown Metrics/01. DSM - DTM = CHM Backups"

# ──────────────────────────────────────────────────────────────────────────────
# 3. Static Spatial Data Loading ####
# ──────────────────────────────────────────────────────────────────────────────
plots_buffered_unsorted <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LAScatalog Boundaries/LAScatalog Plot Boundaries.shp", quiet = TRUE)
plots <- plots_buffered_unsorted[order(plots_buffered_unsorted$id), ]

# NEW: Load the Impact Plot & Compartment Boundaries to mask Top and Bottom DSMs
impact_bounds <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LAScatalog Boundaries/IMPACT Plot & Compartment Boundaries 2048.shp", quiet = TRUE)
if (is.na(st_crs(impact_bounds)$epsg) || st_crs(impact_bounds)$epsg != 2048) {
  impact_bounds <- st_transform(impact_bounds, 2048)
}

# Separate the Top and Bottom vectors using the Name attribute 
top_bound <- terra::vect(impact_bounds[grepl("TOP", impact_bounds$Name, ignore.case = TRUE), ])
bot_bound <- terra::vect(impact_bounds[grepl("BOTTOM", impact_bounds$Name, ignore.case = TRUE), ])

plan(multisession, workers = 6) 

print("================================================================")
print("PART 1: STARTING RASTER MATH PROCESSING PIPELINE...")
print("================================================================")
# ──────────────────────────────────────────────────────────────────────────────
# 4. PROCESSING LOOP (DSM - DTM = CHM) ####
# ──────────────────────────────────────────────────────────────────────────────
for (folder_path in dataset_folders) {
  
  date_folder <- basename(folder_path)
  file_date <- sub("^\\d+\\.\\s*", "", date_folder)
  file_date_safe <- gsub(" ", "_", file_date)
  
  print(paste("----------------------------------------------------------------"))
  print(paste("PROCESSING DATASET:", date_folder))
  
  dsm_folder     <- file.path(folder_path, "02. Digital Surface Models")
  polygons_dir   <- file.path(folder_path, "08. Crown Polygons")
  
  # NEW: Dedicated output directory for Raster Math results
  out_dir <- file.path(folder_path, "10. DSM - DTM = CHM")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  single_chm_path <- file.path(out_dir, paste0("Master_Site_CHM_RasterMath_", file_date_safe, ".tif"))
  
  if (file.exists(single_chm_path)) {
    print(paste("-> SKIPPED: Master RasterMath CHM already exists for", date_folder))
    next
  }
  
  all_dsm_files <- list.files(dsm_folder, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
  if (length(all_dsm_files) == 0) {
    print(paste("-> SKIPPED: No .tif DSM files found in", dsm_folder))
    next
  }
  
  crown_shp_path <- file.path(polygons_dir, paste0("Crown_Polygons_", file_date_safe, ".shp"))
  if (!file.exists(crown_shp_path)) {
    print(paste("-> SKIPPED: No Crown Polygons found for", date_folder))
    next
  }
  
  trees <- st_read(crown_shp_path, quiet = TRUE)
  if (is.na(st_crs(trees)$epsg) || st_crs(trees)$epsg != 2048) {
    trees <- st_transform(trees, 2048)
  }
  trees <- trees %>% select(-any_of(c("group_ulid", "N_GM", "id", "N_FG", "N_BG", "BBox")))
  
  # ────────────────────────────────────────────────────────────────────────────
  # 4.1. DSM Selection, Masking & Mosaicing
  # ────────────────────────────────────────────────────────────────────────────
  tic("DSM Selection & Mosaicing complete")
  
  select_target_dsm <- function(files, sector_keyword) {
    subset_files <- files[grepl(sector_keyword, basename(files), ignore.case = TRUE)]
    bad_cross <- grepl("Cross", basename(subset_files), ignore.case = TRUE) & 
      !grepl("Cross Hatch", basename(subset_files), ignore.case = TRUE)
    subset_files <- subset_files[!bad_cross]
    if (length(subset_files) == 0) return(NULL)
    
    cross_hatch_files <- subset_files[grepl("Cross Hatch", basename(subset_files), ignore.case = TRUE)]
    if (length(cross_hatch_files) > 0) return(cross_hatch_files[1])
    return(subset_files[1])
  }
  
  top_dsm_path <- select_target_dsm(all_dsm_files, "Top")
  bot_dsm_path <- select_target_dsm(all_dsm_files, "Bottom")
  
  r_top <- if(!is.null(top_dsm_path)) rast(top_dsm_path) else NULL
  r_bot <- if(!is.null(bot_dsm_path)) rast(bot_dsm_path) else NULL
  
  if (is.null(r_top) && is.null(r_bot)) {
    print("-> SKIPPED: No valid DSMs matched the criteria.")
    next
  }
  
  print("Cropping and masking DSMs to compartment boundaries...")
  
  # Crop (reduces extent) and mask (turns pixels outside boundary to NA)
  if (!is.null(r_top)) r_top <- terra::crop(r_top, top_bound, mask = TRUE)
  if (!is.null(r_bot)) r_bot <- terra::crop(r_bot, bot_bound, mask = TRUE)
  
  if (!is.null(r_top) && !is.null(r_bot)) {
    print("Aligning raster grids and mosaicing Top and Bottom DSMs...")
    
    # 1. Expand the Top DSM's bounding box to include the Bottom DSM's area
    r_top_extended <- terra::extend(r_top, r_bot)
    
    # 2. Resample the Bottom DSM to fit this new expanded, perfectly aligned grid
    r_bot_aligned <- terra::resample(r_bot, r_top_extended, method = "bilinear")
    
    # 3. Mosaic them together seamlessly
    # Using fun="max" or "mean" is safe here because overlapping pixels are now NAs due to the mask
    site_dsm <- terra::mosaic(r_top_extended, r_bot_aligned, fun="mean")
    
  } else if (!is.null(r_top)) {
    site_dsm <- r_top
  } else if (!is.null(r_bot)) {
    site_dsm <- r_bot
  }
  
  toc()
  
  # ────────────────────────────────────────────────────────────────────────────
  # 4.2. Raster Math (CHM = DSM - DTM)
  # ────────────────────────────────────────────────────────────────────────────
  tic("Raster Math CHM Generation complete")
  print("Aligning Baseline DTM to weekly DSM and calculating CHM...")
  dtm_cropped <- terra::crop(baseline_dtm, site_dsm)
  dtm_aligned <- terra::resample(dtm_cropped, site_dsm, method = "bilinear")
  raw_chm <- site_dsm - dtm_aligned
  toc()
  
  # ────────────────────────────────────────────────────────────────────────────
  # 4.3. Dynamic Metric Extraction & Output to Folder 10
  # ────────────────────────────────────────────────────────────────────────────
  tic("Metric extraction complete")
  print("Extracting metrics and exporting to '10. DSM - DTM = CHM'...")
  
  trees$Tree_Height <- exact_extract(raw_chm, trees, 'max')
  
  max_tree_height <- max(trees$Tree_Height[is.finite(trees$Tree_Height)], na.rm = TRUE)
  dynamic_cap <- ceiling(max_tree_height)
  print(paste("-> Dynamic CHM cap safely set to:", dynamic_cap, "meters"))
  
  site_chm_clamped <- terra::clamp(raw_chm, lower = 0, upper = dynamic_cap)
  
  terra::writeRaster(site_chm_clamped, filename = single_chm_path, overwrite = TRUE)
  st_write(trees, file.path(out_dir, paste0("Crown_Metrics_RasterMath_", file_date_safe, ".shp")), delete_dsn = TRUE, quiet = TRUE)
  write.csv(st_drop_geometry(trees), file.path(out_dir, paste0("Crown_Metrics_RasterMath_", file_date_safe, ".csv")), row.names = FALSE)
  toc()
  
  # --- GARBAGE COLLECTION ---
  rm(dsm_list, site_dsm, dtm_cropped, dtm_aligned, raw_chm, site_chm_clamped, trees, all_dsm_files)
  gc()
}

print("================================================================")
print("PART 2: STARTING DATA EXTRACTION & MERGING PIPELINE...")
print("================================================================")

# ──────────────────────────────────────────────────────────────────────────────
# 5. Extract UAV Data from "10. DSM - DTM = CHM" ####
# ──────────────────────────────────────────────────────────────────────────────
main_folders <- list.dirs(base_dir, recursive = FALSE)
main_folders <- main_folders[!basename(main_folders) %in% exclude_list]
csv_list <- list()

for (folder in main_folders) {
  folder_name <- basename(folder)
  
  # Extract and format the date from the folder name
  date_match <- str_extract(folder_name, "\\d{2} [A-Za-z]+ \\d{4}")
  formatted_date <- NA
  if (!is.na(date_match)) {
    parsed_date <- as.Date(date_match, format="%d %B %Y")
    formatted_date <- format(parsed_date, "%d-%b-%y")
  }
  
  # Look specifically in the newly created RasterMath folder
  target_metrics_path <- file.path(folder, "10. DSM - DTM = CHM")
  
  if (dir.exists(target_metrics_path)) {
    files_to_copy <- list.files(target_metrics_path, 
                                pattern = "\\.(shp|shx|dbf|prj|csv)$", 
                                full.names = TRUE, 
                                ignore.case = TRUE)
    
    if (length(files_to_copy) > 0) {
      current_dest_dir <- file.path(dest_backup_dir, folder_name)
      if (!dir.exists(current_dest_dir)) dir.create(current_dest_dir, recursive = TRUE)
      
      # Copy files to the backup directory
      file.copy(from = files_to_copy, to = current_dest_dir, overwrite = TRUE)
      
      # Identify the CSV file for data extraction
      csv_file <- files_to_copy[grepl("\\.csv$", files_to_copy, ignore.case = TRUE)]
      
      if (length(csv_file) == 1) {
        temp_df <- read_csv(csv_file, show_col_types = FALSE)
        
        # --- Clean Columns ---
        if ("Cmprtmn" %in% names(temp_df)) temp_df <- rename(temp_df, Compartment = Cmprtmn)
        
        if ("Area_m2" %in% names(temp_df) && !"Area" %in% names(temp_df)) {
          temp_df <- rename(temp_df, Crown_Area = Area_m2)
        } else if ("Area" %in% names(temp_df) && !"Area_m2" %in% names(temp_df)) {
          temp_df <- rename(temp_df, Crown_Area = Area)
        } else if ("Area_m2" %in% names(temp_df) && "Area" %in% names(temp_df)) {
          temp_df <- mutate(temp_df, Crown_Area = coalesce(Area_m2, Area))
        }
        
        temp_df$Date <- formatted_date
        
        # Keep only essential columns
        cols_to_keep <- c("Compartment", "Line", "Plot", "Culture", "Spacing", 
                          "Species", "Tree", "Date", "Crown_Area", "Tree_Height")
        temp_df <- select(temp_df, any_of(cols_to_keep))
        
        csv_list[[folder_name]] <- temp_df
        cat(paste("SUCCESS: Extracted CSV for ->", folder_name, "\n"))
      }
    }
  }
}

# ──────────────────────────────────────────────────────────────────────────────
# 6. Merge, Full Join External Data, & Dynamic Cleaning ####
# ──────────────────────────────────────────────────────────────────────────────
if (length(csv_list) > 0) {
  cat("\nMerging RasterMath CSV files...\n")
  master_dataset <- bind_rows(csv_list)
  
  if (file.exists(field_measurements_csv)) {
    cat("Found '02. Field Measurements.csv'. Running FULL JOIN...\n")
    other_data <- read_csv(field_measurements_csv, show_col_types = FALSE)
    
    if ("Tree_Height" %in% names(other_data)) {
      other_data <- other_data %>% rename(Tree_Height_other = Tree_Height)
    }
    if ("Crown_Area" %in% names(other_data)) {
      other_data <- other_data %>% rename(Crown_Area_other = Crown_Area)
    }
    
    master_dataset <- master_dataset %>%
      full_join(other_data, by = c("Compartment", "Line", "Plot", "Culture", "Spacing", "Species", "Tree", "Date"))
    
    if ("Tree_Height_other" %in% names(master_dataset)) {
      master_dataset <- suppressWarnings(
        master_dataset %>%
          mutate(
            Tree_Height = as.numeric(Tree_Height),
            Tree_Height_other = as.numeric(Tree_Height_other),
            Tree_Height = coalesce(Tree_Height, Tree_Height_other)
          ) %>% select(-Tree_Height_other)
      )
    }
    
    if ("Crown_Area_other" %in% names(master_dataset)) {
      master_dataset <- suppressWarnings(
        master_dataset %>%
          mutate(
            Crown_Area = as.numeric(Crown_Area),
            Crown_Area_other = as.numeric(Crown_Area_other),
            Crown_Area = coalesce(Crown_Area, Crown_Area_other)
          ) %>% select(-Crown_Area_other)
      )
    }
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # 7. Outlier Filtering & Chronological Export ####
  # ────────────────────────────────────────────────────────────────────────────
  cat("Cleaning outliers and standardizing numeric formats...\n")
  
  master_dataset <- suppressWarnings(
    master_dataset %>%
      mutate(
        Tree_Height = as.numeric(Tree_Height),
        Ground_Truth_Height = as.numeric(Ground_Truth_Height),
        Crown_Area = as.numeric(Crown_Area),
        Stem_Diameter = as.numeric(Stem_Diameter)
      )
  )
  
  master_dataset <- master_dataset %>%
    group_by(Date) %>%
    mutate(
      Flight_99th = quantile(Tree_Height, probs = 0.99, na.rm = TRUE),
      Tree_Height = ifelse(Tree_Height > (Flight_99th + 5), NA, Tree_Height)
    ) %>%
    select(-Flight_99th) %>% 
    ungroup() 
  
  cat("Organizing dataset chronologically...\n")
  master_dataset <- master_dataset %>%
    mutate(Temp_Date = as.Date(Date, format="%d-%b-%y")) %>%
    arrange(Temp_Date, Compartment, Line, Plot, Tree) %>%
    select(-Temp_Date)
  
  # --- EXPORT ---
  master_csv_folder <- dirname(dest_master_csv)
  if (!dir.exists(master_csv_folder)) dir.create(master_csv_folder, recursive = TRUE)
  
  tryCatch({
    write_csv(master_dataset, dest_master_csv)
    cat(paste("\nDONE: Master RasterMath dataset merged, cleaned, and exported with", nrow(master_dataset), "rows.\n"))
    cat(paste("Location:", dest_master_csv, "\n"))
  }, error = function(e) {
    cat("\nERROR: Could not overwrite the master CSV. Is it open in Excel?\n")
  })
  
} else {
  cat("\nNo CSV files were found in the 10. DSM - DTM = CHM folders to merge.\n")
}