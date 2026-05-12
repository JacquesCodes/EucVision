# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: BATCH SPATIAL DATA MERGING & CSV INTEGRATION PIPELINE
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────
# Description: Automates the batch processing of multiple dates. Aggregates 
#              individual spatial plot shapefiles into a single continuous 
#              spatial dataframe, applies temporal mortality filtering, and 
#              exports the cleaned vectors. Includes shapefile routing overrides.
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
library(sf)      # Core package for spatial vector operations
library(dplyr)   # Core package for data wrangling and piping
library(tictoc)  # For tracking script execution time

tic()
print("Initiating Batch Spatial Data Merging & CSV Integration...")

# Force English locale for consistent date parsing
Sys.setlocale("LC_TIME", "C")

# The strict, exact OGC Well-Known Text (WKT) for EPSG:2048 WKT
pure_epsg_2048_wkt <- 'PROJCS["Hartebeesthoek94 / Lo19",GEOGCS["Hartebeesthoek94",DATUM["Hartebeesthoek94",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6148"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4148"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",19],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Southing",SOUTH],AXIS["Westing",WEST],AUTHORITY["EPSG","2048"]]'

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration & Batch Management ####
# ──────────────────────────────────────────────────────────────────────────────
# Define the base directories
base_input <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/03. QGIS Extracted data"
base_output <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/04. QGIS Combined Output"
csv_path <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data analysis/00. Dataset template.csv"

# --- RUN CONTROLS ---
# Set to a specific folder name to run only that dataset, or set to NULL for full batch.
target_date_override <- "24. 23 April 2026"

# --- EXCLUDE LIST ---
exclude_list <- c("01. 25 February 2025",
                  "07. December 2025 (TLS)")

# Scan the base directory and filter for valid date folders
folders <- list.dirs(base_input, recursive = FALSE)
dataset_folders <- folders[grepl("^\\d{2}\\.", basename(folders)) & !basename(folders) %in% exclude_list]

if (!is.null(target_date_override)) {
  dataset_folders <- dataset_folders[basename(dataset_folders) == target_date_override]
}

# ──────────────────────────────────────────────────────────────────────────────
# 3. CSV Spatial Baseline Locking ####
# ──────────────────────────────────────────────────────────────────────────────
# Load the master CSV baseline once into memory before the loop
master_csv_data <- read.csv(csv_path) %>% mutate(Tree = round(as.numeric(Tree), 2))
master_csv_data$Parsed_Death_Date <- as.Date(master_csv_data$Death_Date, format="%d-%m-%Y")

# --- CRITICAL FIX: Lock the baseline to the 3144 shapefile features ---
shapefile_baseline_date <- as.Date("2025-09-01")
csv_shapefile_baseline <- master_csv_data %>%
  filter(is.na(Parsed_Death_Date) | Parsed_Death_Date > shapefile_baseline_date)

legacy_sf_count <- nrow(csv_shapefile_baseline) # This guarantees exactly 3144!
print(paste("Global Spatial Baseline Set To:", legacy_sf_count, "trees."))


# ──────────────────────────────────────────────────────────────────────────────
# MASTER BATCH LOOP START ####
# ──────────────────────────────────────────────────────────────────────────────
for (folder_path in dataset_folders) {
  
  date_folder <- basename(folder_path)
  file_date <- sub("^\\d+\\.\\s*", "", date_folder)
  file_date_safe <- gsub(" ", "_", file_date)
  
  print("================================================================")
  print(paste("PROCESSING DATASET:", date_folder))
  print("================================================================")
  
  current_flight_date <- as.Date(file_date, format="%d %B %Y")
  if(is.na(current_flight_date)) {
    print(paste("-> WARNING: Date parsing failed. Skipping", date_folder))
    next
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # 4. Shapefile Override Routing ####
  # ────────────────────────────────────────────────────────────────────────────
  shapefile_source_folder <- date_folder
  
  if (date_folder == "05. 14 November 2025") {
    shapefile_source_folder <- "06. 17 November 2025"
    print(paste("-> ROUTING OVERRIDE: Using donor shapefiles from", shapefile_source_folder))
  } else if (date_folder == "19. 16 March 2026") {
    shapefile_source_folder <- "18. 09 March 2026"
    print(paste("-> ROUTING OVERRIDE: Using donor shapefiles from", shapefile_source_folder))
  } else if (date_folder == "22. 08 April 2026") {
    shapefile_source_folder <- "21. 31 March 2026"
    print(paste("-> ROUTING OVERRIDE: Using donor shapefiles from", shapefile_source_folder))
  } else if (date_folder == "23. 13 April 2026") {
    shapefile_source_folder <- "21. 31 March 2026"
    print(paste("-> ROUTING OVERRIDE: Using donor shapefiles from", shapefile_source_folder))
  }
  
  # Dynamically construct paths
  input_folder <- file.path(base_input, shapefile_source_folder)
  output_folder <- file.path(base_output, date_folder)
  
  output_shp <- file.path(output_folder, paste0("Crown_Polygons_", file_date_safe, ".shp")) 
  base_output_file <- paste0("E:/Remote Sensing Media/", date_folder, "/08. Crown Polygons/Crown_Polygons_", file_date_safe, ".shp") 
  
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  if (!dir.exists(dirname(base_output_file))) dir.create(dirname(base_output_file), recursive = TRUE)
  
  # ────────────────────────────────────────────────────────────────────────────
  # 5. Shapefile Loading & Merging ####
  # ────────────────────────────────────────────────────────────────────────────
  shp_files <- list.files(input_folder, pattern = "\\.shp$", full.names = TRUE)
  
  if (length(shp_files) == 0) {
    print(paste("-> WARNING: No shapefiles found in", input_folder, "- Skipping dataset."))
    next
  }
  
  # Sort files numerically
  plot_numbers <- as.numeric(gsub("\\D", "", basename(shp_files)))
  shp_files <- shp_files[order(plot_numbers)]
  
  combined_sf <- lapply(shp_files, function(file) {
    temp_sf <- st_read(file, quiet = TRUE)
    plot_name <- tools::file_path_sans_ext(basename(file))
    temp_sf <- temp_sf %>% mutate(Plot_shp = plot_name) 
    return(temp_sf)
  }) %>% bind_rows() 
  
  # ────────────────────────────────────────────────────────────────────────────
  # 6. Temporal Mortality Filtering & Failsafes ####
  # ────────────────────────────────────────────────────────────────────────────
  # Refresh the CSV baseline for this iteration
  csv_data <- csv_shapefile_baseline
  
  # A tree is alive if it has no death date OR if its death date is AFTER the current flight
  csv_data$Is_Alive <- is.na(csv_data$Parsed_Death_Date) | (csv_data$Parsed_Death_Date > current_flight_date)
  
  csv_alive <- csv_data %>% filter(Is_Alive)
  alive_count <- nrow(csv_alive)
  
  print(paste("   Legacy Shapefile Baseline:", legacy_sf_count))
  print(paste("   Alive on", file_date, ":", alive_count))
  print(paste("   Trees flagged as Dead since baseline:", legacy_sf_count - alive_count))
  
  # ---------------------------------------------------------
  # FAILSAFE 1: Global Count Check & Diagnostic Output
  # ---------------------------------------------------------
  raw_sf_count <- nrow(combined_sf)
  is_legacy <- FALSE
  
  if (raw_sf_count == legacy_sf_count) {
    print("   -> FAILSAFE 1 PASSED: Found legacy shapefiles (3144 features).")
    is_legacy <- TRUE
    csv_for_binding <- csv_data
    
  } else if (raw_sf_count == alive_count) {
    print("   -> FAILSAFE 1 PASSED: Found updated shapefiles. Features match alive trees.")
    is_legacy <- FALSE
    csv_for_binding <- csv_alive
    
  } else {
    # --- NEW DIAGNOSTIC CODE ---
    print("\n--- DIAGNOSTIC: FINDING MISMATCHED PLOTS FOR FAILSAFE 1 ---")
    
    # Calculate expected counts based on alive trees for this date
    diag_csv_counts <- csv_alive %>% 
      group_by(Plot) %>% 
      summarise(Expected = n(), .groups = "drop")
    
    # Calculate actual counts from the shapefiles currently loaded
    diag_sf_counts <- combined_sf %>% 
      st_drop_geometry() %>% 
      mutate(Plot = as.numeric(gsub("\\D", "", Plot_shp))) %>% 
      group_by(Plot) %>% 
      summarise(Actual = n(), .groups = "drop")
    
    # Join and find where the numbers don't match
    diag_mismatch <- full_join(diag_csv_counts, diag_sf_counts, by="Plot") %>% 
      filter(is.na(Expected) | is.na(Actual) | Expected != Actual)
    
    print("Here are the plots with missing or extra features:")
    print(diag_mismatch, n = Inf) # n = Inf ensures all mismatched rows print to console
    print("-----------------------------------------------------------")
    # ---------------------------
    
    stop(paste("\nCRITICAL ERROR - FAILSAFE 1 TRIGGERED FOR", date_folder, 
               "\nShapefile features (", raw_sf_count, ") do NOT match the alive trees (", alive_count, ").",
               "\nBecause it does not match alive trees, it MUST be exactly the legacy baseline (", legacy_sf_count, ").",
               "\n-> Check the console output directly above this error to see which plots are mismatched!"))
  }
  
  # ---------------------------------------------------------
  # FAILSAFE 2: Pre-Bind Plot Consistency Check
  # ---------------------------------------------------------
  # Before gluing the datasets together, verify the per-plot counts match perfectly!
  csv_plot_counts <- csv_for_binding %>% group_by(Plot) %>% summarise(Expected = n(), .groups = "drop")
  
  # Extract the Plot Number dynamically from the 'Plot_shp' attribute
  sf_plot_counts <- combined_sf %>% 
    st_drop_geometry() %>% 
    mutate(Plot = as.numeric(gsub("\\D", "", Plot_shp))) %>% 
    group_by(Plot) %>% 
    summarise(Actual = n(), .groups = "drop")
  
  plot_mismatch <- full_join(csv_plot_counts, sf_plot_counts, by="Plot") %>% 
    filter(is.na(Expected) | is.na(Actual) | Expected != Actual)
  
  if (nrow(plot_mismatch) > 0) {
    print("\n--- MISMATCHED PLOTS DETECTED ---")
    print(plot_mismatch)
    stop(paste("CRITICAL ERROR - FAILSAFE 2 TRIGGERED in", date_folder, 
               "- The shapefile features per plot do not perfectly match the CSV rows per plot! Halting before bad merge."))
  } else {
    print("   -> FAILSAFE 2 PASSED: CSV rows perfectly align with spatial plots.")
  }
  
  # ---------------------------------------------------------
  # Data Merging (Safe to proceed)
  # ---------------------------------------------------------
  if (is_legacy) {
    print("   -> Binding data and filtering out dead geometries...")
    combined_sf <- bind_cols(csv_for_binding, combined_sf) %>% st_as_sf()
    combined_sf <- combined_sf %>% filter(Is_Alive)
  } else {
    print("   -> Binding updated data...")
    combined_sf <- bind_cols(csv_for_binding, combined_sf) %>% st_as_sf()
  }
  
  # ---------------------------------------------------------
  # FAILSAFE 3: Final Spatial Output Check
  # ---------------------------------------------------------
  if (nrow(combined_sf) != alive_count) {
    stop(paste("CRITICAL ERROR - FAILSAFE 3 TRIGGERED in", date_folder, "- Final count mismatch after processing!"))
  } else {
    print("   -> FAILSAFE 3 PASSED: Final geometry count matches alive trees.")
  }
  
  # Clean up temporary structural columns
  combined_sf <- combined_sf %>% select(-Parsed_Death_Date, -Is_Alive, -Plot_shp)
  
  # ────────────────────────────────────────────────────────────────────────────
  # 7. Export Processed Data & Enforce CRS ####
  # ────────────────────────────────────────────────────────────────────────────
  st_write(combined_sf, output_shp, append = FALSE, quiet = TRUE)
  st_write(combined_sf, base_output_file, append = FALSE, quiet = TRUE)
  
  prj_teams <- sub("\\.shp$", ".prj", output_shp, ignore.case = TRUE)
  prj_ssd <- sub("\\.shp$", ".prj", base_output_file, ignore.case = TRUE)
  
  writeLines(pure_epsg_2048_wkt, prj_teams)
  writeLines(pure_epsg_2048_wkt, prj_ssd)
  
  print("-> All Failsafes passed. Geometries filtered and exported.")
  
  # Clean up memory before the next loop
  rm(combined_sf, csv_for_binding, csv_plot_counts, sf_plot_counts, plot_mismatch)
  gc()
}

print("================================================================")
print("BATCH PIPELINE COMPLETE! All datasets safely merged and filtered.")
print("================================================================")
toc()