# ──────────────────────────────────────────────────────────────────────────────
# BATCH AUTOMATED TREE DETECTION (ITD) PIPELINE & PiP VALIDATION
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Project: EucXylo (https://eucxylo.sun.ac.za/) 
# ──────────────────────────────────────────────────────────────────────────────
# Description: Automates the Local Maximum Filter (LMF) tree detection algorithm
#              across all temporal CHMs using dynamically scaled window sizes. 
#              Validates accuracy using a Point-in-Polygon (PiP) spatial join 
#              against manual Geo-SAM reference polygons to explicitly classify 
#              True Positives, Over-segmentation errors, and Missed trees.
#              NEW: Calculates Intersection over Union (IoU) agreement between 
#              Geo-SAM polygons and algorithmic (Watershed) CHM polygons.
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
library(sf)            # Spatial vector data handling
library(terra)         # Spatial raster data handling
library(dplyr)         # Data wrangling and piping logic
library(lidR)          # Individual Tree Detection algorithms
library(ForestTools)   # Marker-Controlled Watershed (mcws) for crown polygons
library(tictoc)        # Script execution timing 
library(stringr)       # String manipulation for dates
library(tidyr)

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration & Batch Management ####
# ──────────────────────────────────────────────────────────────────────────────
base_dir <- "E:/Remote Sensing Media"

# Final output CSV for the complete longitudinal analysis
output_csv <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data Analysis/10. Tree Detection_with_IoU.csv"

# --- RUN CONTROLS ---
# Set this to your TLS or ALS folder name to run the validation on just that date
target_date_override <- "20. 23 March 2026"

# Folders to ignore during the batch processing loop
exclude_list <- c("000. Projects",
                  "00. Baseline DTM",
                  "00. Dataset Template", 
                  #"07. December 2025 (TLS)",
                  "17. 02 March 2026 2.4",
                  "17. 02 March 2026 4.8",
                  "17. 02 March 2026 19.2",
                  "17. 02 March 2026 Double Grid",
                  "20. 23 March 2026 0.6cm",
                  "31. 26 June 2026 Oblique",
                  # "31. 30 June 2026 (ALS)",
                  "40. 12 August 2026",
                  "40. 12 August 2026 Terra")

# Scan the base directory and filter for valid date folders
folders <- list.dirs(base_dir, recursive = FALSE)
dataset_folders <- folders[grepl("^\\d{2}\\.", basename(folders)) & !basename(folders) %in% exclude_list]

if (!is.null(target_date_override)) {
  dataset_folders <- dataset_folders[basename(dataset_folders) == target_date_override]
  if (length(dataset_folders) == 0) {
    stop("CRITICAL ERROR: Target date folder not found! Please check spelling.")
  }
}

# Master list to hold all results before writing to CSV
all_results <- list()

print("Starting Batch ITD PiP & IoU Validation Pipeline...")

# ──────────────────────────────────────────────────────────────────────────────
# MASTER BATCH LOOP START ####
# ──────────────────────────────────────────────────────────────────────────────
for (folder_path in dataset_folders) {
  
  date_folder <- basename(folder_path)
  file_date_safe <- gsub(" ", "_", sub("^\\d+\\.\\s*", "", date_folder))
  
  # Extract true date format for the CSV
  date_match <- str_extract(date_folder, "\\d{2} [A-Za-z]+ \\d{4}")
  formatted_date <- format(as.Date(date_match, format="%d %B %Y"), "%d-%m-%Y")
  
  print("================================================================")
  print(paste("RUNNING TREE DETECTION FOR:", date_folder))
  print("================================================================")
  
  polygons_dir <- file.path(folder_path, "08. Crown Polygons")
  chm_dir      <- file.path(folder_path, "07. Canopy Height Models")
  
  path_trees <- file.path(polygons_dir, paste0("Crown_Polygons_", file_date_safe, ".shp"))
  path_chm   <- file.path(chm_dir, paste0("Master_Site_CHM_Single_", file_date_safe, ".tif"))
  
  if (date_folder == "31. 30 June 2026 (ALS)"){
    
    path_chm   <- file.path(chm_dir, paste0("A_Smoothed_Master_Site_CHM_Single_", file_date_safe, ".tif"))
    
  }
  
  if (date_folder == "20. 23 March 2026"){
    
    path_chm   <- file.path(chm_dir, paste0("A_Smoothed_Master_Site_CHM_Single_", file_date_safe, ".tif"))
    
  }
  
  if (!file.exists(path_trees)) {
    print(paste("-> SKIPPED: No Crown Polygons found for", date_folder))
    next
  }
  if (!file.exists(path_chm)) {
    print(paste("-> SKIPPED: No Master CHM found for", date_folder))
    next
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # 3. Spatial Data Loading & Per-Plot Detection ####
  # ────────────────────────────────────────────────────────────────────────────
  tic("ITD PiP & IoU processing complete")
  
  trees <- st_read(path_trees, quiet = TRUE)
  
  # Standardize the Compartment column name if it exported as Cmprtmn
  if ("Cmprtmn" %in% names(trees)) {
    trees <- trees %>% rename(Compartment = Cmprtmn)
  }
  
  ctg_chm <- rast(path_chm)
  st_crs(trees) <- st_crs(ctg_chm)
  
  unique_plots <- unique(trees$Plot)
  
  # Loop through every plot individually for this specific date
  for (p in unique_plots) {
    
    # 1. Isolate the current plot's trees
    plot_trees <- trees %>% filter(Plot == p)
    actual_count <- nrow(plot_trees)
    
    if(actual_count == 0) next
    
    # Extract structural metadata safely from the first row of the plot
    comp_val    <- if("Compartment" %in% names(plot_trees)) plot_trees$Compartment[1] else NA
    line_val    <- if("Line" %in% names(plot_trees)) plot_trees$Line[1] else NA
    culture_val <- if("Culture" %in% names(plot_trees)) plot_trees$Culture[1] else NA
    spacing_val <- if("Spacing" %in% names(plot_trees)) plot_trees$Spacing[1] else NA
    species_val <- if("Species" %in% names(plot_trees)) plot_trees$Species[1] else NA
    
    # Override Species to "Mix" if Culture is "Mix"
    if (!is.na(culture_val) && culture_val == "Mix") {
      species_val <- "Mix"
    }
    
    # 2. Crop and Mask the Master CHM to the true angled plot boundary
    
    # Crop to the rectangular extent first (to reduce memory load)
    plot_extent <- ext(plot_trees) + 3
    chm_cropped <- crop(ctg_chm, plot_extent)
    
    # Create a tight, angled polygon (Convex Hull) wrapping this specific plot's trees
    plot_boundary_angled <- st_convex_hull(st_union(plot_trees))
    
    # Buffer the angled boundary by 2m so we don't slice off the outer edges of the border crowns
    plot_boundary_buffered <- st_buffer(plot_boundary_angled, 2)
    
    # Mask the CHM: Everything outside your diagonal plot boundary becomes NA
    chm_cropped <- mask(chm_cropped, vect(plot_boundary_buffered))
    
    # 3. Dynamic Algorithm Configuration
    ws_val <- ifelse(is.na(spacing_val) || spacing_val <= 0, 3, spacing_val * 0.8)
    lmf_dynamic <- lmf(ws = ws_val, hmin = 0.5, shape = "circular")
    
    ttops <- tryCatch({
      locate_trees(las = chm_cropped, algorithm = lmf_dynamic)
    }, error = function(e) { NULL })
    
    # 4. Point-in-Polygon (PiP) Validation
    total_detected_points <- 0
    true_positives <- 0
    over_segmented <- 0
    missed_trees <- actual_count 
    
    if (!is.null(ttops) && nrow(ttops) > 0) {
      plot_boundary <- st_union(st_make_valid(plot_trees))
      ttops_in_plot <- st_filter(ttops, plot_boundary)
      total_detected_points <- nrow(ttops_in_plot)
      
      points_per_poly <- lengths(st_intersects(plot_trees, ttops_in_plot))
      
      true_positives <- sum(points_per_poly == 1)
      over_segmented <- sum(points_per_poly > 1)
      missed_trees   <- sum(points_per_poly == 0)
    }
    
    # 5. NEW: Calculate Intersection over Union (IoU) for Trees > 1m
    plot_mean_iou <- NA
    if (!is.null(ttops) && nrow(ttops) > 0) {
      
      # Mask out anything below 0.5m to stop ground-spill of the watershed algorithm
      chm_masked <- chm_cropped
      chm_masked[chm_masked < 0.5] <- NA
      
      # Draw algorithmic polygons
      tls_crowns_spat <- tryCatch({
        mcws(treetops = ttops, CHM = chm_masked, minHeight = 0.5, format = "polygons")
      }, error = function(e) { NULL })
      
      if (!is.null(tls_crowns_spat)) {
        tls_crowns <- st_as_sf(tls_crowns_spat)
        tls_crowns <- st_make_valid(tls_crowns)
        plot_trees_valid <- st_make_valid(plot_trees)
        
        # --- NEW HEIGHT FILTERING LOGIC ---
        # Extract the maximum CHM height for every Geo-SAM polygon
        height_extract <- terra::extract(chm_cropped, vect(plot_trees_valid), fun = max, na.rm = TRUE)
        plot_trees_valid$max_h <- height_extract[, 2]
        
        # Filter Geo-SAM polygons to only keep trees >= 1.0 meter
        plot_trees_filtered <- plot_trees_valid %>% filter(max_h >= 1.0)
        # ----------------------------------
        
        # --- EXCLUDE EDGE-CLIPPED TREES ---
        # 1. Get the actual non-NA footprint of the cropped CHM
        chm_footprint <- as.polygons(chm_cropped > -Inf, dissolve = TRUE)
        chm_footprint_sf <- st_as_sf(chm_footprint)
        
        # 2. Inward buffer by 0.5m so trees touching the edge are excluded
        valid_interior <- st_buffer(chm_footprint_sf, -0.5)
        
        # 3. Keep only Geo-SAM polygons completely inside the valid scan area
        trees_within_extent <- st_contains(valid_interior, plot_trees_filtered, sparse = FALSE)
        plot_trees_filtered <- plot_trees_filtered[trees_within_extent[1, ], ]
        # -----------------------------------
        
        # Only proceed if there are still trees left after filtering
        if (nrow(plot_trees_filtered) > 0) {
          
          # Give filtered Geo-SAM polygons a unique ID and pre-calculate areas
          plot_trees_filtered$uid <- 1:nrow(plot_trees_filtered)
          plot_trees_filtered$area_ref <- st_area(plot_trees_filtered)
          tls_crowns$area_alg <- st_area(tls_crowns)
          
          # Intersect the two polygon layers using the FILTERED trees
          intersections <- suppressWarnings(st_intersection(plot_trees_filtered, tls_crowns))
          
          if (nrow(intersections) > 0) {
            intersections$int_area <- st_area(intersections)
            
            # For every Geo-SAM polygon, find the TLS polygon that overlaps it the most
            best_matches <- intersections %>%
              group_by(uid) %>%
              slice_max(order_by = int_area, n = 1, with_ties = FALSE) %>%
              ungroup()
            
            # Compute IoU: Area of Intersection / Area of Union
            best_matches$union_area <- best_matches$area_ref + best_matches$area_alg - best_matches$int_area
            best_matches$iou <- as.numeric(best_matches$int_area / best_matches$union_area)
            
            # Calculate the mean IoU for this plot (now strictly for trees > 1m)
            plot_mean_iou <- round(mean(best_matches$iou, na.rm = TRUE), 4)
          }
        }
        
        # --- EXPORT ALGORITHMIC POLYGONS FOR QGIS ---
        val_out_dir <- file.path(folder_path, "08.1 Validation Algorithmic Crowns")
        if (!dir.exists(val_out_dir)) dir.create(val_out_dir, recursive = TRUE)
        
        out_shp <- file.path(val_out_dir, paste0("TLS_Watershed_Plot_", p, ".shp"))
        suppressWarnings(st_write(tls_crowns, out_shp, append = FALSE, quiet = TRUE))
        
        # --- QUICK VISUALIZATION IN R ---
        plot(chm_cropped, main = paste("Plot", p, "- IoU (>1m):", plot_mean_iou))
        # Plot filtered Geo-SAM polygons in BLUE (Trees > 1m)
        plot(st_geometry(plot_trees_filtered), border = "blue", lwd = 2, add = TRUE)
        # Plot algorithmic TLS polygons in RED
        plot(st_geometry(tls_crowns), border = "red", lwd = 2, add = TRUE)
      }
    
    }
    
    # 6. Save the plot's metrics
    all_results[[length(all_results) + 1]] <- data.frame(
      Date = formatted_date,
      Compartment = comp_val,
      Line = line_val,
      Plot = p,
      Culture = culture_val,
      Spacing = spacing_val,
      Species = species_val,
      Actual_Trees = actual_count,
      Total_Detected_Points = total_detected_points,
      True_Positives = true_positives,
      Over_Segmented_Trees = over_segmented,
      Missed_Trees = missed_trees,
      Raw_Detection_Rate_Pct = round((total_detected_points / actual_count) * 100, 2),
      True_Detection_Rate_Pct = round((true_positives / actual_count) * 100, 2),
      Mean_IoU = plot_mean_iou # <--- New column added here
    )
  }
  
  toc()
  
  # Clean memory for the next iteration
  rm(trees, ctg_chm)
  gc()
}

# ──────────────────────────────────────────────────────────────────────────────
# 4. Export Final CSV ####
# ──────────────────────────────────────────────────────────────────────────────
if (length(all_results) > 0) {
  final_df <- bind_rows(all_results)
  
  # Ensure the output directory exists
  out_dir <- dirname(output_csv)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  write.csv(final_df, output_csv, row.names = FALSE)
  print("================================================================")
  print(paste("BATCH EXTRACTION COMPLETE! Saved to:", output_csv))
  print("================================================================")
} else {
  print("No data was extracted. Check folder paths.")
}

# ──────────────────────────────────────────────────────────────────────────────
# 5. Generate Summary Table ####
# ──────────────────────────────────────────────────────────────────────────────

# 1. Load your newly generated PiP results
df_pip <- read.csv("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data Analysis/10. Tree Detection_with_IoU.csv")

# 2. Convert Date string to Date object to easily filter the first and last flights
df_pip$Date <- as.Date(df_pip$Date, format="%d-%m-%Y")

# Find your exact baseline and final dates
baseline_date <- min(df_pip$Date, na.rm = TRUE)
final_date <- max(df_pip$Date, na.rm = TRUE)

# 3. Filter to just these two extremes and calculate percentages
table_data <- df_pip %>%
  filter(Date %in% c(baseline_date, final_date)) %>%
  mutate(
    Phase = ifelse(Date == baseline_date, "Baseline", "Final"),
    Omission_Pct = (Missed_Trees / Actual_Trees) * 100,
    Commission_Pct = (Over_Segmented_Trees / Actual_Trees) * 100
  ) %>%
  # 4. Group and summarize by Spacing and Phase
  group_by(Spacing, Phase) %>%
  summarise(
    Mean_True_Detection = round(mean(True_Detection_Rate_Pct, na.rm = TRUE), 1),
    Mean_Omission = round(mean(Omission_Pct, na.rm = TRUE), 1),
    Mean_Commission = round(mean(Commission_Pct, na.rm = TRUE), 1),
    Average_IoU = round(mean(Mean_IoU, na.rm = TRUE), 2), # <--- IoU Summary Added
    .groups = "drop"
  ) %>%
  # Order the factors so Baseline appears above Final for each spacing
  mutate(Phase = factor(Phase, levels = c("Baseline", "Final"))) %>%
  arrange(Spacing, Phase)

# Print the final clean table to the console
print(table_data)

# Optional: Export directly to a CSV to copy-paste into Word/Excel
write.csv(table_data, "C:/Users/jakev/Desktop/ITD_Accuracy_Summary_Table.csv", row.names = FALSE)
