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
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
library(sf)            # Spatial vector data handling
library(terra)         # Spatial raster data handling
library(dplyr)         # Data wrangling and piping logic
library(lidR)          # Individual Tree Detection algorithms
library(tictoc)        # Script execution timing 
library(stringr)       # String manipulation for dates
library(tidyr)

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration & Batch Management ####
# ──────────────────────────────────────────────────────────────────────────────
base_dir <- "E:/Remote Sensing Media"

# Final output CSV for the complete longitudinal analysis
output_csv <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data Analysis/09. Tree Detection.csv"

# --- RUN CONTROLS ---
target_date_override <- NULL 

# Folders to ignore during the batch processing loop
exclude_list <- c("000. Projects",
                  "00. Baseline DTM",
                  "00. Dataset Template", 
                  "07. December 2025 (TLS)",
                  "17. 02 March 2026 2.4",
                  "17. 02 March 2026 4.8",
                  "17. 02 March 2026 19.2",
                  "17. 02 March 2026 Double Grid",
                  "20. 23 March 2026 0.6cm",
                  "31. 26 June 2026 Oblique",
                  "31. 30 June 2026 (ALS)",
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

print("Starting Batch ITD PiP Validation Pipeline...")

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
  tic("ITD PiP processing complete")
  
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
    
    # 2. Crop the Master CHM to just this plot (with a 3m buffer to avoid edge artifacts)
    plot_extent <- ext(plot_trees) + 3
    chm_cropped <- crop(ctg_chm, plot_extent)
    
    # 3. Dynamic Algorithm Configuration
    # Sets window size to 80% of planting spacing to prevent theoretical overlap
    # Minimum height locked to 30cm to avoid weed/slash noise spikes
    ws_val <- ifelse(is.na(spacing_val) || spacing_val <= 0, 3, spacing_val * 0.8)
    lmf_dynamic <- lmf(ws = ws_val, hmin = 0.3, shape = "circular")
    
    ttops <- tryCatch({
      locate_trees(las = chm_cropped, algorithm = lmf_dynamic)
    }, error = function(e) { NULL })
    
    # 4. Point-in-Polygon (PiP) Validation
    total_detected_points <- 0
    true_positives <- 0
    over_segmented <- 0
    missed_trees <- actual_count # Default to all missed if ITD fails entirely
    
    if (!is.null(ttops) && nrow(ttops) > 0) {
      # Dissolve plot boundary and drop external edge-noise points
      plot_boundary <- st_union(st_make_valid(plot_trees))
      ttops_in_plot <- st_filter(ttops, plot_boundary)
      total_detected_points <- nrow(ttops_in_plot)
      
      # Execute spatial intersection: Counts how many tops fall inside each Geo-SAM polygon
      points_per_poly <- lengths(st_intersects(plot_trees, ttops_in_plot))
      
      # Classify the algorithmic accuracy
      true_positives <- sum(points_per_poly == 1)
      over_segmented <- sum(points_per_poly > 1)
      missed_trees   <- sum(points_per_poly == 0)
    }
    
    # 5. Save the plot's metrics
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
      True_Detection_Rate_Pct = round((true_positives / actual_count) * 100, 2)
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
# 4. Generate Summary Table ####
# ──────────────────────────────────────────────────────────────────────────────

# 1. Load your newly generated PiP results
df_pip <- read.csv("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data Analysis/09. Tree detection.csv")

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
    .groups = "drop"
  ) %>%
  # Order the factors so Baseline appears above Final for each spacing
  mutate(Phase = factor(Phase, levels = c("Baseline", "Final"))) %>%
  arrange(Spacing, Phase)

# Print the final clean table to the console
print(table_data)

# Optional: Export directly to a CSV to copy-paste into Word/Excel
write.csv(table_data, "C:/Users/jakev/Desktop/ITD_Accuracy_Summary_Table.csv", row.names = FALSE)





















