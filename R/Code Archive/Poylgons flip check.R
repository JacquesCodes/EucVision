# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: SPATIAL GRID REVERSAL QA/QC PIPELINE (METRICS VERSION)
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/) 
# ──────────────────────────────────────────────────────────────────────────────
# Description: This script loops through temporal SfM datasets and compares the 
#              spatial geometry of each tree polygon against its original 
#              baseline position ("01. 25 February 2025"). If a tree's ID does 
#              not spatially intersect its original planting location, it is 
#              flagged. Plots with a high percentage of non-touching trees are 
#              outputted as confirmed grid reversals (flipped plots).
#              *UPDATED to use final Crown_Metrics shapefiles.
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
library(sf)              # Spatial vector data handling
library(dplyr)           # Data wrangling and piping logic
library(stringr)         # String manipulation
library(tictoc)          # Script execution timing 

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration & Baseline Initialization ####
# ──────────────────────────────────────────────────────────────────────────────
base_dir <- "E:/Remote Sensing Media"
# baseline_date_folder <- "01. 25 February 2025"
baseline_date_folder <- "02. 01 September 2025"

# Safe filename string for baseline
baseline_date_safe <- gsub(" ", "_", sub("^\\d+\\.\\s*", "", baseline_date_folder))

# UPDATED: Point to Crown Metrics instead of Polygons
baseline_shp_path <- file.path(base_dir, baseline_date_folder, "09. Crown Metrics", 
                               paste0("Crown_Metrics_", baseline_date_safe, ".shp"))

if (!file.exists(baseline_shp_path)) {
  stop("CRITICAL ERROR: Baseline shapefile not found! Run master pipeline for Feb 2025 first.")
}

print("Loading Master Baseline Metrics (25 Feb 2025)...")

# Load baseline and create a Unique ID (UID) for exact matching
baseline_trees <- st_read(baseline_shp_path, quiet = TRUE) %>%
  # Safely handle standard DBF column name truncations (like "Cmprtmn" or "Compartmen")
  rename(any_of(c("Compartment" = "Cmprtmn", "Compartment" = "Compartmen"))) %>% 
  mutate(
    Tree = round(Tree, 2),
    UID = paste(Compartment, Plot, Tree, sep = "_")
  )

# --- EXCLUDE LIST ---
exclude_list <- c("000. Projects",
                  "00. Baseline DTM",
                  "00. Dataset Template", 
                  "01. 25 February 2025", # Exclude baseline from being checked against itself
                  "07. December 2025 (TLS)",
                  "17. 03 March 2026 (Multispectral)",
                  "20. 24 March 2026 (Multispectral)")

folders <- list.dirs(base_dir, recursive = FALSE)
dataset_folders <- folders[grepl("^\\d{2}\\.", basename(folders)) & !basename(folders) %in% exclude_list]

# Output dataframe to store all flagged plots across the whole season
all_suspect_plots <- data.frame()

# ──────────────────────────────────────────────────────────────────────────────
# 3. SPATIAL INTERSECTION BATCH LOOP ####
# ──────────────────────────────────────────────────────────────────────────────
tic("Total Spatial QA/QC Time")

for (folder_path in dataset_folders) {
  date_folder <- basename(folder_path)
  file_date_safe <- gsub(" ", "_", sub("^\\d+\\.\\s*", "", date_folder))
  
  # UPDATED: Point to Crown Metrics
  target_shp_path <- file.path(folder_path, "09. Crown Metrics", 
                               paste0("Crown_Metrics_", file_date_safe, ".shp"))
  
  if (!file.exists(target_shp_path)) {
    next # Skip if this date hasn't been processed by the main pipeline yet
  }
  
  # Load target date polygons
  target_trees <- st_read(target_shp_path, quiet = TRUE) %>%
    rename(any_of(c("Compartment" = "Cmprtmn", "Compartment" = "Compartmen"))) %>%
    mutate(
      Tree = round(Tree, 2),
      UID = paste(Compartment, Plot, Tree, sep = "_")
    )
  
  # Find intersecting UIDs (trees that exist in both baseline and target date)
  common_uids <- intersect(baseline_trees$UID, target_trees$UID)
  
  if(length(common_uids) == 0) next
  
  # Subset and align the dataframes so row 1 in baseline exactly matches row 1 in target
  base_sub <- baseline_trees[match(common_uids, baseline_trees$UID), ]
  targ_sub <- target_trees[match(common_uids, target_trees$UID), ]
  
  # Calculate pairwise distance element-by-element 
  # (st_distance == 0 means the polygons are physically touching or overlapping)
  distances <- st_distance(base_sub, targ_sub, by_element = TRUE)
  
  # Add a tiny 0.1m tolerance just in case minor edge smoothing pulled a crown slightly away
  targ_sub$Touching_Baseline <- as.numeric(distances) <= 0.1
  
  # Aggregate stats per plot
  plot_stats <- st_drop_geometry(targ_sub) %>%
    group_by(Compartment, Plot) %>%
    summarise(
      Total_Trees = n(),
      Trees_Moved = sum(Touching_Baseline == FALSE),
      Pct_Flipped = (Trees_Moved / Total_Trees) * 100,
      .groups = 'drop'
    ) %>%
    # If more than 25% of the trees in a plot aren't touching their original spot, it's flipped!
    filter(Pct_Flipped >= 25) %>%
    mutate(Date = file_date_safe)
  
  if(nrow(plot_stats) > 0) {
    all_suspect_plots <- bind_rows(all_suspect_plots, plot_stats)
    print(paste("-> ⚠️ Found", nrow(plot_stats), "flipped plot(s) in:", date_folder))
  }
}

toc()

# ──────────────────────────────────────────────────────────────────────────────
# 4. Final Output ####
# ──────────────────────────────────────────────────────────────────────────────
if(nrow(all_suspect_plots) > 0) {
  print("================================================================")
  print("CRITICAL SPATIAL GRID REVERSALS DETECTED:")
  print("================================================================")
  
  # Sort by worst offenders
  final_report <- all_suspect_plots %>% 
    arrange(desc(Pct_Flipped), Compartment, Plot) %>%
    select(Date, Compartment, Plot, Total_Trees, Trees_Moved, Pct_Flipped)
  
  print(final_report)
  
  # Optional: Write to CSV so you can tick them off in QGIS
  write.csv(final_report, file.path(base_dir, "Flipped_Plots_Report.csv"), row.names = FALSE)
  print("-> Saved report to E:/Remote Sensing Media/Flipped_Plots_Report.csv")
} else {
  print("All plots spatially aligned with baseline. No grid reversals found!")
}