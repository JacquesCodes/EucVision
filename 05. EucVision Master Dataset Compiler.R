# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: DATA EXTRACTION, MERGING, & CLEANING PIPELINE ####
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────
# Description: Automates the extraction, consolidation, and spatial alignment of 
#              temporal crown metrics across all UAV flight datasets. It seamlessly 
#              joins drone-derived data with ground-truth field measurements and 
#              dynamically pads the dataset using a static master baseline template. 
#              This ensures all dead and unmeasured trees are explicitly tracked 
#              across time with properly aligned NA values and mortality dates. 
#              Finally, the pipeline applies statistical outlier filtering to remove 
#              anomalous height spikes and exports a clean, chronologically sorted 
#              Master Dataset for downstream longitudinal analysis.
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
library(dplyr)
library(readr)
library(stringr)

Sys.setlocale("LC_TIME", "C")

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration & Path Management ####
# ──────────────────────────────────────────────────────────────────────────────
src_base_dir <- "E:/Remote Sensing Media"
dest_backup_dir <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/05. Crown Metrics"

dest_master_csv <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data analysis/01. Master Dataset.csv"
field_measurements_csv <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data analysis/02. Field Measurements.csv"

# --- ADDED: Path to your Master Template to track the Dead Trees ---
template_csv <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data analysis/00. Dataset template.csv"

# Load the master template once into memory
if (file.exists(template_csv)) {
  master_template <- read_csv(template_csv, show_col_types = FALSE) %>%
    mutate(Tree = round(as.numeric(Tree), 2))
} else {
  stop("CRITICAL ERROR: '00. Dataset template.csv' not found. Cannot pad dead trees.")
}

# ──────────────────────────────────────────────────────────────────────────────
# 3. Run Script & Extract UAV Data ####
# ──────────────────────────────────────────────────────────────────────────────
# --- EXCLUDE LIST ---
# Folders to ignore during the batch processing loop
exclude_list <- c("000. Projects",
                  "00. Baseline DTM",
                  "00. Dataset Template", 
                  "07. December 2025 (TLS)",
                  "17. 02 March 2026 0.6",
                  "17. 02 March 2026 2.4",
                  "17. 02 March 2026 4.8",
                  "17. 02 March 2026 19.2",
                  "17. 03 March 2026 (Multispectral)",
                  "20. 23 March 2026 0.6cm",
                  "20. 24 March 2026 (Multispectral)")

# Scan the base directory and filter for valid date folders
folders <- list.dirs(src_base_dir, recursive = FALSE)
main_folders <- folders[grepl("^\\d{2}\\.", basename(folders)) & !basename(folders) %in% exclude_list]

csv_list <- list()

cat("\n--- Starting Teams backup and data extraction ---\n")

for (folder in main_folders) {
  folder_name <- basename(folder)
  
  date_match <- str_extract(folder_name, "\\d{2} [A-Za-z]+ \\d{4}")
  formatted_date <- NA
  if (!is.na(date_match)) {
    parsed_date <- as.Date(date_match, format="%d %B %Y")
    formatted_date <- format(parsed_date, "%d-%m-%Y") 
  }
  
  crown_metrics_path <- file.path(folder, "09. Crown Metrics")
  
  if (dir.exists(crown_metrics_path)) {
    files_to_copy <- list.files(crown_metrics_path, 
                                pattern = "\\.(shp|shx|dbf|prj|csv)$", 
                                full.names = TRUE, ignore.case = TRUE)
    
    if (length(files_to_copy) > 0) {
      current_dest_dir <- file.path(dest_backup_dir, folder_name)
      if (!dir.exists(current_dest_dir)) dir.create(current_dest_dir, recursive = TRUE)
      
      file.copy(from = files_to_copy, to = current_dest_dir, overwrite = TRUE)
      
      csv_file <- files_to_copy[grepl("\\.csv$", files_to_copy, ignore.case = TRUE)]
      
      if (length(csv_file) == 1) {
        # Load drone data and ensure Tree column matches template precision
        temp_df <- read_csv(csv_file, show_col_types = FALSE) %>%
          mutate(Tree = round(as.numeric(Tree), 2))
        
        # Clean Columns
        if ("Cmprtmn" %in% names(temp_df)) temp_df <- rename(temp_df, Compartment = Cmprtmn)
        
        if ("Area_m2" %in% names(temp_df) && !"Area" %in% names(temp_df)) {
          temp_df <- rename(temp_df, Crown_Area = Area_m2)
        } else if ("Area" %in% names(temp_df) && !"Area_m2" %in% names(temp_df)) {
          temp_df <- rename(temp_df, Crown_Area = Area)
        } else if ("Area_m2" %in% names(temp_df) && "Area" %in% names(temp_df)) {
          temp_df <- mutate(temp_df, Crown_Area = coalesce(Area_m2, Area))
        }
        
        cols_to_keep <- c("Compartment", "Line", "Plot", "Culture", "Spacing", 
                          "Species", "Tree", "Crown_Area", "Tree_Height")
        temp_df <- select(temp_df, any_of(cols_to_keep))
        
        # --- INJECT DEAD TREES: Left join the master template with the drone data ---
        # Alive trees get their drone metrics. Missing/Dead trees get NA.
        temp_df <- master_template %>%
          left_join(temp_df, by = c("Compartment", "Line", "Plot", "Culture", "Spacing", "Species", "Tree"))
        
        # Assign the current flight date to all rows (alive and dead)
        temp_df$Date <- formatted_date
        
        csv_list[[folder_name]] <- temp_df
        cat(paste("SUCCESS: Processed and padded ->", folder_name, "\n"))
      }
    }
  }
}

# ──────────────────────────────────────────────────────────────────────────────
# 4. Merge, Full Join External Data, & Dynamic Cleaning ####
# ──────────────────────────────────────────────────────────────────────────────
if (length(csv_list) > 0) {
  cat("\nMerging UAV CSV files...\n")
  master_dataset <- bind_rows(csv_list)
  
  if (file.exists(field_measurements_csv)) {
    cat("Found '02. Field Measurements.csv'. Padding field dates and running FULL JOIN...\n")
    other_data <- read_csv(field_measurements_csv, show_col_types = FALSE) %>%
      mutate(Tree = round(as.numeric(Tree), 2))
    
    # --- PAD FIELD DATES: Ensure field-only dates also include the dead trees ---
    padded_other_list <- list()
    for (f_date in unique(other_data$Date)) {
      df_date <- other_data %>% filter(Date == f_date)
      
      # Remove Death_Date if it accidentally exists in field data to prevent .x/.y duplication
      df_date <- select(df_date, -any_of("Death_Date"))
      
      # Pad the field date with the template
      df_padded <- master_template %>%
        left_join(df_date, by = c("Compartment", "Line", "Plot", "Culture", "Spacing", "Species", "Tree"))
      df_padded$Date <- f_date
      
      padded_other_list[[f_date]] <- df_padded
    }
    other_data <- bind_rows(padded_other_list)
    
    # Safely rename columns in the external dataset
    if ("Tree_Height" %in% names(other_data)) {
      other_data <- other_data %>% rename(Tree_Height_other = Tree_Height)
    }
    if ("Crown_Area" %in% names(other_data)) {
      other_data <- other_data %>% rename(Crown_Area_other = Crown_Area)
    }
    
    # Include 'Death_Date' in the join to perfectly align the two padded datasets
    master_dataset <- master_dataset %>%
      full_join(other_data, by = c("Compartment", "Line", "Plot", "Culture", "Spacing", "Species", "Tree", "Date", "Death_Date"))
    
    # Neatly merge the external heights
    if ("Tree_Height_other" %in% names(master_dataset)) {
      master_dataset <- suppressWarnings(
        master_dataset %>%
          mutate(
            Tree_Height = as.numeric(Tree_Height),
            Tree_Height_other = as.numeric(Tree_Height_other),
            Tree_Height = coalesce(Tree_Height, Tree_Height_other)
          ) %>%
          select(-Tree_Height_other)
      )
    }
    
    # Neatly merge the external Crown_Area
    if ("Crown_Area_other" %in% names(master_dataset)) {
      master_dataset <- suppressWarnings(
        master_dataset %>%
          mutate(
            Crown_Area = as.numeric(Crown_Area),
            Crown_Area_other = as.numeric(Crown_Area_other),
            Crown_Area = coalesce(Crown_Area, Crown_Area_other)
          ) %>%
          select(-Crown_Area_other)
      )
    }
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # 5. Outlier Filtering & Data Type Standardization ####
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
  
  # Dynamic Outlier Filtering
  master_dataset <- master_dataset %>%
    group_by(Date) %>%
    mutate(
      # Safely handle dates that have purely field data and zero drone heights
      Flight_99th = ifelse(all(is.na(Tree_Height)), NA, quantile(Tree_Height, probs = 0.99, na.rm = TRUE)),
      Tree_Height = ifelse(!is.na(Flight_99th) & Tree_Height > (Flight_99th + 5), NA, Tree_Height)
    ) %>%
    select(-Flight_99th) %>% 
    ungroup()
  
  # ────────────────────────────────────────────────────────────────────────────
  # 6. Chronological & Spatial Sorting & Export ####
  # ────────────────────────────────────────────────────────────────────────────
  cat("Organizing dataset chronologically...\n")
  master_dataset <- master_dataset %>%
    mutate(Temp_Date = as.Date(Date, format="%d-%m-%Y")) %>% 
    arrange(Temp_Date, Compartment, Line, Plot, Tree) %>%
    select(-Temp_Date)
  
  # --- EXPORT ---
  master_csv_folder <- dirname(dest_master_csv)
  if (!dir.exists(master_csv_folder)) dir.create(master_csv_folder, recursive = TRUE)
  
  tryCatch({
    write_csv(master_dataset, dest_master_csv)
    cat(paste("\nDONE: Master dataset updated and sorted chronologically with", nrow(master_dataset), "rows.\n"))
  }, error = function(e) {
    cat("\nERROR: Could not overwrite the master CSV. Is it open in Excel?\n")
  })
  
} else {
  cat("\nNo CSV files were found in the E: drive to merge.\n")
}