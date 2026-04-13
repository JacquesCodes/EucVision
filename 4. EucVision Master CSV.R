# 1. Install missing packages if necessary
if (!require("dplyr")) install.packages("dplyr")
if (!require("readr")) install.packages("readr")
if (!require("stringr")) install.packages("stringr")

library(dplyr)
library(readr)
library(stringr)

# Force R to use English for date parsing
Sys.setlocale("LC_TIME", "C")

# =====================================================================
# 1. DEFINE DIRECTORIES
# =====================================================================
src_base_dir <- "E:/Remote Sensing Media"
dest_backup_dir <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/05. Crown Metrics"

dest_master_csv <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data analysis/01. Main dataset CSV.csv"
other_dataset_csv <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data analysis/01. Other dataset.csv"

# =====================================================================
# 2. RUN SCRIPT & EXTRACT UAV DATA
# =====================================================================
main_folders <- list.dirs(src_base_dir, recursive = FALSE)
csv_list <- list()

cat("\n--- Starting Teams backup and data extraction ---\n")

for (folder in main_folders) {
  folder_name <- basename(folder)
  
  date_match <- str_extract(folder_name, "\\d{2} [A-Za-z]+ \\d{4}")
  formatted_date <- NA
  if (!is.na(date_match)) {
    parsed_date <- as.Date(date_match, format="%d %B %Y")
    formatted_date <- format(parsed_date, "%d-%b-%y")
  }
  
  crown_metrics_path <- file.path(folder, "09. Crown Metrics")
  
  if (dir.exists(crown_metrics_path)) {
    files_to_copy <- list.files(crown_metrics_path, 
                                pattern = "\\.(shp|shx|dbf|prj|csv)$", 
                                full.names = TRUE, 
                                ignore.case = TRUE)
    
    if (length(files_to_copy) > 0) {
      current_dest_dir <- file.path(dest_backup_dir, folder_name)
      if (!dir.exists(current_dest_dir)) dir.create(current_dest_dir, recursive = TRUE)
      
      file.copy(from = files_to_copy, to = current_dest_dir, overwrite = TRUE)
      
      csv_file <- files_to_copy[grepl("\\.csv$", files_to_copy, ignore.case = TRUE)]
      
      if (length(csv_file) == 1) {
        temp_df <- read_csv(csv_file, show_col_types = FALSE)
        
        # Clean Columns
        if ("Cmprtmn" %in% names(temp_df)) temp_df <- rename(temp_df, Compartment = Cmprtmn)
        
        if ("Area_m2" %in% names(temp_df) && !"Area" %in% names(temp_df)) {
          temp_df <- rename(temp_df, Crown_Area = Area_m2)
        } else if ("Area" %in% names(temp_df) && !"Area_m2" %in% names(temp_df)) {
          temp_df <- rename(temp_df, Crown_Area = Area)
        } else if ("Area_m2" %in% names(temp_df) && "Area" %in% names(temp_df)) {
          temp_df <- mutate(temp_df, Crown_Area = coalesce(Area_m2, Area))
        }
        
        temp_df$Date <- formatted_date
        
        cols_to_keep <- c("Compartment", "Line", "Plot", "Culture", "Spacing", 
                          "Species", "Tree", "Date", "Crown_Area", "Tree_Height")
        temp_df <- select(temp_df, any_of(cols_to_keep))
        
        csv_list[[folder_name]] <- temp_df
        cat(paste("SUCCESS: Processed ->", folder_name, "\n"))
      }
    }
  }
}

# =====================================================================
# 3. MERGE, FULL JOIN EXTERNAL DATA, & DYNAMIC CLEANING
# =====================================================================
if (length(csv_list) > 0) {
  cat("\nMerging UAV CSV files...\n")
  master_dataset <- bind_rows(csv_list)
  
  if (file.exists(other_dataset_csv)) {
    cat("Found '01. Other dataset.csv'. Running FULL JOIN...\n")
    other_data <- read_csv(other_dataset_csv, show_col_types = FALSE)
    
    # Safely rename columns in the other dataset to prevent .x and .y duplicates
    if ("Tree_Height" %in% names(other_data)) {
      other_data <- other_data %>% rename(Tree_Height_other = Tree_Height)
    }
    if ("Crown_Area" %in% names(other_data)) {
      other_data <- other_data %>% rename(Crown_Area_other = Crown_Area)
    }
    
    master_dataset <- master_dataset %>%
      full_join(other_data, by = c("Compartment", "Line", "Plot", "Culture", "Spacing", "Species", "Tree", "Date"))
    
    # Neatly merge the Photo/Other heights into the main Tree_Height column
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
    
    # Neatly merge the manual Crown_Area if it exists
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
  
  # --- OUTLIER CLEANING & DATA TYPE FIXING ---
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
  
  # 2. DYNAMIC OUTLIER FILTERING
  master_dataset <- master_dataset %>%
    group_by(Date) %>%
    mutate(
      Flight_99th = quantile(Tree_Height, probs = 0.99, na.rm = TRUE),
      Tree_Height = ifelse(Tree_Height > (Flight_99th + 5), NA, Tree_Height)
    ) %>%
    select(-Flight_99th) %>% 
    ungroup() 
  
  # --- CHRONOLOGICAL & SPATIAL SORTING ---
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
    cat(paste("\nDONE: Master dataset updated and sorted chronologically with", nrow(master_dataset), "rows.\n"))
  }, error = function(e) {
    cat("\nERROR: Could not overwrite the master CSV. Is it open in Excel?\n")
  })
  
} else {
  cat("\nNo CSV files were found in the E: drive to merge.\n")
}