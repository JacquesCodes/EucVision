# 1. Install missing packages if necessary
if (!require("dplyr")) install.packages("dplyr")
if (!require("readr")) install.packages("readr")
if (!require("stringr")) install.packages("stringr")

library(dplyr)
library(readr)
library(stringr)

# Force R to use English for date parsing (fixes South African locale issues)
Sys.setlocale("LC_TIME", "C")

# =====================================================================
# 1. DEFINE DIRECTORIES
# =====================================================================
src_base_dir <- "E:/Remote Sensing Media"
dest_backup_dir <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/05. Crown Metrics"
dest_master_csv <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data analysis/01. Main dataset CSV.csv"

# =====================================================================
# 2. RUN SCRIPT
# =====================================================================
main_folders <- list.dirs(src_base_dir, recursive = FALSE)
csv_list <- list()

cat("\n--- Starting Teams backup and data merge ---\n")

for (folder in main_folders) {
  folder_name <- basename(folder)
  
  # Extract the date from the folder name
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
        
        # Clean Compartment
        if ("Cmprtmn" %in% names(temp_df)) temp_df <- rename(temp_df, Compartment = Cmprtmn)
        
        # Clean Area
        if ("Area_m2" %in% names(temp_df) && !"Area" %in% names(temp_df)) {
          temp_df <- rename(temp_df, Crown_Area = Area_m2)
        } else if ("Area" %in% names(temp_df) && !"Area_m2" %in% names(temp_df)) {
          temp_df <- rename(temp_df, Crown_Area = Area)
        } else if ("Area_m2" %in% names(temp_df) && "Area" %in% names(temp_df)) {
          temp_df <- mutate(temp_df, Crown_Area = coalesce(Area_m2, Area))
        }
        
        # Apply Date
        temp_df$Date <- formatted_date
        
        # Filter Columns strictly matching your Template
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
# 3. EXPORT MASTER DATASET
# =====================================================================
if (length(csv_list) > 0) {
  cat("\nMerging CSV files...\n")
  master_dataset <- bind_rows(csv_list)
  
  master_csv_folder <- dirname(dest_master_csv)
  if (!dir.exists(master_csv_folder)) dir.create(master_csv_folder, recursive = TRUE)
  
  # Try Catch block to check if the file is locked by Excel
  tryCatch({
    write_csv(master_dataset, dest_master_csv)
    cat(paste("\nDONE: Master dataset updated with", nrow(master_dataset), "rows.\n"))
    cat("Please check the '01. Main dataset CSV.csv' file now!\n")
  }, error = function(e) {
    cat("\nERROR: Could not overwrite the master CSV.\n")
    cat("Is '01. Main dataset CSV.csv' currently open in Excel? If yes, please close it and run again.\n")
    cat("Technical error message:", conditionMessage(e), "\n")
  })
  
} else {
  cat("\nNo CSV files were found in the E: drive to merge.\n")
}