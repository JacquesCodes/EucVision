# List of all main project directories
main_folders <- c(
  "E:/Remote Sensing Media/01. 25 February 2025 (DJI Mavic)",
  "E:/Remote Sensing Media/02. 01 September 2025 (DJI M300)",
  "E:/Remote Sensing Media/03. 30 October 2025",
  "E:/Remote Sensing Media/04. 07 November 2025",
  "E:/Remote Sensing Media/05. 14 November 2025",
  "E:/Remote Sensing Media/06. 17 November 2025",
  "E:/Remote Sensing Media/07. 28 November 2025",
  "E:/Remote Sensing Media/08. 08 December 2025",
  "E:/Remote Sensing Media/09. 11 December 2025",
  "E:/Remote Sensing Media/10. 22 December 2025",
  "E:/Remote Sensing Media/11. 13 January 2026",
  "E:/Remote Sensing Media/12. 22 January 2026",
  "E:/Remote Sensing Media/13. 29 January 2026",
  "E:/Remote Sensing Media/14. 06 February 2026",
  "E:/Remote Sensing Media/15. 16 February 2026",
  "E:/Remote Sensing Media/16. 23 February 2026",
  "E:/Remote Sensing Media/17. 02 March 2026",
  "E:/Remote Sensing Media/17. 03 March 2026 (Multispectral)",
  "E:/Remote Sensing Media/18. 09 March 2026",
  "E:/Remote Sensing Media/19. 16 March 2026",
  "E:/Remote Sensing Media/20. 23 March 2026",
  "E:/Remote Sensing Media/20. 24 March 2026 (Multispectral)"
)

# Loop through each main flight folder
for (folder in main_folders) {
  
  # 1. Check if the main folder actually exists before proceeding
  if (!dir.exists(folder)) {
    message(paste("Skipping, directory not found:", folder))
    next
  }
  
  message(paste("\nProcessing:", basename(folder)))
  
  # 2. Extract and format the date for the file names
  base_name <- basename(folder)
  clean_date <- sub("^\\d+\\.\\s*", "", base_name)
  clean_date <- sub("\\s*\\(.*?\\)$", "", clean_date)
  file_date_safe <- gsub(" ", "_", clean_date)
  
  # 3. Process the "08. Crown shape file" folder specifically
  dir_08_old <- file.path(folder, "08. Crown shape file")
  dir_08_new <- file.path(folder, "08. Crown Polygons")
  
  if (dir.exists(dir_08_old)) {
    files_08 <- list.files(dir_08_old, pattern = "^All[_ ]Plots\\.", full.names = TRUE)
    for (f in files_08) {
      ext <- tools::file_ext(f)
      new_file_name <- paste0("Crown_Polygons_", file_date_safe, ".", ext)
      file.rename(f, file.path(dir_08_old, new_file_name))
    }
    file.rename(dir_08_old, dir_08_new)
    message("  -> Renamed folder to: 08. Crown Polygons")
  }
  
  # 4. Process the "09. Tree heights" folder specifically
  dir_09_old <- file.path(folder, "09. Tree heights")
  dir_09_new <- file.path(folder, "09. Crown Metrics")
  
  if (dir.exists(dir_09_old)) {
    files_09 <- list.files(dir_09_old, pattern = "^All[_ ]Plots\\.", full.names = TRUE)
    for (f in files_09) {
      ext <- tools::file_ext(f)
      new_file_name <- paste0("Crown_Metrics_", file_date_safe, ".", ext)
      file.rename(f, file.path(dir_09_old, new_file_name))
    }
    file.rename(dir_09_old, dir_09_new)
    message("  -> Renamed folder to: 09. Crown Metrics")
  }
  
  # 5. Global Title Case Standardizer for ALL subdirectories
  # Grab all subdirectories inside the current flight folder
  subdirs <- list.dirs(folder, full.names = TRUE, recursive = FALSE)
  
  for (subdir in subdirs) {
    old_name <- basename(subdir)
    
    # Check if the folder starts with numbers and a dot (e.g., "04. Point clouds clipped")
    if (grepl("^\\d+\\.\\s+.*", old_name)) {
      
      # Split the prefix ("04. ") from the text ("Point clouds clipped")
      prefix <- sub("^(\\d+\\.\\s+).*", "\\1", old_name)
      text_part <- sub("^\\d+\\.\\s+(.*)", "\\1", old_name)
      
      # Capitalize the first letter of every word (Base R Regex Magic)
      title_case_text <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", text_part, perl = TRUE)
      
      # Recombine them
      new_name <- paste0(prefix, title_case_text)
      
      # Rename the folder if the new Title Case name is different from the old name
      if (old_name != new_name) {
        file.rename(subdir, file.path(folder, new_name))
        message(paste("  -> Standardized casing:", old_name, "==>", new_name))
      }
    }
  }
}

message("\nBatch renaming and casing standardization complete!")