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

for (folder in main_folders) {
  
  if (!dir.exists(folder)) next
  
  base_name <- basename(folder)
  clean_date <- sub("^\\d+\\.\\s*", "", base_name)
  clean_date <- sub("\\s*\\(.*?\\)$", "", clean_date)
  file_date_safe <- gsub(" ", "_", clean_date)
  
  message(paste("\nCleaning up point clouds in:", base_name))
  
  dir_04 <- file.path(folder, "04. Point Clouds Clipped")
  dir_05 <- file.path(folder, "05. Point Clouds Ground Classified")
  dir_06 <- file.path(folder, "06. Point Clouds Normalised")
  dir_07 <- file.path(folder, "07. Canopy Height Models")
  
  rename_plot_files <- function(target_dir, folder_type) {
    if (!dir.exists(target_dir)) return()
    files <- list.files(target_dir, full.names = TRUE)
    
    for (f in files) {
      file_name <- basename(f)
      ext <- tools::file_ext(file_name)
      
      if (grepl("^Plot_\\d+", file_name)) {
        plot_id <- sub("^(Plot_\\d+).*", "\\1", file_name)
        
        # Build the new name with the date occurring ONLY ONCE
        if (folder_type == "04") {
          new_name <- paste0(plot_id, "_", file_date_safe, ".", ext)
        } else if (folder_type == "05") {
          new_name <- paste0(plot_id, "_", file_date_safe, "_classified.", ext)
        } else if (folder_type == "06") {
          new_name <- paste0(plot_id, "_", file_date_safe, "_classified_normalised.", ext)
        } else if (folder_type == "07") {
          new_name <- paste0(plot_id, "_", file_date_safe, "_classified_normalised_chm.", ext)
        }
        
        new_file_path <- file.path(target_dir, new_name)
        
        if (file_name != new_name) {
          file.rename(f, new_file_path)
          message(paste("  -> Cleaned:", new_name))
        }
      }
    }
  }
  
  rename_plot_files(dir_04, "04")
  rename_plot_files(dir_05, "05")
  rename_plot_files(dir_06, "06")
  rename_plot_files(dir_07, "07")
}
message("\nPoint cloud names cleaned and standardized!")