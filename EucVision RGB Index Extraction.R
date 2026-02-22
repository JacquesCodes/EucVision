# Load necessary libraries
library(terra)
library(sf)
library(exactextractr)
library(dplyr)
library(stringr)
library(tidyr)

# --- MEMORY & STORAGE CONFIGURATION ---
base_dir <- "E:/Remote Sensing Media"

# Force Terra to use your E: drive for temp files to prevent C: drive crashes
# memfrac = 0.5 tells it to only use up to 75% of your RAM before writing to disk
temp_dir <- file.path(base_dir, "terra_temp")
dir.create(temp_dir, showWarnings = FALSE)
terraOptions(memfrac = 0.75, tempdir = temp_dir)

folders <- list.dirs(base_dir, recursive = FALSE)
exclude_list <- c("00. Dataset template", "000. Projects", "terra_temp","14. 06 February 2026")
folders <- folders[!basename(folders) %in% exclude_list]
dataset_folders <- folders[grepl("^\\d{2}\\.", basename(folders))]

# For only one folder at a time
# dataset_folders <- "14. 06 February 2026"

all_extractions <- list()

print("Starting memory-efficient batch extraction...")

# --- BATCH PROCESSING LOOP ---
for (folder_path in dataset_folders) {
  
  folder_name <- basename(folder_path)
  date_str <- str_replace(folder_name, "^\\d{2}\\.\\s*", "")
  
  print("========================================")
  print(paste("Processing Date:", date_str))
  
  ortho_dir <- file.path(folder_path, "01. Orthomosaics")
  shp_dir <- file.path(folder_path, "08. Crown shape file", folder_name)
  
  all_tifs <- list.files(ortho_dir, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
  shp_files <- list.files(shp_dir, pattern = "\\.shp$", full.names = TRUE, ignore.case = TRUE)
  valid_tifs <- all_tifs[!grepl("Cross", basename(all_tifs), ignore.case = TRUE)]
  
  if (length(valid_tifs) == 0 || length(shp_files) == 0) next
  
  for (tif_file in valid_tifs) {
    print(paste("  -> Loading Raster:", basename(tif_file)))
    
    is_top <- grepl("Top", basename(tif_file), ignore.case = TRUE)
    is_bottom <- grepl("Bottom", basename(tif_file), ignore.case = TRUE)
    
    # Load the massive raster lazily (doesn't read pixels into memory yet)
    ortho <- rast(tif_file)
    
    for (shp_file in shp_files) {
      
      plot_num <- as.numeric(str_extract(basename(shp_file), "\\d+"))
      should_process <- TRUE
      
      if (!is.na(plot_num)) {
        if (is_top && plot_num > 21) should_process <- FALSE
        else if (is_bottom && plot_num < 22) should_process <- FALSE
      }
      
      if (!should_process) next
      
      print(paste("    -> Cropping and Extracting for:", basename(shp_file)))
      
      crowns <- st_read(shp_file, quiet = TRUE)
      if (st_crs(crowns) != st_crs(ortho)) crowns <- st_transform(crowns, st_crs(ortho))
      
      # ---------------------------------------------------------
      # THE MAGIC BULLET: Crop the massive TIF to the plot extent
      # ---------------------------------------------------------
      # We use a tryCatch block. If the shapefile is completely outside the raster 
      # (e.g., trying to extract a bottom plot on a top raster), crop() will fail. 
      # We catch that expected error, return NULL, and gracefully skip to the next.
      plot_ortho <- tryCatch({
        crop(ortho, crowns)
      }, error = function(e) {
        NULL
      })
      
      # If it returned NULL, it means no overlap. Skip to next!
      if (is.null(plot_ortho)) {
        print("    -> Plot out of bounds for this TIF. Skipping.")
        next
      }
      
      # Calculate indices ONLY on the tiny cropped raster
      R <- plot_ortho[[1]]
      G <- plot_ortho[[2]]
      B <- plot_ortho[[3]]
      
      VARI <- (G - R) / (G + R - B)
      GLI  <- ((2 * G) - R - B) / ((2 * G) + R + B)
      TGI  <- G - 0.39 * R - 0.61 * B
      
      index_stack <- c(VARI, GLI, TGI)
      names(index_stack) <- c("VARI", "GLI", "TGI")
      
      # Extract values
      ext_vals <- exact_extract(index_stack, crowns, fun = 'mean', progress = FALSE)
      colnames(ext_vals) <- c("Mean_VARI", "Mean_GLI", "Mean_TGI")
      
      crowns_data <- crowns %>%
        st_drop_geometry() %>%
        bind_cols(ext_vals) %>%
        mutate(Date = date_str, Source_Raster = basename(tif_file), Plot_Number = plot_num)
      
      all_extractions[[length(all_extractions) + 1]] <- crowns_data
      
      # Free up memory immediately after this plot is done
      rm(plot_ortho, R, G, B, VARI, GLI, TGI, index_stack, ext_vals, crowns)
      gc() # Force R garbage collection
    }
  }
}

# --- FINALIZE AND EXPORT ---
print("========================================")
print("All folders processed. Merging datasets...")

final_master_df <- bind_rows(all_extractions)
output_path <- file.path(base_dir, "Master_RGB_Indices.csv")
write.csv(final_master_df, output_path, row.names = FALSE)

# Clean up the temporary directory we made
unlink(temp_dir, recursive = TRUE)

print(paste("Success! File saved to:", output_path))