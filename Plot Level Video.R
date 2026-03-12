# Load necessary libraries
library(terra)
library(sf)
library(magick)
library(stringr)

# --- MEMORY & STORAGE CONFIGURATION ---
base_dir <- "E:/Remote Sensing Media"

# Allocate memory aggressively but safely for 32GB RAM
temp_dir <- file.path(base_dir, "terra_temp")
dir.create(temp_dir, showWarnings = FALSE)
terraOptions(memfrac = 0.75, tempdir = temp_dir)

# --- INPUT CONFIGURATION ---
# Point this to where you saved the uploaded Plot 37-40.shp
shp_path <- "E:/Remote Sensing Media/000. Projects/0. Plot boundaries for cropping/Plot 37-40.shp"
# shp_path <- "E:/Remote Sensing Media/000. Projects/0. Plot boundaries for cropping/Plots shape files/id_19.shp" 
plot_shp <- st_read(shp_path, quiet = TRUE)

# Directory to save the temporary clipped PNGs and final outputs
output_dir <- file.path(base_dir, "TimeLapse_Outputs")
dir.create(output_dir, showWarnings = FALSE)

# Get all dataset folders (assuming the same structure as your previous script)
folders <- list.dirs(base_dir, recursive = FALSE)
exclude_list <- c("00. Dataset template", "000. Projects", "terra_temp", "01. 25 February 2025", "02. 01 September 2025")
dataset_folders <- folders[grepl("^\\d{2}\\.", basename(folders)) & !basename(folders) %in% exclude_list]

# List to store the paths of the generated images
clipped_images <- c()

print("Starting batch crop and image generation...")

# --- BATCH CROP LOOP ---
for (folder_path in dataset_folders) {
  
  folder_name <- basename(folder_path)
  date_str <- str_replace(folder_name, "^\\d{2}\\.\\s*", "")
  ortho_dir <- file.path(folder_path, "01. Orthomosaics")
  
  if (!dir.exists(ortho_dir)) next
  
  all_tifs <- list.files(ortho_dir, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
  
  # Filter OUT "Top" tifs. Keep Single, Bottom, and Cross Bottom.
  valid_tifs <- all_tifs[!grepl("Top", basename(all_tifs), ignore.case = TRUE)]
  
  if (length(valid_tifs) == 0) next
  
  for (tif_file in valid_tifs) {
    print(paste("Processing:", basename(tif_file), "for date:", date_str))
    
    ortho <- rast(tif_file)
    
    # --- CRS & SOUTH AFRICAN AXIS FIX ---
    plot_shp_proj <- plot_shp
    
    # 1. Transform if the CRS strings are technically different
    if (st_crs(plot_shp_proj) != st_crs(ortho)) {
      print("  -> CRS mismatch detected. Transforming shapefile...")
      plot_shp_proj <- suppressWarnings(st_transform(plot_shp_proj, st_crs(ortho)))
    }
    
    # 2. Check for the South African coordinate sign flip
    rast_xmin <- ext(ortho)[1]
    plot_xmin <- st_bbox(plot_shp_proj)[1]
    
    # If the signs are opposite (e.g., raster is negative, plot is positive)
    if (sign(rast_xmin) != sign(plot_xmin)) {
      print("  -> SA Lo axis flip detected! Inverting plot coordinates (* -1)...")
      st_geometry(plot_shp_proj) <- st_geometry(plot_shp_proj) * -1
      st_crs(plot_shp_proj) <- st_crs(ortho) # Re-assign the CRS after math
    }
    
    # Convert to terra format for the crop
    plot_vect_proj <- vect(plot_shp_proj)
    
    # --- CROPPING ---
    cropped_raster <- tryCatch({
      crop(ortho, plot_vect_proj)
    }, error = function(e) {
      NULL
    })
    
    if (is.null(cropped_raster)) {
      print("  -> Plot STILL out of bounds. Skipping.")
      next
    }

    
    # --- EXPORT AS IMAGE FOR ANIMATION ---
    # Save the cropped raster as a high-res PNG
    # We use a naming convention that ensures they sort chronologically
    out_png <- file.path(output_dir, paste0("Plot37_40_", date_str, "_", basename(tif_file), ".png"))
    
    # 1. Make the canvas rectangular (wide) to match the shape of your plot
    png(out_png, width = 3600, height = 1600, res = 300)
    
    # 2. Crush the margins: c(bottom, left, top, right)
    # This leaves just enough space at the top for the title, and almost none elsewhere
    par(mar = c(0.5, 0.5, 2.5, 0.5)) 
    
    # Plot the RGB without a title first so R doesn't auto-space it
    plotRGB(cropped_raster, r = 1, g = 2, b = 3, stretch = "lin", axes = FALSE)
    
    # 3. Add the title manually
    # cex.main = 2.5 makes it bigger. line = 0.5 pulls it down closer to the photo
    title(main = paste("Date:", date_str), cex.main = 2.5, line = 0.5)
    
    dev.off()
    
    clipped_images <- c(clipped_images, out_png)
    
    # Cleanup memory
    rm(ortho, cropped_raster)
    gc()
  }
}

print("Cropping complete! Generating time-lapse and collage...")

library(av)

# --- AV: HIGH-RES MP4 ANIMATION ---
if (length(clipped_images) > 0) {
  
  video_path <- file.path(output_dir, "Plot37_40_TimeLapse.mp4")
  
  # framerate = 1 means 1 image per second
  av_encode_video(clipped_images, output = video_path, framerate = 1)
  
  print(paste("High-Res Time-lapse video saved to:", video_path))
  
} else {
  print("No overlapping images found to animate.")
}

# Clean up temporary directory
unlink(temp_dir, recursive = TRUE)
print("All tasks complete!")