# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: TIME-LAPSE GENERATION & SPATIAL OVERLAY PIPELINE ####
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
# Load required libraries for spatial operations, raster manipulation, and image/video editing
library(terra)
library(sf)
library(magick)
library(stringr)
library(av)

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration & Memory Management ####
# ──────────────────────────────────────────────────────────────────────────────
# --- MEMORY & STORAGE CONFIGURATION ---
base_dir <- "E:/Remote Sensing Media"

# Allocate memory aggressively but safely for 32GB RAM environments
temp_dir <- file.path(base_dir, "terra_temp")
dir.create(temp_dir, showWarnings = FALSE)
terraOptions(memfrac = 0.75, tempdir = temp_dir)

# --- INPUT/OUTPUT CONFIGURATION ---
# 1. Define the target plot boundaries for cropping the viewing extent
shp_path <- "E:/Remote Sensing Media/000. Projects/0. Plot boundaries for cropping/Plot 37-40.shp"
plot_shp <- st_read(shp_path, quiet = TRUE)

# 2. Define the static Normal Plot Boundaries that will be overlaid
normal_plot_path <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LAScatalog Boundaries/Normal Plot Boundaries.shp"
normal_plots <- st_read(normal_plot_path, quiet = TRUE)

# Directory to save the temporary clipped PNGs and the final video outputs
output_dir <- file.path(base_dir, "TimeLapse_Outputs")
dir.create(output_dir, showWarnings = FALSE)

# --- COLOR PALETTES ---
# Color palette mapped to specific Eucalyptus species (Crowns)
species_colors <- c(
  "Cladocalyx"    = "#336998",
  "Grandis"       = "#97dde3",
  "Cloeziana"     = "#ffffff",
  "Urophylla"     = "#e3acff",
  "Grandis clone" = "#ff7da0"
)

# Color palette mapped to the Normal Plot Spacings
plot_colors <- c(
  "1m_Spacing"   = "#eedf13",
  "2m_Spacing"   = "#f10070",
  "3m_Spacing"   = "#15e94e",
  "3x2m_Spacing" = "#404040",
  "5m_Spacing"   = "#ff0004"
)

# ──────────────────────────────────────────────────────────────────────────────
# 3. Batch Processing Setup ####
# ──────────────────────────────────────────────────────────────────────────────
# Fetch all dataset folders and exclude specific utility/baseline directories
folders <- list.dirs(base_dir, recursive = FALSE)
exclude_list <- c("00. Dataset template", 
                  "000. Projects", "terra_temp",
                  "00. Baseline DTM and Plot Cropping",
                  "01. 25 February 2025", 
                  "02. 01 September 2025",
                  "03. 30 October 2025",
                  "17. 03 March 2026 (Multispectral)",
                  "20. 24 March 2026 (Multispectral)")

# Filter for standard dated folders matching the specific nomenclature
dataset_folders <- folders[grepl("^\\d{2}\\.", basename(folders)) & !basename(folders) %in% exclude_list]

# Initialize a list to store the absolute paths of the generated images
clipped_images <- c()

print("Starting batch crop, overlay, and image generation...")

# ──────────────────────────────────────────────────────────────────────────────
# 4. Image Extraction & Overlay Loop ####
# ──────────────────────────────────────────────────────────────────────────────
# Iterate through all valid dataset folders to extract and process imagery
for (folder_path in dataset_folders) {
  
  folder_name <- basename(folder_path)
  date_str <- str_replace(folder_name, "^\\d{2}\\.\\s*", "")
  ortho_dir <- file.path(folder_path, "01. Orthomosaics")
  
  # --- DYNAMIC CROWN POLYGON PATH ---
  crown_dir <- file.path(folder_path, "08. Crown Polygons")
  
  # Search for any shapefile in this specific week's folder
  crown_shp_files <- list.files(crown_dir, pattern = "\\.shp$", full.names = TRUE, ignore.case = TRUE)
  
  # Assign the shapefile if found, otherwise keep an empty string
  if (length(crown_shp_files) > 0) {
    crown_shp_path <- crown_shp_files[1] 
  } else {
    crown_shp_path <- "" 
  }
  
  # Skip iteration if no Orthomosaic folder exists
  if (!dir.exists(ortho_dir)) next
  
  all_tifs <- list.files(ortho_dir, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
  
  # Filter OUT specific unneeded TIFs ("Top" and "Cross Hatch")
  valid_tifs <- all_tifs[!grepl("Top|Cross Hatch", basename(all_tifs), ignore.case = TRUE)]
  
  if (length(valid_tifs) == 0) next
  
  for (tif_file in valid_tifs) {
    print(paste("Processing:", basename(tif_file), "for date:", date_str))
    
    ortho <- rast(tif_file)
    
    # --- CRS & SOUTH AFRICAN AXIS FIX FOR PLOT CROPPING SHAPE ---
    plot_shp_proj <- plot_shp
    
    # Ensure the bounding box shape matches the raster's CRS
    if (st_crs(plot_shp_proj) != st_crs(ortho)) {
      plot_shp_proj <- suppressWarnings(st_transform(plot_shp_proj, st_crs(ortho)))
    }
    
    rast_xmin <- ext(ortho)[1]
    plot_xmin <- st_bbox(plot_shp_proj)[1]
    
    # Fix spatial flipping (Common issue with South African projected coordinates)
    if (sign(rast_xmin) != sign(plot_xmin)) {
      st_geometry(plot_shp_proj) <- st_geometry(plot_shp_proj) * -1
      st_crs(plot_shp_proj) <- st_crs(ortho) 
    }
    
    plot_vect_proj <- vect(plot_shp_proj)
    
    # --- CROPPING RASTER ---
    cropped_raster <- tryCatch({
      crop(ortho, plot_vect_proj)
    }, error = function(e) {
      NULL
    })
    
    # Handle cases where the plot boundary is outside the current raster extent
    if (is.null(cropped_raster)) {
      print("  -> Plot out of bounds. Skipping.")
      next
    }
    
    # --- PROCESS STATIC NORMAL PLOT POLYGONS ---
    normal_plots_cropped <- NULL
    if (exists("normal_plots") && nrow(normal_plots) > 0) {
      normal_plots_proj <- normal_plots
      
      # 1. Match CRS
      if (st_crs(normal_plots_proj) != st_crs(ortho)) {
        normal_plots_proj <- suppressWarnings(st_transform(normal_plots_proj, st_crs(ortho)))
      }
      
      # 2. Check for SA Axis Flip
      normal_xmin <- st_bbox(normal_plots_proj)[1]
      if (sign(rast_xmin) != sign(normal_xmin)) {
        st_geometry(normal_plots_proj) <- st_geometry(normal_plots_proj) * -1
        st_crs(normal_plots_proj) <- st_crs(ortho)
      }
      
      # 3. Intersect with the viewing extent to only draw what is on-screen
      normal_plots_cropped <- suppressWarnings(st_intersection(normal_plots_proj, plot_shp_proj))
    }
    
    # --- PROCESS DYNAMIC CROWN POLYGONS ---
    crowns_cropped <- NULL
    if (file.exists(crown_shp_path)) {
      crowns <- st_read(crown_shp_path, quiet = TRUE)
      
      # 1. Match CRS
      if (st_crs(crowns) != st_crs(ortho)) {
        crowns <- suppressWarnings(st_transform(crowns, st_crs(ortho)))
      }
      
      # 2. Check for SA Axis Flip
      crown_xmin <- st_bbox(crowns)[1]
      if (sign(rast_xmin) != sign(crown_xmin)) {
        st_geometry(crowns) <- st_geometry(crowns) * -1
        st_crs(crowns) <- st_crs(ortho)
      }
      
      # 3. Intersect crowns with the plot boundary
      crowns_cropped <- suppressWarnings(st_intersection(crowns, plot_shp_proj))
    } else {
      print("  -> No crown polygons found for this date. Proceeding without crown overlay.")
    }
    
    # --- EXPORT AS IMAGE FOR ANIMATION ---
    out_png <- file.path(output_dir, paste0("Timelapse_", date_str, "_", basename(tif_file), ".png"))
    
    # Setup high-resolution canvas
    png(out_png, width = 3600, height = 1600, res = 300)
    
    # Crush all margins to 0 so the imagery fills the entire canvas seamlessly
    par(mar = c(0, 0, 0, 0)) 
    
    # Render the RGB raster
    plotRGB(cropped_raster, r = 1, g = 2, b = 3, stretch = "lin", axes = FALSE)
    
    # Overlay 1: Normal Plot Boundaries
    if (!is.null(normal_plots_cropped) && nrow(normal_plots_cropped) > 0) {
      # Dynamically find which column contains the spacing names (e.g., "1m_Spacing")
      color_col <- names(normal_plots_cropped)[sapply(normal_plots_cropped, function(col) any(col %in% names(plot_colors)))][1]
      
      # Assign colors based on the matched column, fallback to white if names don't perfectly match
      if (!is.na(color_col)) {
        border_colors_plot <- plot_colors[as.character(normal_plots_cropped[[color_col]])]
      } else {
        border_colors_plot <- "white" 
      }
      
      # Plot only the borders (col = NA ensures transparent fill, lwd controls thickness)
      plot(st_geometry(normal_plots_cropped), add = TRUE, border = border_colors_plot, col = NA, lwd = 3)
    }
    
    # Overlay 2: Dynamic Crown Polygons
    if (!is.null(crowns_cropped) && nrow(crowns_cropped) > 0) {
      # Map the pre-defined colors to the respective Species
      border_colors_crowns <- species_colors[crowns_cropped$Species]
      
      # Plot only the borders
      plot(st_geometry(crowns_cropped), add = TRUE, border = border_colors_crowns, col = NA, lwd = 3)
    }
    
    dev.off()
    
    # --- ADD OUTLINED TITLE WITH MAGICK ---
    img <- image_read(out_png)
    img <- image_annotate(img, paste("Date:", date_str), 
                          gravity = "north", location = "+0+50", 
                          size = 120, weight = 700, 
                          color = "black", 
                          strokecolor = "white",
                          strokewidth = 3)
    image_write(img, out_png)
    
    # Append the new image path to the master list for video generation
    clipped_images <- c(clipped_images, out_png)
    
    # Execute garbage collection to keep RAM usage stable during long loops
    rm(ortho, cropped_raster, crowns_cropped, normal_plots_cropped)
    gc()
  }
}

print("Cropping and overlays complete! Generating time-lapse...")

# ──────────────────────────────────────────────────────────────────────────────
# 5. Time-Lapse Video Generation ####
# ──────────────────────────────────────────────────────────────────────────────
# --- AV: HIGH-RES MP4 ANIMATION ---
if (length(clipped_images) > 0) {
  
  video_path <- file.path(output_dir, "Timelapse.mp4")
  
  # Compile all generated PNGs into a 1 FPS video
  av_encode_video(clipped_images, output = video_path, framerate = 1)
  
  print(paste("High-Res Time-lapse video saved to:", video_path))
  
} else {
  print("No overlapping images found to animate.")
}

# Clean up the aggressive terra temporary directory
unlink(temp_dir, recursive = TRUE)
print("All tasks complete!")