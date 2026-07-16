# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: TEMPERATURE HEAT MAP TIME-LAPSE (ALL PLOTS) ####
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Description: Generates a heat map time-lapse of daily maximum soil temperatures
#              for the entire Lourensford site, overlaid on the merged orthomosaic.
# ──────────────────────────────────────────────────────────────────────────────

# 1. Setup and Imports ####
library(terra)
library(sf)
library(magick)
library(stringr)
library(av)
library(dplyr)       
library(lubridate)   
library(ggplot2)     
library(tidyterra)   

# --- CONFIGURATION ---
base_dir <- "E:/Remote Sensing Media"
output_dir <- file.path(base_dir, "TimeLapse_Outputs", "Temperature_Heatmaps_All")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Load and Prepare Spatial Data ####
print("Loading spatial data for the entire site...")

# Load static Normal Plot Boundaries
normal_plot_path <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. QGIS Shapefiles/02. Plot Boundaries/Normal_Plot_Boundaries_EPSG_2048.shp"
target_plots <- st_read(normal_plot_path, quiet = TRUE)

# Load the newly merged Base TIFF
merged_tif_path <- file.path(base_dir, "20. 23 March 2026/01. Orthomosaics/Merged_23_March_2026.tiff")
base_raster <- rast(merged_tif_path)

# Ensure CRS matches
if (st_crs(target_plots) != st_crs(base_raster)) {
  target_plots <- st_transform(target_plots, crs(base_raster))
}

# --- EXTENT MODIFICATION ---
# Buffer the plots by 5 meters to create a slightly larger viewing area
expanded_extent <- st_buffer(target_plots, dist = 5)

# Crop the heavy raster to this new expanded extent
base_raster_cropped <- crop(base_raster, ext(expanded_extent))

# 3. Load and Prepare Wilo Temperature Data ####
print("Processing Wilo CSV data...")

wilo_csv_path <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/09. Wilo/Wilo_Cleaned.csv"
wilo_data <- read.csv(wilo_csv_path)

# Clean and summarize the data for Daily Maximum Soil Temperature 
daily_max_temp <- wilo_data %>%
  mutate(
    Date = make_date(year, month, day),
    Plot_ID = as.numeric(str_extract(device_id, "\\d+")) 
  ) %>%
  # Filter out hardware glitches and keep all plot data from Sept 1 onwards
  filter(Date >= as.Date("2025-09-01"), 
         SoilTemp > -10 & SoilTemp < 60) %>% 
  group_by(Date, Plot_ID) %>%
  summarize(Max_Temp = max(SoilTemp, na.rm = TRUE), .groups = "drop")

# Get unique dates to loop through
dates_to_process <- sort(unique(daily_max_temp$Date))
clipped_images <- c()

# 4. Rendering Loop ####
print("Starting frame generation...")

for (d in dates_to_process) {
  d_date <- as.Date(d, origin = "1970-01-01")
  date_str <- format(d_date, "%Y-%m-%d")
  print(paste("Rendering frame for:", date_str))
  
  # Filter data for the specific day
  day_data <- daily_max_temp %>% filter(Date == d_date)
  
  # Join temperature data to the spatial plot boundaries 
  plot_sf_daily <- target_plots %>%
    left_join(day_data, by = c("id" = "Plot_ID"))
  
  out_png <- file.path(output_dir, paste0("TempMap_All_", date_str, ".png"))
  
  # Create the plot
  p <- ggplot() +
    # Draw the cropped RGB raster base map
    geom_spatraster_rgb(data = base_raster_cropped, maxcell = 5e6) + 
    # Draw the plot boundaries, filled by the Soil Temp (Max_Temp)
    geom_sf(data = plot_sf_daily, aes(fill = Max_Temp), color = "white", linewidth = 0.8, alpha = 0.75) +
    # Apply the constant color scale (5 to 40 degrees) with dark gray for offline sensors
    scale_fill_viridis_c(
      option = "inferno", 
      limits = c(5, 40), 
      na.value = "#303030",
      name = "Max Soil Temp (°C)"
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(t = 20, b = 10))
    ) +
    ggtitle(paste("Daily Maximum Soil Temperature:", date_str))
  
  # Save the frame 
  # Note: You may want to tweak the width/height depending on the full site's aspect ratio
  ggsave(out_png, plot = p, width = 12, height = 10, dpi = 300, bg = "white")
  
  clipped_images <- c(clipped_images, out_png)
}

# 5. Video Compilation ####
print("Frames complete! Compiling MP4...")

if (length(clipped_images) > 0) {
  video_path <- file.path(output_dir, "All_Plots_Max_Soil_Temp.mp4")
  av_encode_video(clipped_images, output = video_path, framerate = 5) 
  print(paste("Video successfully saved to:", video_path))
} else {
  print("No frames were generated. Check if the CSV date/sensor filtering returned data.")
}

print("Pipeline complete! 🌳")