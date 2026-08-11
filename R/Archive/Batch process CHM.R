library(lidR)
library(RCSF)
library(RMCC)
library(gstat)
library(sf)
library(tictoc)
library(geometry)
library(dplyr)
library(future)
library(sp)
library(terra)
library(exactextractr)

# ==============================================================================
# 1. SETUP & STATIC DATA (Runs only once)
# ==============================================================================

# Limit the amount of workers (threads) if you don't have enough RAM. 
# It is better to initialize the plan once outside the loop.
plan(multisession, workers = 8)

# Read in shape files for individual plot boundaries (Does not change per date)
plots_buffered_unsorted <- st_read("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LidR Boundaries/EucVision LidR Boundaries.shp")
plots <- plots_buffered_unsorted[order(plots_buffered_unsorted$id), ]

# List of all date folders to process
date_folders <- c(
  "03. 30 October 2025",
  "04. 07 November 2025",
  "05. 14 November 2025",
  "06. 17 November 2025",
  "07. 28 November 2025",
  "08. 08 December 2025",
  "09. 11 December 2025",
  "10. 22 December 2025",
  "11. 13 January 2026",
  "12. 22 January 2026",
  "13. 29 January 2026",
  "14. 06 February 2026",
  "15. 16 February 2026",
  "16. 23 February 2026",
  "17. 02 March 2026"
)

# ==============================================================================
# 2. BATCH PROCESSING LOOP
# ==============================================================================

for (date_folder in date_folders) {
  
  cat("\n======================================================\n")
  cat("Starting processing for:", date_folder, "\n")
  cat("======================================================\n")

  # Read in tree shape files for height extraction
  trees <- st_read(paste0("E:/Remote Sensing Media/", date_folder, "/08. Crown shape file/All_Plots.shp"))
  
  # Automatically check and transform to EPSG: 2048 if it doesn't match
  if (is.na(st_crs(trees)$epsg) || st_crs(trees)$epsg != 2048) {
    trees <- st_transform(trees, 2048)
    print("Transformed CRS to 2048 successfully.")
  }
  
  trees <- trees %>%
    select(-any_of(c("group_ulid", "N_GM", "id", "N_FG", "N_BG", "BBox")))

  # Rasterize plots ####
  ctg_normalised <- readLAScatalog(paste0("E:/Remote Sensing Media/", date_folder, "/06. Point clouds normalised"))
  
  plan(multisession, workers = 8)
  opt_independent_files(ctg_normalised) <- TRUE
  opt_select(ctg_normalised) <- "xyz"
  
  # Drop ground points and sub-surface noise
  opt_filter(ctg_normalised) <- "-drop_class 2 -drop_z_below 0"
  
  tic("Rasterize canopy time")
  # Write to disk rather than memory:
  opt_output_files(ctg_normalised) <- paste0("E:/Remote Sensing Media/", date_folder, "/07. Canopy Height Models/", "{*}_chm")
  
  # Rasterize canopy with interpolation:
  ctg_chm <- rasterize_canopy(ctg_normalised, 
                              res = 0.05, 
                              algorithm = p2r(na.fill = tin()))
  toc()
  
  # Extract tree heights ####
  tic("Tree height extraction and saving time")
  
  # Calculate metrics using exact_extract (Outputs directly as a vector)
  trees$Tree_Height <- exact_extract(ctg_chm, trees, 'max')
  
  # Save to shapefile (completely overwriting old files to prevent schema errors)
  st_write(trees, paste0("E:/Remote Sensing Media/", date_folder, "/09. Tree heights/All Plots.shp"), delete_dsn = TRUE, quiet = TRUE)
  
  # Save lightweight CSV without the messy spatial geometry text
  write.csv(st_drop_geometry(trees), paste0("E:/Remote Sensing Media/", date_folder, "/09. Tree heights/All Plots.csv"), row.names = FALSE)
  toc()
  
  cat("Finished:", date_folder, "\n")
}

cat("\n======================================================\n")
cat("ALL BATCH PROCESSING COMPLETE!\n")
cat("======================================================\n")