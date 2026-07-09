# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: VISUAL INSPECTION OF GROUND VS CANOPY SPLIT
# ──────────────────────────────────────────────────────────────────────────────
# Description: Loads a single normalized plot, applies the 0.3m threshold, 
#              classifies the points, and opens an interactive 3D viewer.
# ──────────────────────────────────────────────────────────────────────────────

library(lidR)

# 1. Configuration
base_dir <- "E:/Remote Sensing Media"
target_folder <- "30. 30 June 2026 (ALS)"
threshold <- 0.3

# 2. Locate the normalized files
normalised_dir <- file.path(base_dir, target_folder, "06. Point Clouds Normalised")
las_files <- list.files(normalised_dir, pattern = "\\.(las|laz)$", full.names = TRUE)

if (length(las_files) == 0) {
  stop("No normalized point clouds found! Check the directory path.")
}

# 3. Load just the first plot for a quick visual test
# (You can change the index [1] to another number if you want to inspect a different plot)
test_file <- las_files[1]
print(paste("Loading for inspection:", basename(test_file)))

las <- readLAS(test_file)

# 4. Apply Binary Classification based on the threshold
# Standard LAS codes: 2 = Ground, 5 = High Vegetation
# First, reset all points to unclassified (code 1) just to be safe
las@data$Classification <- 1L

# Assign Ground (<= 0.3m)
las@data$Classification[las@data$Z <= threshold] <- 2L

# Assign Canopy (> 0.3m)
las@data$Classification[las@data$Z > threshold] <- 5L

# 5. Interactive 3D Plot
print(paste("Rendering 3D Viewer. Threshold set at", threshold, "meters."))
print("You can rotate, zoom, and pan in the pop-up window.")

# Plot colored by our new classification
plot(las, color = "Classification", bg = "white", legend = TRUE)