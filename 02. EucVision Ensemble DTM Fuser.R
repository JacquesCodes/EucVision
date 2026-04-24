library(terra)
library(tictoc)

tic()
print("Fusing multiple temporal DTMs into an Ultimate Baseline DTM...")

# 1. Define the paths to your individual, smoothed baseline DTMs
dtm_1 <- terra::rast("E:/Remote Sensing Media/03. 30 October 2025/05b. Baseline Plot DTMs/Master_Baseline_DTM_Smoothed_30_October_2025.tif")
dtm_2 <- terra::rast("E:/Remote Sensing Media/04. 07 November 2025/05b. Baseline Plot DTMs/Master_Baseline_DTM_Smoothed_07_November_2025.tif")
dtm_3 <- terra::rast("E:/Remote Sensing Media/05. 14 November 2025/05b. Baseline Plot DTMs/Master_Baseline_DTM_Smoothed_14_November_2025.tif")

# 2. Ensure they all perfectly align (resample using dtm_1 as the master reference grid)
if (!ext(dtm_2) == ext(dtm_1)) dtm_2 <- terra::resample(dtm_2, dtm_1, method = "bilinear")
if (!ext(dtm_3) == ext(dtm_1)) dtm_3 <- terra::resample(dtm_3, dtm_1, method = "bilinear")

# 3. Stack them into a multi-layer SpatRaster
dtm_stack <- c(dtm_1, dtm_2, dtm_3)

# 4. Apply the Pixel-wise Temporal MAXIMUM
# Because PTD SfM errors are inherently biased downward (sinkholes), taking the 
# absolute highest value across all 3 dates ensures that if even one flight caught 
# the true solid ground, it overwrites the sinkholes from the other dates.
ultimate_dtm <- terra::app(dtm_stack, fun = "max", na.rm = TRUE)

# 5. Export the God-Tier Baseline DTM
output_path <- "E:/Remote Sensing Media/00. Baseline DTM/Ultimate_Ensemble_Baseline_DTM.tif"

# Ensure the output directory exists
if (!dir.exists(dirname(output_path))) dir.create(dirname(output_path), recursive = TRUE)

terra::writeRaster(ultimate_dtm, filename = output_path, overwrite = TRUE)

print(paste("Ultimate Ensemble DTM saved to:", output_path))
toc()