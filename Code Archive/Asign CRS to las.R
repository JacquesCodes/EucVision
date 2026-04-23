# Load the necessary spatial libraries
library(lidR)
library(sf)

# 1. Define file path (using the original, un-altered CloudCompare file)
las_original <- "E:/Remote Sensing Media/20. 16 March 2026/03. Point clouds/Bottom_16 March 2026.las"

cat("Loading the raw 28-million point cloud...\n")
las <- readLAS(las_original)

# 2. Stamp the raw, positive points with ESRI:102562
cat("Assigning ESRI:102562 to the raw point cloud...\n")
st_crs(las) <- st_crs("ESRI:102562")

# 3. Transform the points to EPSG:2048 to force the mathematical axis flip
cat("Transforming the points to EPSG:2048 to match the negative orientation...\n")
cat("(This will take a few seconds, but it safely recalculates the LAS header offsets)\n")
las_aligned <- st_transform(las, 2048)

# 4. Save the fixed file
final_file <- "E:/Remote Sensing Media/20. 16 March 2026/03. Point clouds/Bottom_16 March 2026_FinalAligned.las"

cat("Writing the correctly projected file to disk...\n")
writeLAS(las_aligned, final_file)

cat("Success! The point coordinates are permanently aligned.\n")