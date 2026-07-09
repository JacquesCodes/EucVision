# Load the necessary spatial libraries
library(lidR)
library(sf)

# Disable scientific notation for clear coordinate reading
options(scipen = 999)

# 1. Define all file paths
las_top_file <- "E:/Remote Sensing Media/18. 03 March 2026/03. Point clouds/Top Sector 50m V2.las"
las_bot_file <- "E:/Remote Sensing Media/18. 03 March 2026/03. Point clouds/Bottom Sector 50m V2.las"
shp_bounds_file <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LidR Boundaries/EucVision LidR Boundaries.shp"
shp_plots_file <- "E:/Remote Sensing Media/20. 16 March 2026/08. Crown shape file/All_Plots.shp"

# 2. Read the data (Headers only for LAS)
cat("Loading spatial metadata...\n")
las_top <- readLASheader(las_top_file)
las_bot <- readLASheader(las_bot_file)
bounds <- st_read(shp_bounds_file, quiet = TRUE)
plots <- st_read(shp_plots_file, quiet = TRUE)

# =========================================
# FIX: Invert Shapefile Coordinates
# =========================================
cat("Inverting Shapefile coordinates to match LAS positive Lo domain...\n")

# Multiply the geometry by -1 to flip the signs across the 0,0 origin
st_geometry(bounds) <- st_geometry(bounds) * -1
st_geometry(plots) <- st_geometry(plots) * -1

# Note: Performing arithmetic on sf geometry sometimes drops or invalidates the CRS string. 
# We re-assign the CRS from the LAS header to ensure they perfectly match moving forward.
st_crs(bounds) <- st_crs(las_top)
st_crs(plots) <- st_crs(las_top)

cat("Fix applied. Re-run bounding box checks.\n")

# 3. CRS Verification
cat("\n=========================================\n")
cat("          CRS CONFORMITY CHECK           \n")
cat("=========================================\n")
cat("Top LAS:        ", st_crs(las_top)$input, "\n")
cat("Bottom LAS:     ", st_crs(las_bot)$input, "\n")
cat("Euc Boundaries: ", st_crs(bounds)$input, "\n")
cat("Crown Plots:    ", st_crs(plots)$input, "\n")
cat("=========================================\n\n")

# 4. Extract bounding boxes and convert to polygons for plotting
poly_top <- st_as_sfc(st_bbox(las_top))
poly_bot <- st_as_sfc(st_bbox(las_bot))

# 5. Visual Verification Plot
cat("Generating visual overlay...\n")

# Setup the base plot area using the larger EucVision boundaries
plot(st_geometry(bounds), border = "darkgreen", lwd = 3, col = "white",
     main = "Master Alignment Verification")

# Add the LAS flight extents (using semi-transparent fills to see overlaps)
plot(poly_top, add = TRUE, border = "blue", lwd = 2, col = adjustcolor("lightblue", alpha.f = 0.3))
plot(poly_bot, add = TRUE, border = "purple", lwd = 2, col = adjustcolor("plum", alpha.f = 0.3))

# Overlay the specific crown plots
plot(st_geometry(plots), add = TRUE, border = "red", lwd = 1.5, col = NA)

# Add a comprehensive legend
legend("topright", 
       legend = c("EucVision Boundaries", "Top LAS Extent", "Bottom LAS Extent", "Crown Plots"), 
       fill = c("white", "lightblue", "plum", NA), 
       border = c("darkgreen", "blue", "purple", "red"),
       bty = "n")

# 6. Raw Coordinate Printout (For final mathematical certainty)
cat("--- Raw Bounding Boxes ---\n")
print(st_bbox(bounds))
print(st_bbox(las_top))
print(st_bbox(las_bot))
print(st_bbox(plots))