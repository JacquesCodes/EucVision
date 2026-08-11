library(lidR)
library(sf)
library(tictoc)

tic()

# 1. Load the catalog (uses minimal RAM)
ctg <- readLAScatalog("E:/Remote Sensing Media/17. 02 March 2026/03. Point clouds/2.4cm Bottom Sector Demo_02 March 2026_group1_densified_point_cloud.las")

# 2. Load your plot boundaries (e.g., an sf object from a shapefile)
plots <- st_read("E:/Remote Sensing Media/000. Projects/0. Plot boundaries for cropping/Plot 37-40.shp")

# 3. Tell the catalog to write the outputs directly to disk 
# This ensures your RAM never fills up, no matter how many plots you have
opt_select(ctg) <- "xyz"
opt_output_files(ctg) <- "E:/Output/Plot_37-40."

# 4. Extract the plots
# The function will process chunks and save them to the output folder
clip_roi(ctg, plots)

toc()