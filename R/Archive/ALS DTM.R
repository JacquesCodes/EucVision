library(lidR)
library(RCSF)
library(RMCC)
library(sf)
library(tictoc)
library(dplyr)
library(future)
library(terra)

# Reload the explicitly classified chunks from disk
ctg_classified <- readLAScatalog("C:/Users/jakev/Downloads/cloud0.las")

# --- CRITICAL SITE-LEVEL ENGINE SETTINGS ---
# opt_independent_files(ctg_classified) <- FALSE
# opt_chunk_size(ctg_classified) <- 100          
# opt_chunk_buffer(ctg_classified) <- 2          # TIN only needs a 2m buffer to stitch seams

opt_select(ctg_classified) <- "xyzc"
opt_filter(ctg_classified) <- "-keep_class 2"  # Only use ground points

opt_output_files(ctg_classified) <- "C:/Users/jakev/Downloads/ALS_DTM_30_June_2026"

tic()
print("Generating continuous site DTMs at 0.05m resolution...")
site_dtms <- rasterize_terrain(ctg_classified, 
                               res = 0.05, 
                               algorithm = tin())
toc()
