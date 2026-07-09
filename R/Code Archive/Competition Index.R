# 1. Install and load necessary spatial packages
# install.packages(c("sf", "dplyr")) # Run this if you haven't installed them yet
library(sf)
library(dplyr)

# Change this single variable for each new batch!
date_folder <- "17. 02 March 2026"

# 2. Read the Shapefile
# Make sure All Plots.shp, .shx, .dbf, and .prj are in your working directory
trees_sf <- st_read(paste0("E:/Remote Sensing Media/",date_folder,"/09. Tree heights/All Plots.shp"))

# 3. BULLETPROOFING: Force variables to numeric
trees_sf$Area_Num <- as.numeric(trees_sf$Area)
trees_sf$Height_Num <- as.numeric(trees_sf$Tr_Hght)
# Clean the spacing column just in case it has text like "2m"
trees_sf$Spacing_Num <- as.numeric(gsub("[^0-9.]", "", as.character(trees_sf$Spacing)))

# Convert tree polygons to points (centroids) to make distance and buffer math much faster
trees_pts <- st_centroid(trees_sf)

# --- GENERATE PLOT BOUNDARIES FOR EDGE CORRECTION ---
# This automatically draws a tight "fence" (convex hull) around the outermost trees of each plot
plot_boundaries <- trees_pts %>%
  group_by(Plt_shp) %>%
  summarise(geometry = st_convex_hull(st_union(geometry)))

# 4. Calculate a pairwise distance matrix between all trees
dist_matrix <- st_distance(trees_pts)
dist_matrix <- as.numeric(dist_matrix) 
dim(dist_matrix) <- c(nrow(trees_sf), nrow(trees_sf)) 

# 5. Initialize empty vectors
raw_ci_values <- numeric(nrow(trees_sf))
adj_ci_values <- numeric(nrow(trees_sf))
prop_inside_values <- numeric(nrow(trees_sf))

print("Calculating Competition Index with Edge Correction. This may take a minute...")

# 6. Loop through each tree
for (i in 1:nrow(trees_sf)) {
  
  CPA_s <- trees_sf$Area_Num[i]
  H_s   <- trees_sf$Height_Num[i] 
  
  # --- DYNAMIC SEARCH RADIUS LOGIC ---
  current_spacing <- trees_sf$Spacing_Num[i]
  if (is.na(current_spacing) || current_spacing <= 0) {
    dynamic_radius <- 8.0 
  } else {
    dynamic_radius <- current_spacing * 1.5
  }
  
  # Find competitors
  comp_indices <- which(dist_matrix[i, ] <= dynamic_radius & dist_matrix[i, ] > 0)
  
  raw_ci <- 0 # Default to 0
  
  # --- CALCULATE RAW COMPETITION ---
  if (length(comp_indices) > 0) {
    CPA_c   <- trees_sf$Area_Num[comp_indices]
    H_c     <- trees_sf$Height_Num[comp_indices]
    Dist_sc <- dist_matrix[i, comp_indices]
    
    if (isTRUE(CPA_s > 0) && isTRUE(H_s > 0)) {
      valid_comps <- which(CPA_c > 0 & H_c > 0)
      if (length(valid_comps) > 0) {
        competition_terms <- (CPA_c[valid_comps] * H_c[valid_comps]) / (CPA_s * H_s * Dist_sc[valid_comps])
        raw_ci <- sum(competition_terms, na.rm = TRUE)
      }
    } else {
      raw_ci <- NA 
    }
  }
  
  raw_ci_values[i] <- raw_ci
  
  # --- AREA-EXPANSION EDGE CORRECTION ---
  if (!is.na(raw_ci) && raw_ci > 0) {
    
    # Draw the search circle around the subject tree
    tree_buffer <- st_buffer(trees_pts[i, ], dist = dynamic_radius)
    
    # Get the specific boundary for this tree's plot
    current_plot_name <- trees_sf$Plt_shp[i]
    plot_poly <- plot_boundaries %>% filter(Plt_shp == current_plot_name)
    
    # Calculate the overlap (suppress warnings about spatial attributes)
    suppressWarnings({
      overlap <- st_intersection(tree_buffer, plot_poly)
    })
    
    area_inside <- as.numeric(st_area(overlap))
    area_total <- as.numeric(st_area(tree_buffer))
    
    prop_inside <- area_inside / area_total
    
    # Safety catch: cap proportion to avoid breaking the math
    if (is.na(prop_inside) || prop_inside < 0.1) prop_inside <- 1.0
    if (prop_inside > 1.0) prop_inside <- 1.0 
    
    # Mathematically inflate the CI based on the missing slice of the circle
    adj_ci_values[i] <- raw_ci / prop_inside
    prop_inside_values[i] <- prop_inside
    
  } else {
    adj_ci_values[i] <- raw_ci
    prop_inside_values[i] <- 1.0
  }
}

# 7. Append the calculated spatial values back to our dataset
trees_sf$Castagneri_CI_Raw <- raw_ci_values
trees_sf$Castagneri_CI_Adj <- adj_ci_values
trees_sf$Edge_Overlap_Pct <- round(prop_inside_values * 100, 1)

# --- 8. DATA TRANSFORMATION & PRESENTATION ---
print("Applying statistical transformations...")

# We use the Edge-Adjusted CI for all these calculations
valid_ci <- trees_sf$Castagneri_CI_Adj

# Log Transformation (Adding 1 so we don't take the log of 0)
trees_sf$CI_Log <- log(valid_ci + 1)

# Z-Score Standardization (Mean = 0, Standard Deviation = 1)
trees_sf$CI_Zscore <- as.numeric(scale(valid_ci, center = TRUE, scale = TRUE))

# Categorical Stress Bands (Using Quartiles: 0-25%, 25-75%, 75-100%)
breaks <- quantile(valid_ci, probs = c(0, 0.25, 0.75, 1.0), na.rm = TRUE)
trees_sf$CI_Class <- cut(
  valid_ci, 
  breaks = breaks, 
  labels = c("Free Growing", "Moderate Competition", "Severe Stress"),
  include.lowest = TRUE
)

# 9. Export the results to a clean CSV
final_df <- st_drop_geometry(trees_sf)
write.csv(final_df, paste0("E:/Remote Sensing Media/",date_folder,"/10. Competition indices/All_Plots_with_Competition_Index.csv"), row.names = FALSE)

print("Pipeline Complete! Check your folder for the final CSV.")