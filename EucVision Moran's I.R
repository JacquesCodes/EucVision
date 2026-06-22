library(sf)
library(spdep)
library(dplyr)
library(tibble)

# 1. Attach residuals from all 6 models to their respective dataframes
# Height models
df_h$resid_sp <- residuals(models_h$species)
df_h$resid_sc <- residuals(models_h$spacing)

# Crown models
df_c$resid_sp <- residuals(models_c$species)
df_c$resid_sc <- residuals(models_c$spacing)

# CA:H Ratio models
df_r$resid_sp <- residuals(models_r$species)
df_r$resid_sc <- residuals(models_r$spacing)

# 2. Define your dates and fixed file paths
date_start <- "2025-09-01"
shp_start  <- "E:/Remote Sensing Media/02. 01 September 2025/09. Crown Metrics/Crown_Metrics_01_September_2025.shp"

date_end   <- "2026-05-25"
shp_end    <- "E:/Remote Sensing Media/27. 25 May 2026/09. Crown Metrics/Crown_Metrics_25_May_2026.shp"


# 3. Create a custom function to run the Moran's I test
run_moran <- function(data, target_date, shp_path, resid_col, response_name, model_type) {
  
  # Load the shapefile quietly
  crowns <- st_read(shp_path, quiet = TRUE) %>%
    mutate(Tree = round(as.numeric(Tree), 2))
  
  # Filter the dataset for the target date
  df_sub <- data %>%
    filter(Date == as.Date(target_date)) %>%
    mutate(Tree = round(as.numeric(Tree), 2))
  
  # Grab the specific residual column we are testing
  df_sub$target_resid <- df_sub[[resid_col]]
  
  # Join spatial data and drop any trees not included in the model
  joined <- crowns %>%
    left_join(df_sub %>% select(Compartment, Line, Plot, Tree, target_resid), 
              by = c("Cmprtmn" = "Compartment", "Line", "Plot", "Tree")) %>%
    filter(!is.na(target_resid))
  
  # Build the spatial weights matrix (6 closest neighbors)
  coords <- st_coordinates(st_centroid(joined))
  nb <- knn2nb(knearneigh(coords, k = 6))
  lw <- nb2listw(nb, style = "W")
  
  # Run the test
  m_test <- moran.test(joined$target_resid, lw)
  
  # Return a tidy row of results
  tibble(
    Date = target_date,
    Response = response_name,
    Model = model_type,
    Moran_I_Stat = round(m_test$estimate[1], 4),
    p_value = format.pval(m_test$p.value, eps = 0.001)
  )
}

# 4. Run the function for all 12 combinations and combine into one table
cat("Calculating Moran's I for all models across both dates. This may take a few seconds...\n")

moran_results <- bind_rows(
  # --- September 1, 2025 (Start) ---
  run_moran(df_h, date_start, shp_start, "resid_sp", "Height", "Species"),
  run_moran(df_h, date_start, shp_start, "resid_sc", "Height", "Spacing"),
  run_moran(df_c, date_start, shp_start, "resid_sp", "Crown Area", "Species"),
  run_moran(df_c, date_start, shp_start, "resid_sc", "Crown Area", "Spacing"),
  run_moran(df_r, date_start, shp_start, "resid_sp", "CA:H Ratio", "Species"),
  run_moran(df_r, date_start, shp_start, "resid_sc", "CA:H Ratio", "Spacing"),
  
  # --- May 25, 2026 (End) ---
  run_moran(df_h, date_end, shp_end, "resid_sp", "Height", "Species"),
  run_moran(df_h, date_end, shp_end, "resid_sc", "Height", "Spacing"),
  run_moran(df_c, date_end, shp_end, "resid_sp", "Crown Area", "Species"),
  run_moran(df_c, date_end, shp_end, "resid_sc", "Crown Area", "Spacing"),
  run_moran(df_r, date_end, shp_end, "resid_sp", "CA:H Ratio", "Species"),
  run_moran(df_r, date_end, shp_end, "resid_sc", "CA:H Ratio", "Spacing")
)

# 5. Print a clean table to the console
cat("\n================ MORAN'S I RESULTS ================\n")
print(as.data.frame(moran_results), row.names = FALSE)
cat("===================================================\n")