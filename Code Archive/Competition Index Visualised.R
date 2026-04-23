# 1. Load required libraries
library(sf)
library(dplyr)
library(ggplot2)
library(viridis)

# 2. Define your folder path (matching your previous script setup)
date_folder <- "17. 02 March 2026"

# 3. Read the Original Shapefile AND the newly calculated CSV
trees_sf <- st_read(paste0("E:/Remote Sensing Media/", date_folder, "/09. Tree heights/All Plots.shp"))
ci_data <- read.csv(paste0("E:/Remote Sensing Media/", date_folder, "/10. Competition indices/All_Plots_with_Competition_Index.csv"))

# 4. Merge the Competition Index math back into the spatial geometries
# We use Plt_shp (Plot Name) and Tree (Tree Number) to link the data perfectly
trees_spatial <- trees_sf %>%
  left_join(ci_data %>% select(Plt_shp, Tree, Castagneri_CI_Adj), by = c("Plt_shp", "Tree"))

# 5. Filter for ONLY Plot 18
plot18_sf <- trees_spatial %>% filter(Plt_shp == "Plot 18")

# 6. Build the Visual Map
p_plot18 <- ggplot(data = plot18_sf) +
  
  # Map Color to Competition Index, and Size to Crown Area
  geom_sf(aes(color = Castagneri_CI_Adj, size = as.numeric(Area)), alpha = 0.85) +
  
  # Use the 'magma' palette. 'direction = -1' makes High Stress = Dark Purple/Black, Low Stress = Bright Yellow
  scale_color_viridis_c(option = "magma", name = "Competition\nStress (CI)", direction = -1) +
  
  # Scale the bubble sizes so crowns look proportional
  scale_size_continuous(range = c(2, 12), name = "Crown Area (m²)") +
  
  labs(
    title = "Spatial Competition Map: Plot 18",
    subtitle = "Visualizing the 'Crown Crush' effect: Size = Canopy, Color = Neighborhood Stress",
    x = "Easting",
    y = "Northing"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray40"),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 7. Display the map in your R console
print(p_plot18)

# 8. Save the map as a high-res image directly to your folder
ggsave(paste0("E:/Remote Sensing Media/", date_folder, "/10. Competition indices/Plot18_Competition_Map.png"), 
       plot = p_plot18, width = 10, height = 8, dpi = 300)

print("Plot 18 Spatial Map successfully generated and saved!")