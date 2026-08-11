import os
import cv2
import numpy as np
import geopandas as gpd
import pandas as pd
import rasterio
from rasterio.windows import Window
from shapely.geometry import Polygon
from ultralytics import YOLO

# ==========================================================
# 1. SETUP & PATHS
# ==========================================================
# YOLO Weights
model_path = "C:/Users/jakev/Downloads/runs_segment_Eucalyptus_L40S_Final_weights_best.pt"

# The new dataset you want to process
ortho_path = "E:/Remote Sensing Media/14. 06 February 2026/01. Orthomosaics/Top Section_6 Feb 2026_transparent_mosaic_group1.tif"

# Your original manual master shapefile to pull IDs from
manual_master_path = "E:/Remote Sensing Media/13. 29 January 2026/08. Crown shape file/All_Plots.shp"

# The final, tracked output
output_shp = "E:/Remote Sensing Media/14. 06 February 2026/08. Crown shape file/Top_Section_Tracked.shp"

# Plot filter for this specific flight (e.g., Top Section = Plots 1-21)
MIN_PLOT = 1
MAX_PLOT = 21

TILE_SIZE = 640
OVERLAP = 100

os.makedirs(os.path.dirname(output_shp), exist_ok=True)

# ==========================================================
# 2. LOAD MODEL & DATA
# ==========================================================
print("Loading Custom YOLOv8 Model...")
model = YOLO(model_path)

print("Loading Manual Master Shapefile...")
manual_trees = gpd.read_file(manual_master_path)

# Filter for the target plots
target_trees = manual_trees[(manual_trees['Plot'] >= MIN_PLOT) & (manual_trees['Plot'] <= MAX_PLOT)].copy()
print(f"Tracking {len(target_trees)} specific trees for Plots {MIN_PLOT}-{MAX_PLOT}...")

# ==========================================================
# 3. SLIDING WINDOW INFERENCE
# ==========================================================
all_polygons = []
all_confidences = []

print("Scanning orthomosaic and predicting canopies...")
with rasterio.open(ortho_path) as src:
    crs = src.crs
    
    for y in range(0, src.height, TILE_SIZE - OVERLAP):
        for x in range(0, src.width, TILE_SIZE - OVERLAP):
            window = Window(x, y, TILE_SIZE, TILE_SIZE)
            tile_transform = src.window_transform(window)
            
            img = src.read([1, 2, 3], window=window)
            if np.mean(img) < 15: continue  
            
            img = np.moveaxis(img, 0, -1)
            img_bgr = cv2.cvtColor(img, cv2.COLOR_RGB2BGR)
            
            results = model.predict(img_bgr, imgsz=TILE_SIZE, conf=0.25, verbose=False)
            result = results[0]
            
            if result.masks is not None:
                for idx, mask_coords in enumerate(result.masks.xy):
                    if len(mask_coords) < 3: continue
                    
                    spatial_coords = []
                    for px, py in mask_coords:
                        global_x, global_y = tile_transform * (px, py)
                        spatial_coords.append((global_x, global_y))
                        
                    poly = Polygon(spatial_coords)
                    
                    if poly.is_valid and poly.area > 0.5: 
                        all_polygons.append(poly)
                        all_confidences.append(float(result.boxes.conf[idx]))

# ==========================================================
# 4. CLEANUP OVERLAPS (NON-MAXIMUM SUPPRESSION)
# ==========================================================
print(f"Raw detection found {len(all_polygons)} potential crowns. Cleaning overlaps...")
yolo_raw_gdf = gpd.GeoDataFrame({'confidence': all_confidences}, geometry=all_polygons, crs=crs)
yolo_raw_gdf = yolo_raw_gdf.sort_values('confidence', ascending=False).reset_index(drop=True)

centroids = yolo_raw_gdf.geometry.centroid
spatial_join = gpd.sjoin(gpd.GeoDataFrame(geometry=centroids, crs=crs), yolo_raw_gdf, how='inner', predicate='within')

# FIX 1: Drop duplicate indices to prevent a single YOLO polygon from being copied twice
clean_indices = spatial_join[~spatial_join.index.duplicated(keep='first')]['index_right']
clean_indices = clean_indices.drop_duplicates()

yolo_clean_gdf = yolo_raw_gdf.iloc[clean_indices].copy()
yolo_clean_gdf['YOLO_Area'] = yolo_clean_gdf.geometry.area

# FIX 2: Reset the index so .loc[] lookups always return a single row
yolo_clean_gdf = yolo_clean_gdf.reset_index(drop=True)
print(f"Cleaned YOLO output contains {len(yolo_clean_gdf)} unique crowns.")

# ==========================================================
# 5. THE SPATIAL JOIN (TIME-SERIES TRACKING)
# ==========================================================
print("Matching YOLO predictions to original tree IDs...")

# Ensure CRS matches before joining
if yolo_clean_gdf.crs != target_trees.crs:
    yolo_clean_gdf = yolo_clean_gdf.to_crs(target_trees.crs)

# FIX 3: Reset the master shapefile index to ensure completely clean lookups
target_trees = target_trees.reset_index(drop=True)

old_centroids = target_trees.copy()
old_centroids.geometry = old_centroids.geometry.centroid

# Match old centers to new YOLO polygons
joined = gpd.sjoin(old_centroids, yolo_clean_gdf, how='left', predicate='intersects')
joined = joined[~joined.index.duplicated(keep='first')]

matched_geoms = []
matched_areas = []
ai_detected_flag = []

for idx in target_trees.index:
    yolo_idx = joined.loc[idx, 'index_right']
    
    # SUCCESS: YOLO found the tree
    if pd.notna(yolo_idx):
        yolo_idx = int(yolo_idx) # Convert the float back to an integer index
        matched_geoms.append(yolo_clean_gdf.loc[yolo_idx].geometry)
        matched_areas.append(yolo_clean_gdf.loc[yolo_idx, 'YOLO_Area'])
        ai_detected_flag.append("Yes")
        
    # FAILSAFE: YOLO missed the tree, fallback to manual
    else:
        matched_geoms.append(target_trees.loc[idx].geometry)
        matched_areas.append(target_trees.loc[idx].geometry.area)
        ai_detected_flag.append("No - Manual Fallback")

# ==========================================================
# 6. EXPORT FINAL DATASET
# ==========================================================
final_tracked = target_trees.copy()
final_tracked.geometry = matched_geoms
final_tracked['Area_m2'] = matched_areas
final_tracked['YOLO_Found'] = ai_detected_flag 

final_tracked.to_file(output_shp)
print("\n==========================================================")
print("TASK COMPLETE!")
print(f"Total Trees Tracked: {len(final_tracked)}")
print(f"Successfully updated by YOLO: {final_tracked['YOLO_Found'].value_counts().get('Yes', 0)}")
print(f"Saved to: {output_shp}")
print("==========================================================")