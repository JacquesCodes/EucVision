import os
import cv2
import numpy as np
import geopandas as gpd
import rasterio
import random
from rasterio.windows import Window
from shapely.geometry import Polygon, box
import glob

# --- CONFIGURATION ---
BASE_PATH = "E:/Remote Sensing Media"
OUTPUT_BASE = "E:/EucVision_YOLO"
TILE_SIZE = 640
OVERLAP = 150 # Increased overlap for better training coverage
CLASS_ID = 0 

TARGET_FOLDERS = ["03", "04", "06", "07", "08", "09", "10", "11", "12", "13"]

# Setup directories
for split in ["train", "val"]:
    os.makedirs(os.path.join(OUTPUT_BASE, split, "images"), exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_BASE, split, "labels"), exist_ok=True)

def get_plot_filter(tif_name, all_tifs_in_folder):
    """
    Determines which plots to extract based on the filename logic.
    """
    name_lower = tif_name.lower()
    
    # If there is only one TIF in the folder, take all plots
    if len(all_tifs_in_folder) == 1:
        return 1, 74
    
    # Logic for Top section
    if "top" in name_lower:
        return 1, 21
    
    # Logic for Bottom/Lower section
    if "bottom" in name_lower or "lower" in name_lower:
        return 22, 74
        
    # Fallback: if filename is ambiguous but there are multiple, 
    # we assume it's a full-site flight unless "top/bottom" exists in others.
    return 1, 74

def export_tile(img_tile, polygons, tile_name, split):
    img_tile = cv2.cvtColor(img_tile, cv2.COLOR_RGB2BGR)
    cv2.imwrite(os.path.join(OUTPUT_BASE, split, "images", f"{tile_name}.jpg"), img_tile)
    
    label_path = os.path.join(OUTPUT_BASE, split, "labels", f"{tile_name}.txt")
    with open(label_path, "w") as f:
        for poly in polygons:
            coords = np.array(poly.exterior.coords)
            norm_coords = [f"{x/TILE_SIZE:.6f} {y/TILE_SIZE:.6f}" for x, y in coords]
            f.write(f"{CLASS_ID} {' '.join(norm_coords)}\n")

# --- MAIN LOOP ---
for folder_name in os.listdir(BASE_PATH):
    prefix = folder_name.split('.')[0].strip()
    if prefix not in TARGET_FOLDERS:
        continue

    print(f"\nProcessing Folder: {folder_name}")
    ortho_dir = os.path.join(BASE_PATH, folder_name, "01. Orthomosaics")
    shape_dir = os.path.join(BASE_PATH, folder_name, "08. Crown shape file")
    
    tif_files = glob.glob(os.path.join(ortho_dir, "*.tif"))
    shp_files = glob.glob(os.path.join(shape_dir, "*.shp"))
    
    if not shp_files or not tif_files:
        continue
    
    master_gdf = gpd.read_file(shp_files[0])

    for tif_path in tif_files:
        tif_filename = os.path.basename(tif_path)
        min_p, max_p = get_plot_filter(tif_filename, tif_files)
        
        # Filter the shapefile for the specific plots needed for this TIF
        # This handles Top (1-21), Bottom (22-74), and Cross flights accordingly
        gdf_filtered = master_gdf[(master_gdf['Plot'] >= min_p) & (master_gdf['Plot'] <= max_p)]
        
        print(f" -> TIF: {tif_filename} | Plots: {min_p}-{max_p} ({len(gdf_filtered)} trees)")
        
        with rasterio.open(tif_path) as src:
            gdf = gdf_filtered.to_crs(src.crs)
            
            for y in range(0, src.height - TILE_SIZE, TILE_SIZE - OVERLAP):
                for x in range(0, src.width - TILE_SIZE, TILE_SIZE - OVERLAP):
                    split = "val" if random.random() < 0.10 else "train"
                    window = Window(x, y, TILE_SIZE, TILE_SIZE)
                    tile_transform = src.window_transform(window)
                    
                    # Read only if the window isn't mostly empty/nodata
                    img = src.read([1, 2, 3], window=window)
                    if np.mean(img) < 15: continue 
                    
                    img = np.moveaxis(img, 0, -1)
                    tile_box = box(*src.window_bounds(window))
                    intersecting_trees = gdf[gdf.intersects(tile_box)]
                    
                    tile_polys = []
                    for _, tree in intersecting_trees.iterrows():
                        clipped_poly = tree.geometry.intersection(tile_box)
                        if clipped_poly.is_empty: continue
                        
                        # Handle both Polygon and MultiPolygon
                        parts = [clipped_poly] if isinstance(clipped_poly, Polygon) else list(clipped_poly.geoms)
                        for p in parts:
                            if not isinstance(p, Polygon) or p.area < 0.01: continue
                            local_coords = []
                            for px, py in p.exterior.coords:
                                row, col = ~tile_transform * (px, py)
                                local_coords.append((col, row))
                            tile_polys.append(Polygon(local_coords))
                    
                    if tile_polys:
                        clean_name = tif_filename.replace(".tif", "").replace(" ", "_")
                        unique_name = f"{prefix}_{clean_name}_{x}_{y}"
                        export_tile(img, tile_polys, unique_name, split)

print("\nDone! Dataset created.")
