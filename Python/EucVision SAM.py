import os
import urllib.request
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.mask import mask
from rasterio.transform import Affine
import numpy as np
import cv2
from samgeo import SamGeo
import torch

# ==========================================================
# 1. SETUP & LOAD FILES
# ==========================================================
date_folder = "14. 06 February 2026"
previous_trees_folder = "13. 29 January 2026"
orthomosaic_path = f"E:/Remote Sensing Media/{date_folder}/01. Orthomosaics/Top Section_6 Feb 2026_transparent_mosaic_group1.tif"

checkpoint_path = 'C:/Users/jakev/Documents/sam_vit_b_01ec64.pth' 

if not os.path.exists(checkpoint_path):
    print("Downloading SAM Base model (~375MB) to your project folder...")
    url = "https://dl.fbaipublicfiles.com/segment_anything/sam_vit_b_01ec64.pth"
    urllib.request.urlretrieve(url, checkpoint_path)
    print("Download complete!")

print("Loading shapefiles...")
plots = gpd.read_file("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LidR Boundaries/EucVision LidR Boundaries.shp")
old_trees = gpd.read_file(f"E:/Remote Sensing Media/{previous_trees_folder}/08. Crown shape file/All_Plots.shp")

# Set your target plot number here (just the number!)
target_plot_id = 21         

# Automatically generate the correctly formatted string (e.g., 19 -> "Plot 19", 1 -> "Plot 01")
target_plot_name = f"Plot {target_plot_id:02d}"

# Filter your shapefiles
single_plot_boundary = plots[plots['id'] == target_plot_id] 
trees_in_plot = old_trees[old_trees['Plt_shp'] == target_plot_name]

print(f"Found {len(trees_in_plot)} trees in {target_plot_name}.")

# ==========================================================
# 2. INITIALIZE MODEL
# ==========================================================
print("Loading SAM Model onto the GTX 1060...")
device = 'cuda' if torch.cuda.is_available() else 'cpu'

sam = SamGeo(
    model_type="vit_b",
    checkpoint=checkpoint_path,
    device=device,
    sam_kwargs=None,
)

if not hasattr(sam, 'predictor'):
    from segment_anything import SamPredictor
    sam.predictor = SamPredictor(sam.sam)

# ==========================================================
# 3. CLIP & DOWNSAMPLE IMAGE (1cm)
# ==========================================================
plot_geom = single_plot_boundary.geometry.iloc[0]
target_resolution = 0.01  

print("Clipping and downsampling orthomosaic to 1cm/pixel...")
with rasterio.open(orthomosaic_path) as src:
    out_image, out_transform = mask(src, [plot_geom], crop=True)
    out_meta = src.meta
    
    current_res = src.res[0] 
    scale_factor = current_res / target_resolution
    
    new_height = int(out_image.shape[1] * scale_factor)
    new_width = int(out_image.shape[2] * scale_factor)
    
    resampled_image = np.empty((out_image.shape[0], new_height, new_width), dtype=out_image.dtype)
    for i in range(out_image.shape[0]):
        resampled_image[i] = cv2.resize(out_image[i], (new_width, new_height), interpolation=cv2.INTER_LINEAR)
        
    scale_x = out_image.shape[2] / new_width
    scale_y = out_image.shape[1] / new_height
    new_transform = out_transform * Affine.scale(scale_x, scale_y)
    
temp_tif = f"temp_{target_plot_name.replace(' ', '_')}.tif"
out_meta.update({
    "driver": "GTiff", 
    "height": new_height, 
    "width": new_width, 
    "transform": new_transform
})

with rasterio.open(temp_tif, "w", **out_meta) as dest:
    dest.write(resampled_image)

# ==========================================================
# 4. CONVERT GEOMETRY TO PIXELS
# ==========================================================
geo_boxes = trees_in_plot.bounds.values.tolist()
inv_transform = ~new_transform 

pixel_boxes = []
for minx, miny, maxx, maxy in geo_boxes:
    col_min, row_max = inv_transform * (minx, miny)
    col_max, row_min = inv_transform * (maxx, maxy)
    
    px_minx, px_maxx = min(col_min, col_max), max(col_min, col_max)
    px_miny, px_maxy = min(row_min, row_max), max(row_min, row_max)
    pixel_boxes.append([px_minx, px_miny, px_maxx, px_maxy])

# ==========================================================
# 5. SAFE MONKEY-PATCH
# ==========================================================
print("Generating image embeddings on the GPU...")
sam.set_image(temp_tif)

if not hasattr(sam.predictor, "_original_predict_torch"):
    sam.predictor._original_predict_torch = sam.predictor.predict_torch

    def patched_predict_torch(point_coords=None, point_labels=None, boxes=None, mask_input=None, multimask_output=True, return_logits=False):
        import torch 
        if boxes is not None and not isinstance(boxes, torch.Tensor):
            boxes = torch.tensor(boxes, dtype=torch.float32, device=sam.predictor.device)
            boxes = sam.predictor.transform.apply_boxes_torch(boxes, sam.predictor.original_size)
            multimask_output = False 
            
        return sam.predictor._original_predict_torch(
            point_coords=point_coords, 
            point_labels=point_labels, 
            boxes=boxes, 
            mask_input=mask_input, 
            multimask_output=multimask_output, 
            return_logits=return_logits
        )

    sam.predictor.predict_torch = patched_predict_torch

# ==========================================================
# 6. RUN PREDICTION & EXPORT WITH ORDER MATCHING
# ==========================================================
print("Running SAM prediction...")
temp_mask_tif = f"temp_mask_{target_plot_name.replace(' ', '_')}.tif"
output_gpkg = f"E:/Remote Sensing Media/{date_folder}/08. Crown shape file/{target_plot_name}_SAM_Test.gpkg"
output_shp = f"E:/Remote Sensing Media/{date_folder}/08. Crown shape file/{target_plot_name}_SAM_Test.shp"

sam.predict(
    boxes=pixel_boxes, 
    output=temp_mask_tif,
    crs=old_trees.crs,
    index=0  
)

# ==========================================================
# NEW ADDITION: Trim external shadows but KEEP canopies solid
# ==========================================================
print("Trimming shadows and filling canopy holes...")
with rasterio.open(temp_tif) as src_rgb:
    r = src_rgb.read(1).astype(float)
    g = src_rgb.read(2).astype(float)
    b = src_rgb.read(3).astype(float)

# Calculate Excess Green Index
exg = 2 * g - r - b
# If you run it and find that the script is still including some shadows, increase the threshold to 20 or 25
# If you find that the script is being too aggressive and cutting off the dark green edges of the actual leaves, lower the threshold to 5 or 10.
exg_threshold = 5 
vegetation_mask = exg > exg_threshold

# Load SAM's raw predicted pixels
with rasterio.open(temp_mask_tif) as src_mask:
    sam_mask = src_mask.read(1)
    mask_meta = src_mask.meta

# 1. Initial Trim: Keep only vegetation pixels
raw_trimmed = np.where(vegetation_mask, sam_mask, 0)

# 2. Convert to a binary image so OpenCV can read it (0 = background, 255 = tree)
binary_trimmed = (raw_trimmed > 0).astype(np.uint8) * 255

# 3. HOLE FILLER: Find ONLY the external boundaries of the trees
contours, _ = cv2.findContours(binary_trimmed, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

# Create a blank canvas and draw those outer boundaries completely solid
solid_binary = np.zeros_like(binary_trimmed)
cv2.drawContours(solid_binary, contours, -1, 255, thickness=cv2.FILLED)

# 4. Apply the solid, hole-free mask back to SAM's original data
final_trimmed_mask = np.where(solid_binary > 0, sam_mask, 0).astype(sam_mask.dtype)

# Overwrite SAM's mask with our clean, solid version
with rasterio.open(temp_mask_tif, 'w', **mask_meta) as dest:
    dest.write(final_trimmed_mask, 1)
# ==========================================================

print("Converting pixel masks to vector polygons...")
sam.raster_to_vector(temp_mask_tif, output_gpkg)

if os.path.exists(output_gpkg):
    test_crowns = gpd.read_file(output_gpkg)
    
    # 1. Calculate the area of the new SAM polygons (in square meters)
    test_crowns['Area_m2'] = test_crowns.geometry.area
    
    # 2. Extract the exact center points of your old trees
    old_centroids = trees_in_plot.copy()
    old_centroids.geometry = old_centroids.geometry.centroid
    
    # 3. Match the old centers to the new SAM polygons
    joined = gpd.sjoin(old_centroids, test_crowns, how='left', predicate='intersects')
    
    # Drop duplicates just in case a center point falls exactly on a boundary between two SAM crowns
    joined = joined[~joined.index.duplicated(keep='first')]
    
    matched_geoms = []
    matched_areas = []
    
    # 4. Loop through the exact original order
    for idx in trees_in_plot.index:
        sam_idx = joined.loc[idx, 'index_right']
        
        # If SAM successfully generated a polygon here, grab it
        if pd.notna(sam_idx):
            matched_geoms.append(test_crowns.loc[sam_idx].geometry)
            matched_areas.append(test_crowns.loc[sam_idx, 'Area_m2'])
        # FAILSAFE: If SAM failed to draw a tree, keep the original manual polygon so you don't lose the row
        else:
            matched_geoms.append(trees_in_plot.loc[idx].geometry)
            matched_areas.append(trees_in_plot.loc[idx].geometry.area)
            
    # 5. Overwrite the geometries and add the new Area column to your perfectly ordered original dataset
    final_crowns = trees_in_plot.copy()
    final_crowns.geometry = matched_geoms
    final_crowns['Area_m2'] = matched_areas
    
    # Save it
    final_crowns.to_file(output_shp)
    os.remove(output_gpkg)
    print(f"Success! Segmented {len(final_crowns)} crowns with perfectly matching order. Saved to: {output_shp}")
    
# Clean up temporary raster files
if os.path.exists(temp_tif):
    os.remove(temp_tif)
if os.path.exists(temp_mask_tif):
    os.remove(temp_mask_tif)