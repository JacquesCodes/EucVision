import os
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

print("Loading shapefiles...")
plots = gpd.read_file("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/02. Templates/EucVision LidR Boundaries/EucVision LidR Boundaries.shp")
old_trees = gpd.read_file(f"E:/Remote Sensing Media/{previous_trees_folder}/08. Crown shape file/All_Plots.shp")

# Set your target plot number here (just the number!)
target_plot_id = 17       

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
# 4. CONVERT GEOMETRY TO PIXELS (ORIGINAL BOXES)
# ==========================================================
# THE FIX: We are using the EXACT original bounding boxes to capture the red leaves
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
# 6. RUN PREDICTION
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
# 7. PER-TREE BRIGHTNESS TRIMMING & SAFE ZONES
# ==========================================================
print("Trimming ground shadows while protecting canopy self-shadows...")
with rasterio.open(temp_tif) as src_rgb:
    r = src_rgb.read(1).astype(float)
    g = src_rgb.read(2).astype(float)
    b = src_rgb.read(3).astype(float)

# Calculate Brightness to trim ground shadows
brightness = 0.299 * r + 0.587 * g + 0.114 * b
shadow_threshold = 40 
non_shadow_mask = brightness > shadow_threshold

# NEW ADDITION: Create a "Safe Zone" mask from your old manual polygons
import rasterio.features
old_trees_safe_zone = rasterio.features.rasterize(
    [(geom, 1) for geom in trees_in_plot.geometry],
    out_shape=(new_height, new_width),
    transform=new_transform,
    fill=0,
    dtype=np.uint8
)

with rasterio.open(temp_mask_tif) as src_mask:
    sam_mask = src_mask.read(1)
    mask_meta = src_mask.meta

# Create a blank canvas to store the independent trees
final_trimmed_mask = np.zeros_like(sam_mask)

# Loop through each tree individually
for tree_id in np.unique(sam_mask):
    if tree_id == 0: 
        continue # Skip the empty background
        
    single_tree_mask = (sam_mask == tree_id)
    
    # THE FIX: Keep the pixel if it is inside the old tree (Safe Zone) 
    # OR if it is bright new growth spilling outside the old tree!
    valid_pixels = (old_trees_safe_zone == 1) | non_shadow_mask
    
    # Apply the logic
    trimmed_pixels = np.where(single_tree_mask & valid_pixels, 255, 0).astype(np.uint8)
    
    # Fill the holes for JUST this tree
    contours, _ = cv2.findContours(trimmed_pixels, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if contours:
        # Draw the solid tree onto the final mask using its original ID
        cv2.drawContours(final_trimmed_mask, contours, -1, int(tree_id), thickness=cv2.FILLED)

with rasterio.open(temp_mask_tif, 'w', **mask_meta) as dest:
    dest.write(final_trimmed_mask, 1)

# ==========================================================
# 8. EXPORT & ORDER MATCHING
# ==========================================================
print("Converting pixel masks to vector polygons...")
sam.raster_to_vector(temp_mask_tif, output_gpkg)

if os.path.exists(output_gpkg):
    test_crowns = gpd.read_file(output_gpkg)
    test_crowns['Area_m2'] = test_crowns.geometry.area
    
    old_centroids = trees_in_plot.copy()
    old_centroids.geometry = old_centroids.geometry.centroid
    
    joined = gpd.sjoin(old_centroids, test_crowns, how='left', predicate='intersects')
    joined = joined[~joined.index.duplicated(keep='first')]
    
    matched_geoms = []
    matched_areas = []
    
    for idx in trees_in_plot.index:
        sam_idx = joined.loc[idx, 'index_right']
        if pd.notna(sam_idx):
            matched_geoms.append(test_crowns.loc[sam_idx].geometry)
            matched_areas.append(test_crowns.loc[sam_idx, 'Area_m2'])
        else:
            matched_geoms.append(trees_in_plot.loc[idx].geometry)
            matched_areas.append(trees_in_plot.loc[idx].geometry.area)
            
    final_crowns = trees_in_plot.copy()
    final_crowns.geometry = matched_geoms
    final_crowns['Area_m2'] = matched_areas
    
    final_crowns.to_file(output_shp)
    os.remove(output_gpkg)
    print(f"Success! Segmented {len(final_crowns)} crowns with perfectly matching order. Saved to: {output_shp}")
    
if os.path.exists(temp_tif):
    os.remove(temp_tif)
if os.path.exists(temp_mask_tif):
    os.remove(temp_mask_tif)