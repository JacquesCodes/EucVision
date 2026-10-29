import laspy
import geopandas as gpd
import rasterio

# Define file paths using raw strings (r"") to handle Windows backslashes
las_path = r"E:\Remote Sensing Media\50. 21 September 2026 (ALS)\03. Point Clouds\21 September 2026 ALS Point cloud.las"
shp_path = r"E:\Remote Sensing Media\50. 21 September 2026 (ALS)\08. Crown Polygons\Crown_Polygons_21_September_2026_(ALS).shp"
dtm_path = r"E:\Remote Sensing Media\00. Baseline DTM\DTM_from_UAV_EPSG_2048.tif"

print("--- Point Cloud (.las) ---")
with laspy.open(las_path) as las:
    try:
        las_crs = las.header.parse_crs()
        print(f"CRS: {las_crs}")
    except Exception:
        print("CRS: Not defined in header VLRs (very common for raw ALS data)")
    
    print(f"Extent (minX, minY, maxX, maxY): {las.header.x_min:.3f}, {las.header.y_min:.3f}, {las.header.x_max:.3f}, {las.header.y_max:.3f}\n")

print("--- Crown Polygons (.shp) ---")
gdf = gpd.read_file(shp_path)
print(f"CRS: {gdf.crs}")
bounds = gdf.total_bounds
print(f"Extent (minX, minY, maxX, maxY): {bounds[0]:.3f}, {bounds[1]:.3f}, {bounds[2]:.3f}, {bounds[3]:.3f}\n")

print("--- Baseline DTM (.tif) ---")
with rasterio.open(dtm_path) as src:
    print(f"CRS: {src.crs}")
    print(f"Extent (minX, minY, maxX, maxY): {src.bounds.left:.3f}, {src.bounds.bottom:.3f}, {src.bounds.right:.3f}, {src.bounds.top:.3f}")
    

import laspy
import numpy as np

las_path = r"E:\Remote Sensing Media\50. 21 September 2026 (ALS)\03. Point Clouds\21 September 2026 ALS Point cloud.las"
out_path = r"E:\Remote Sensing Media\50. 21 September 2026 (ALS)\03. Point Clouds\21 September 2026 ALS Point cloud_Negative.las"

print("Reading point cloud...")
las = laspy.read(las_path)

print("Inverting X and Y coordinates...")
# 1. Extract coordinates as independent numpy arrays
new_x = np.array(las.x) * -1.0
new_y = np.array(las.y) * -1.0

print("Updating header offsets...")
# 2. Update the header offsets FIRST to accommodate the new negative space.
# Setting the offset to the minimum value prevents integer overflow.
las.header.offsets = [np.min(new_x), np.min(new_y), las.header.offsets[2]]

print("Applying new coordinates...")
# 3. Assign the new coordinates back to the las object
las.x = new_x
las.y = new_y

print("Saving corrected point cloud...")
las.write(out_path)
print("Done! Use the new _Negative.las file for your normalisation.")