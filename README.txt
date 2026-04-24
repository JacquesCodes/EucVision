EucVision: Plot-Level UAV SfM Processing Pipeline
Author: Jacques Vermeulen

Project: EucXylo (https://eucxylo.sun.ac.za/)

Overview
The EucVision pipeline is a modular, automated workflow designed to process high-resolution RGB UAV imagery (Structure from Motion point clouds and orthomosaics) to extract individual tree metrics for juvenile Eucalyptus stands. The pipeline dynamically handles varying plot boundaries, normalizes terrain using temporal ensemble DTMs, extracts 3D canopy models, and compiles the data for statistical analysis.

The R scripts are numbered chronologically. Scripts 01 through 05 must be run in sequence as they build the foundational spatial data. Scripts 06 through 08 are independent utility and visualization tools that rely on the outputs of the core engine.

Phase 1: The Core Processing Engine
Run these sequentially to process raw data into statistical datasets.

01. EucVision Baseline DTM Generator.R
Purpose: Processes raw, baseline photogrammetry point clouds across the entire site. It applies Statistical Outlier Removal (SOR) to prune noise, utilizes Progressive TIN Densification (PTD) to classify ground points beneath the canopy, and applies morphological sinkhole filling to output a smoothed, continuous Digital Terrain Model (DTM).

Execution: Run individually on the dates with the best ground visibility.

02. EucVision Ensemble DTM Fuser.R
Purpose: Stacks the individual temporal DTMs created in Script 01. It takes the pixel-wise maximum across multiple flights to confidently overwrite SfM sinkholes, producing a single, highly accurate "Ultimate Baseline DTM" for the entire study site.

03. EucVision Crown Polygons Merger.R
Purpose: Acts as the bridge between manual QGIS extraction and the automated R pipeline. It reads numerical plot shapefiles, merges them into a continuous spatial dataframe, and binds them to the master CSV template.

Execution: Must be run after manually generating QGIS crown polygons for a new flight.

04. EucVision SfM Pipeline.R
Purpose: The heavy-lifting batch processor. It loops through the designated dataset folders, crops the point clouds to the plot boundaries, normalizes height against the Ultimate Baseline DTM, generates 2D Canopy Height Models (CHMs), and dynamically extracts individual tree heights.

Note: This script utilizes parallel processing (future). Ensure sufficient RAM/CPU cooling before executing on large batches.

05. EucVision Master Dataset Compiler.R
Purpose: The final data aggregator. It loops through all processed dataset folders, extracts the plot-level CSV metrics, performs a full join with the physical field measurements, applies dynamic outlier filtering (99th percentile + 5m), and outputs the chronological 01. Master Dataset.csv.

Next Step: This CSV output is the direct input for the Python machine learning and metric derivation pipeline.

Phase 2: Visualization & Utility Tools
These can be run independently at any time after the core processing is complete.

06. EucVision Plot Visualizer.R
Purpose: An interactive 3D QA/QC sandbox. Uses the rgl library to render specific cropped point clouds, ground classifications, and CHMs. It includes commented-out Local Maximum Filter (LMF) algorithms for testing Individual Tree Detection (ITD) parameters.

07. EucVision Publication Maps.R
Purpose: Generates 2D, publication-ready cartography. Uses ggplot2 and tidyterra to render top-down smoothed CHM maps complete with Viridis color scaling, north arrows, and spatial scale bars for thesis figures.

08. EucVision Timelapse Generator.R
Purpose: Creates temporal MP4 animations. It loops through the chronologically sorted orthomosaics, strictly crops them to a target extent, overlays the static spacing plot boundaries and the dynamic, species-color-coded crown polygons, and renders the result as a high-resolution video.