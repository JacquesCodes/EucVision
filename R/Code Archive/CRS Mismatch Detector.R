# ──────────────────────────────────────────────────────────────────────────────
# EUCVISION: SPATIAL AUDIT & CRS MISMATCH DETECTOR ####
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Email: Jacques.Stellies@gmail.com
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────
# Description: Recursively scans a target directory for all .shp and .tif files.
#              Audits the assigned Coordinate Reference System (CRS) against the 
#              actual geometric extent (X-minimum) to detect South African 
#              axis-flipping mismatches. Generates a diagnostic summary report.
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup and Imports ####
# ──────────────────────────────────────────────────────────────────────────────
library(terra)
library(sf)

# ──────────────────────────────────────────────────────────────────────────────
# 2. Configuration ####
# ──────────────────────────────────────────────────────────────────────────────
# base_dir <- "E:/Remote Sensing Media"

base_dir <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/04. QGIS Combined Output"

print("Scanning directories for .shp and .tif files... This may take a moment.")

# Find all shapefiles and tiffs recursively
all_files <- list.files(base_dir, pattern = "\\.(shp|tif)$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)

# Filter out script-generated temporary folders to avoid double-counting
files_to_audit <- all_files[!grepl("terra_temp|TimeLapse_Outputs", all_files, ignore.case = TRUE)]

print(paste("Found", length(files_to_audit), "files to audit. Beginning inspection..."))

# Initialize an empty list to store the results
audit_results <- list()

# ──────────────────────────────────────────────────────────────────────────────
# 3. Spatial Audit Loop ####
# ──────────────────────────────────────────────────────────────────────────────
for (i in seq_along(files_to_audit)) {
  
  file_path <- files_to_audit[i]
  file_ext <- tolower(tools::file_ext(file_path))
  file_name <- basename(file_path)
  
  # Default fallback values
  crs_label <- "Unknown/Missing"
  extent_sign <- "Unknown"
  status <- "Error/Unreadable"
  
  # Process TIFF files
  if (file_ext == "tif") {
    tryCatch({
      r <- rast(file_path)
      crs_wkt <- crs(r)
      xmin <- ext(r)[1]
      
      # Determine assigned CRS label
      if (grepl("2048", crs_wkt)) crs_label <- "EPSG:2048"
      else if (grepl("102562", crs_wkt)) crs_label <- "ESRI:102562"
      else if (crs_wkt != "") crs_label <- "Other"
      
      # Determine true geometry sign
      extent_sign <- ifelse(xmin < 0, "Negative", "Positive")
      
    }, error = function(e) {
      status <<- "Corrupted/Unreadable TIFF"
    })
  }
  
  # Process SHP files
  if (file_ext == "shp") {
    tryCatch({
      v <- st_read(file_path, quiet = TRUE)
      
      # Ensure geometry exists before checking bounding box
      if (nrow(v) > 0) {
        crs_wkt <- st_crs(v)$wkt
        xmin <- st_bbox(v)[1]
        
        # Check if CRS is missing completely
        if (is.na(crs_wkt)) {
          crs_label <- "Missing .prj file"
        } else {
          # Determine assigned CRS label
          if (grepl("2048", crs_wkt)) crs_label <- "EPSG:2048"
          else if (grepl("102562", crs_wkt)) crs_label <- "ESRI:102562"
          else crs_label <- "Other"
        }
        
        # Determine true geometry sign
        extent_sign <- ifelse(xmin < 0, "Negative", "Positive")
        
      } else {
        crs_label <- "Empty Shapefile"
        extent_sign <- "None"
      }
      
    }, error = function(e) {
      status <<- "Corrupted/Unreadable SHP"
    })
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # 4. Diagnostic Classification Logic ####
  # ────────────────────────────────────────────────────────────────────────────
  # Only classify if the file was successfully read
  if (status != "Corrupted/Unreadable TIFF" && status != "Corrupted/Unreadable SHP" && crs_label != "Empty Shapefile") {
    
    if (crs_label == "EPSG:2048" && extent_sign == "Negative") {
      status <- "Valid EPSG:2048"
      
    } else if (crs_label == "EPSG:2048" && extent_sign == "Positive") {
      status <- "FAKE EPSG:2048 (Disguised ESRI:102562)"
      
    } else if (crs_label == "ESRI:102562" && extent_sign == "Positive") {
      status <- "Valid ESRI:102562"
      
    } else if (crs_label == "ESRI:102562" && extent_sign == "Negative") {
      status <- "FAKE ESRI:102562 (Disguised EPSG:2048)"
      
    } else if (crs_label == "Missing .prj file") {
      status <- paste("Missing CRS - Extent is", extent_sign)
      
    } else {
      status <- paste("Other CRS - Extent is", extent_sign)
    }
  }
  
  # Save the result for this file
  audit_results[[i]] <- data.frame(
    File_Name = file_name,
    Type = toupper(file_ext),
    Assigned_CRS = crs_label,
    Geometry_Extent = extent_sign,
    Diagnostic_Status = status,
    Path = file_path,
    stringsAsFactors = FALSE
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# 5. Compile and Report Results ####
# ──────────────────────────────────────────────────────────────────────────────
# Combine all results into a single data frame
df_results <- do.call(rbind, audit_results)

# Print a clean summary table to the console
cat("\n\n──────────────────────────────────────────────────────────────────────────────\n")
cat("                     SPATIAL CRS AUDIT SUMMARY\n")
cat("──────────────────────────────────────────────────────────────────────────────\n\n")

# Create a cross-tabulation of File Types vs Diagnostic Status
summary_table <- table(df_results$Diagnostic_Status, df_results$Type)
print(summary_table)

cat("\n──────────────────────────────────────────────────────────────────────────────\n")
cat("Note: To view the specific files causing issues, you can filter the 'df_results'\n")
cat("dataframe. For example:\n")
cat("View(df_results[df_results$Diagnostic_Status == 'FAKE EPSG:2048 (Disguised ESRI:102562)', ])\n")