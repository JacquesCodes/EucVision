#!/usr/bin/env python3
"""
Pix4D Quality Report Standalone Bulk Extractor (v3 - Optimized for Pix4D)
----------------------------------------------------------------------
Instructions:
1. Link to folder with all PDF reports.
2. In your terminal/command prompt, run:
   pip install pypdf
3. Execute the script:
   python local_pix4d_extractor_v3.py --dir . --out final_all_flights_v3.csv
"""

import os
import re
import csv
import argparse
from pypdf import PdfReader

# Bulletproof extraction patterns specifically calibrated for standard Pix4D Quality Reports
PATTERNS = {
    # Summary Section
    'camera_model': re.compile(r'Camera Model Name\(s\)\s+(.*?)(?=\s*(?:Average Ground Sampling|GSD|Average GSD|$))', re.IGNORECASE),
    # Note: pypdf often extracts the "²" in "km²" as a separate glyph on its own
    # line (e.g. "0.077 km\n2\n / 7.7326 ha"), so the superscript digit must be
    # allowed to appear *after* whitespace/newlines, not just immediately after "km".
    'area_covered': re.compile(r'Area Covered\s*([\d.]+)\s*km\s*[²2]?\s*/\s*([\d.]+)\s*ha', re.IGNORECASE),
    'gsd_cm': re.compile(r'Average Ground Sampling Distance \(GSD\)\s+([\d.]+)\s*cm', re.IGNORECASE),
    
    # Image counts (calibrated vs total)
    'calib_ratio': re.compile(r'(\d+)\s+out of\s+(\d+)\s+images calibrated\s*\(([\d.]+%?)\)', re.IGNORECASE),
    
    # Matching and Keypoint extraction
    'median_keypoints': re.compile(r'median of\s+(\d+)\s+keypoints per image', re.IGNORECASE),
    'median_matched': re.compile(r'median of\s+([\d.]+)\s+matches per calibrated image', re.IGNORECASE),
    'keypoint_table_median': re.compile(r'Median\s+(\d+)\s+(\d+)', re.IGNORECASE),
    
    # Precision Metrics
    'mean_reprojection': re.compile(r'Mean Reprojection Error\s*\[pixels\]\s*([\d.]+)', re.IGNORECASE),
    'gcp_rms_summary': re.compile(r'mean RMS error\s*=\s*([\d.]+)\s*m', re.IGNORECASE),
    'gcp_rms_table': re.compile(r'RMS Error\s*\[m\]\s*([\d.]+)\s*([\d.]+)\s*([\d.]+)', re.IGNORECASE),
    
    # 3D Point Cloud Densification and Density
    'densified_points': re.compile(r'Number of 3D Densified Points\s+([\d,]+)', re.IGNORECASE),
    'point_density': re.compile(r'Average Density\s*\(\s*per\s+m\s*[³3]\s*\)\s*([\d,.]+)', re.IGNORECASE),
    
    # Operational processing times
    'densification_time': re.compile(r'Time for Point Cloud Densification\s+([\w\d\s:]+?)(?=\s*(?:Time for|$))', re.IGNORECASE),
    'dsm_time': re.compile(r'Time for DSM Generation\s+([\w\d\s:]+?)(?=\s*(?:Time for|$))', re.IGNORECASE),
    'orthomosaic_time': re.compile(r'Time for Orthomosaic Generation\s+([\w\d\s:]+?)(?=\s*(?:Time for|$))', re.IGNORECASE),
}

def clean_date_from_filename(filename):
    cleaned = filename.replace("Lourensford_", "").replace("Bottom Sector_", "").replace("Top Sector_", "")
    cleaned = cleaned.replace("Bottom Section_", "").replace("Top Section_", "")
    cleaned = cleaned.replace("Bottom Cross_", "").replace("Top Cross_", "")
    cleaned = cleaned.replace("Bottom cross_", "").replace("Top cross_", "")
    cleaned = cleaned.replace("Bottom Section Cross_", "").replace("Top Section Cross_", "")
    cleaned = cleaned.replace("Bottom Sector Cross Hatch3_", "").replace("Top Sector Cross Hatch3__", "")
    cleaned = re.sub(r'_(?:0\.6|0\.6cm|0\.6_report|report.*)', '', cleaned, flags=re.IGNORECASE)
    cleaned = cleaned.replace("_report.pdf", "").replace("_report", "").replace(".pdf", "")
    cleaned = cleaned.replace("Cross Hatch", "").replace("CrossHatch", "").replace("Cross_Hatch", "").replace("Cross", "")
    cleaned = cleaned.replace("Combined", "").replace("section", "").replace("Section", "").replace("Sector", "")
    cleaned = cleaned.replace("__", "_").replace("_", " ").strip()
    
    m = re.search(r'(\d{1,2}\s+[A-Za-z]+(?:\s+\d{4})?)', cleaned)
    if m:
        return m.group(1).strip()
    m = re.search(r'([A-Za-z]+\s+\d{4})', cleaned)
    if m:
        return m.group(1).strip()
    return cleaned

def parse_report_pdf(pdf_path):
    filename = os.path.basename(pdf_path)
    compartment = 'Unknown'
    if 'top' in filename.lower():
        compartment = 'Top'
    elif 'bottom' in filename.lower():
        compartment = 'Bottom'
        
    flight_path = 'Normal'
    if 'cross hatch' in filename.lower() or 'crosshatch' in filename.lower() or 'cross_hatch' in filename.lower():
        flight_path = 'Cross Hatch'
    elif 'cross' in filename.lower():
        flight_path = 'Cross'
        
    date_val = clean_date_from_filename(filename)
    if date_val == "SU Lourensford 1 September 2025":
        date_val = "1 September 2025"
        compartment = "Unknown"
        
    metrics = {
        'Filename': filename, 
        'Date': date_val, 
        'Compartment': compartment, 
        'Flight path': flight_path,
        'Camera Model': None,
        'Area Covered': None,
        'GSD (cm)': None, 
        'Total Images': None, 
        'Calibrated Images': None, 
        '% images calibrated': None,
        'Median keypoints per image': None, 
        'Median matched keypoints per image': None,
        'Mean reprojection error (px)': None, 
        'Total 3D Densified points': None, 
        'Point density (pts/m³)': None, 
        'GCP RMS error (m)': None,
        'Densification time': None, 
        'DSM time': None, 
        'Orthomosaic time': None
    }
    
    try:
        reader = PdfReader(pdf_path)
        full_text = ""
        # Read the first 12 pages of the PDF to extract summary details
        for i in range(min(12, len(reader.pages))):
            page_text = reader.pages[i].extract_text()
            if page_text:
                full_text += "\n" + page_text
                
        # 1. Camera Model
        cam_m = PATTERNS['camera_model'].search(full_text)
        if cam_m:
            metrics['Camera Model'] = cam_m.group(1).strip()
            
        # 2. Area Covered (Standardizing as e.g., "0.031 km² / 3.05 ha")
        area_m = PATTERNS['area_covered'].search(full_text)
        if area_m:
            metrics['Area Covered'] = f"{area_m.group(1)} km² / {area_m.group(2)} ha"
            
        # 3. Ground Sampling Distance (GSD)
        gsd_m = PATTERNS['gsd_cm'].search(full_text)
        if gsd_m:
            metrics['GSD (cm)'] = float(gsd_m.group(1))
            
        # 4. Image Calibration Counts (Calibrated / Total)
        calib_m = PATTERNS['calib_ratio'].search(full_text)
        if calib_m:
            metrics['Calibrated Images'] = int(calib_m.group(1))
            metrics['Total Images'] = int(calib_m.group(2))
            metrics['% images calibrated'] = calib_m.group(3).strip()
            
        # 5. Median Keypoints & Matches
        kp_m = PATTERNS['median_keypoints'].search(full_text)
        if kp_m:
            metrics['Median keypoints per image'] = int(kp_m.group(1))
            
        match_m = PATTERNS['median_matched'].search(full_text)
        if match_m:
            metrics['Median matched keypoints per image'] = float(match_m.group(1))
            
        if metrics['Median keypoints per image'] is None or metrics['Median matched keypoints per image'] is None:
            tbl_m = PATTERNS['keypoint_table_median'].search(full_text)
            if tbl_m:
                if metrics['Median keypoints per image'] is None:
                    metrics['Median keypoints per image'] = int(tbl_m.group(1))
                if metrics['Median matched keypoints per image'] is None:
                    metrics['Median matched keypoints per image'] = float(tbl_m.group(2))
                    
        # 6. Mean Reprojection Error (pixels)
        reproj_m = PATTERNS['mean_reprojection'].search(full_text)
        if reproj_m:
            metrics['Mean reprojection error (px)'] = float(reproj_m.group(1))
            
        # 7. Total 3D Densified Points
        dense_m = PATTERNS['densified_points'].search(full_text)
        if dense_m:
            metrics['Total 3D Densified points'] = int(dense_m.group(1).replace(',', ''))
            
        # 8. Point Density (pts/m³)
        density_m = PATTERNS['point_density'].search(full_text)
        if density_m:
            metrics['Point density (pts/m³)'] = float(density_m.group(1).replace(',', ''))
            
        # 9. Georeferencing GCP RMS Error
        gcp_m = PATTERNS['gcp_rms_summary'].search(full_text)
        if gcp_m:
            metrics['GCP RMS error (m)'] = float(gcp_m.group(1))
        else:
            gcp_tbl_m = PATTERNS['gcp_rms_table'].search(full_text)
            if gcp_tbl_m:
                metrics['GCP RMS error (m)'] = float(gcp_tbl_m.group(3))
                
        # 10. Operational Processing Times (Fixed split/strip typo)
        dens_time_m = PATTERNS['densification_time'].search(full_text)
        if dens_time_m:
            metrics['Densification time'] = dens_time_m.group(1).strip()
            
        dsm_time_m = PATTERNS['dsm_time'].search(full_text)
        if dsm_time_m:
            metrics['DSM time'] = dsm_time_m.group(1).strip()
            
        ortho_time_m = PATTERNS['orthomosaic_time'].search(full_text)
        if ortho_time_m:
            metrics['Orthomosaic time'] = ortho_time_m.group(1).strip()
                
    except Exception as e:
        print(f"Error parsing {filename}: {e}")
        
    return metrics

def main():
    parser = argparse.ArgumentParser(description="Bulk Extract Pix4D mapper Quality Report metrics.")
    parser.add_argument("--dir", default=r"C:\Users\jakev\Stellenbosch University\JacquesV B.Sc. skripsie M.Sc. project - Documents\Processed Data\EucVision\11. Flight Report Summary\All Flight Reports", help="Directory containing PDF reports.")
    parser.add_argument("--out", default=r"C:\Users\jakev\Stellenbosch University\JacquesV B.Sc. skripsie M.Sc. project - Documents\Processed Data\EucVision\01. Data Analysis\07. Flight Reports Summary.csv", help="Output CSV name.")
    args = parser.parse_args()
    
    pdf_files = [f for f in os.listdir(args.dir) if f.lower().endswith('.pdf')]
    if not pdf_files:
        print(f"No PDFs found in directory '{args.dir}'")
        return
        
    print(f"Found {len(pdf_files)} PDF reports. Extracting all metrics...")
    results = [parse_report_pdf(os.path.join(args.dir, f)) for f in pdf_files]
    
    fields = [
        "Filename", "Date", "Compartment", "Flight path", "Camera Model", "Area Covered", "GSD (cm)", 
        "Total Images", "Calibrated Images", "% images calibrated",
        "Median keypoints per image", "Median matched keypoints per image", 
        "Mean reprojection error (px)", "Total 3D Densified points", 
        "Point density (pts/m³)", "GCP RMS error (m)",
        "Densification time", "DSM time", "Orthomosaic time"
    ]
    with open(args.out, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)
    print(f"\nSuccessfully compiled '{args.out}'!")

if __name__ == '__main__':
    main()