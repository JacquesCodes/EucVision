import pandas as pd

# 1. Setup file paths
field_path = r"C:\Users\jakev\Stellenbosch University\JacquesV B.Sc. skripsie M.Sc. project - Documents\Processed Data\EucVision\01. Data analysis\02. Field Measurements.csv"
drone_path = r"C:\Users\jakev\Downloads\02. Processed Master Dataset.csv"

print("Loading datasets...")
field_df = pd.read_csv(field_path, low_memory=False)
proc_df = pd.read_csv(drone_path, low_memory=False)

# 2. Identify the 'Suspects' from the March 2026 Field Data
march_field = field_df[field_df['Date'] == '23-03-2026'].copy()

# Generate unique Tree_IDs
march_field['Tree_ID'] = (march_field['Compartment'].astype(str) + '_' +
                          march_field['Line'].astype(str) + '_' +
                          march_field['Plot'].astype(str) + '_' +
                          pd.to_numeric(march_field['Tree']).map('{:.2f}'.format))

# Trees that have a Stem_Diameter are definitively ALIVE
alive_confirmed = march_field.dropna(subset=['Stem_Diameter'])['Tree_ID'].unique()

# Trees that are in the baseline drone dataset but missing from the alive list are our SUSPECTS
all_baseline_trees = proc_df['Tree_ID'].unique()
suspects = [t for t in all_baseline_trees if t not in alive_confirmed]

print(f"Total Baseline Trees: {len(all_baseline_trees)}")
print(f"Confirmed Alive in Field: {len(alive_confirmed)}")
print(f"Unmeasured Suspects to Investigate: {len(suspects)}\n")

# 3. Analyze Temporal Drone Data for the Suspects
proc_df['Date'] = pd.to_datetime(proc_df['Date'], errors='coerce')
proc_df = proc_df.sort_values(by=['Tree_ID', 'Date'])

# Filter drone data to only look at our unmeasured suspects
suspect_data = proc_df[proc_df['Tree_ID'].isin(suspects)].copy()

# Calculate week-to-week absolute drops (m2) and percentage drops (%)
suspect_data['Area_Drop'] = suspect_data.groupby('Tree_ID')['Crown_Area'].diff()
suspect_data['Area_Pct_Change'] = suspect_data.groupby('Tree_ID')['Crown_Area'].pct_change()

# 4. Classify Suspects based on Canopy Collapse
hit_list = []
for tree_id in suspects:
    tree_history = suspect_data[suspect_data['Tree_ID'] == tree_id].dropna(subset=['Area_Drop'])
    
    if not tree_history.empty:
        # Find the row where the Area Drop was the most severe (minimum percentage change)
        worst_row = tree_history.loc[tree_history['Area_Pct_Change'].idxmin()]
        
        pct_drop = worst_row['Area_Pct_Change']
        abs_drop = worst_row['Area_Drop']
        
        # THRESHOLD LOGIC: If a tree loses more than 15% of its canopy, OR more than 0.15m2 in a week, it collapsed.
        if pct_drop < -0.15 or abs_drop < -0.15:
            status = "Likely Dead (Catastrophic Collapse)"
        else:
            status = "Likely Alive (Field Measurement Discarded)"
            
        hit_list.append({
            'Tree_ID': tree_id,
            'Plot': int(worst_row['Plot']),
            'Tree': worst_row['Tree'],
            'Species': worst_row['Species'],
            'Status_Prediction': status,
            'Max_Area_Loss (m2)': round(abs_drop, 3),
            'Max_Area_Loss (%)': f"{round(pct_drop * 100, 1)}%",
            'Date_To_Check_QGIS': worst_row['Date'].strftime('%d-%m-%Y')
        })

# 5. Format and Export the Results
results_df = pd.DataFrame(hit_list)

# Sort so the "Likely Dead" ones are at the top, grouped by Date for easy QGIS checking
results_df = results_df.sort_values(by=['Status_Prediction', 'Date_To_Check_QGIS', 'Plot', 'Tree'])

output_path = r"C:\Users\jakev\Desktop\Intelligent_Mortality_List.csv"
results_df.to_csv(output_path, index=False)

print(f"SUCCESS! Intelligent list generated and saved to: {output_path}")

# Print a quick summary of what the algorithm found
summary = results_df['Status_Prediction'].value_counts()
print("\n--- ALGORITHM SUMMARY ---")
for status, count in summary.items():
    print(f"{status}: {count} trees")