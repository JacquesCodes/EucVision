import pandas as pd
import numpy as np

# 1. Load the files
# CRITICAL: Read the tree columns as strings right from the start to prevent zero-stripping
template_df = pd.read_csv('C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data analysis/00. Dataset template.csv', dtype={'Tree': str})
impact_df = pd.read_csv('C:/Users/jakev/Desktop/03. IMPACT Tree layout and Codes.csv', dtype={'TreeRelCoord': str})

# 2. Build the cross-reference dataframe from IMPACT
impact_df['Compartment_ref'] = impact_df['Comp'].str.replace('C', '').astype(float)
impact_df['Line_ref'] = impact_df['Line_L'].str.replace('L', '').astype(float)
impact_df['Plot'] = impact_df['Plot_P'].str.replace('P', '').astype(float)

impact_df['Spacing_ref'] = impact_df['Spacing'].str.replace('m', '', regex=False)
impact_df['Spacing_ref'] = pd.to_numeric(impact_df['Spacing_ref'], errors='coerce')

culture_mapping = {'Mono': 'Single', 'Mix': 'Mix'}
impact_df['Culture_ref'] = impact_df['Culture'].map(culture_mapping)

# Restoring the requested 5-group taxonomy 
species_mapping = {
    'EU': 'Urophylla',
    'EG': 'Grandis',
    'GC': 'Grandis clone',
    'CX': 'Cladocalyx',
    'CA': 'Cloeziana'
}
impact_df['Species_ref'] = impact_df['Var_V'].map(species_mapping)

# --- THE COORDINATE STRING FIX ---

# Rule for IMPACT: Pad with a leading zero if the column digit is single (e.g., "3.3" -> "3.03")
def fix_impact_coord(val):
    if pd.isnull(val): return val
    val_str = str(val).strip()
    if '.' in val_str:
        row, col = val_str.split('.', 1)
        if len(col) == 1:
            col = '0' + col  
        return f"{row}.{col}"
    return val_str

# Rule for Template: Pad with a trailing zero if Excel stripped it (e.g., "3.1" -> "3.10")
def fix_template_coord(val):
    if pd.isnull(val): return val
    val_str = str(val).strip()
    if '.' in val_str:
        row, col = val_str.split('.', 1)
        if len(col) == 1:
            col = col + '0'  
        return f"{row}.{col}"
    return val_str

# Apply the custom string formatting to create an exact merge key
impact_df['Tree_Merge_Key'] = impact_df['TreeRelCoord'].apply(fix_impact_coord)
template_df['Tree_Merge_Key'] = template_df['Tree'].apply(fix_template_coord)

# Isolate reference columns
ref_cols = ['Plot', 'Tree_Merge_Key', 'Compartment_ref', 'Line_ref', 'Culture_ref', 'Spacing_ref', 'Species_ref']
cross_ref_df = impact_df[ref_cols].drop_duplicates(subset=['Plot', 'Tree_Merge_Key'])

# Format Template Plot column for merging
template_df['Plot'] = pd.to_numeric(template_df['Plot'], errors='coerce')

# 3. Merge based on Plot and the new exact String Key
merged_df = template_df.merge(cross_ref_df, on=['Plot', 'Tree_Merge_Key'], how='left')

# Replace any blank empty strings with true NaNs so fillna() operates correctly
merged_df = merged_df.replace(r'^\s*$', np.nan, regex=True)

# 4. Fill missing values in original columns with reference data
for col in ['Compartment', 'Line', 'Culture', 'Spacing', 'Species']:
    merged_df[col] = merged_df[col].fillna(merged_df[f'{col}_ref'])

# Update the status column to the agreed 10-character shapefile limit
if 'Death_Date' in merged_df.columns:
    merged_df = merged_df.rename(columns={'Death_Date': 'TreeStatus'})
elif 'Status' in merged_df.columns:
    merged_df = merged_df.rename(columns={'Status': 'TreeStatus'})

# 5. Clean up and export
final_cols = ['Compartment', 'Line', 'Plot', 'Culture', 'Spacing', 'Species', 'Tree', 'TreeStatus']
final_template = merged_df[final_cols].copy()

# Ensure clean integers where possible without crashing on decimals
for col in ['Compartment', 'Line', 'Plot']:
    final_template[col] = pd.to_numeric(final_template[col], errors='coerce').astype('Int64')

final_template['Spacing'] = pd.to_numeric(final_template['Spacing'], errors='coerce')

final_template.to_csv('C:/Users/jakev/Desktop/00. Dataset template_UPDATED.csv', index=False)
print("SUCCESS: Template updated utilizing corrected 5-group taxonomy and split-string coordinate matching.")