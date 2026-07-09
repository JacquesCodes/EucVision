import pandas as pd
import numpy as np

# 1. Load Data
df = pd.read_csv('C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data analysis/01. Master Dataset.csv')
df['Date'] = pd.to_datetime(df['Date'], format='%d-%m-%Y')
df['Tree_Height'] = pd.to_numeric(df['Tree_Height'], errors='coerce')

# 2. Sort and Calculate Height Change
df = df.sort_values(by=['Compartment', 'Plot', 'Tree', 'Date'])
df['Height_Change'] = df.groupby(['Compartment', 'Plot', 'Tree'])['Tree_Height'].diff()

# 3. Define the Chaos Signature Function
def chaotic_signature(group):
    # Count how many trees shrunk or grew significantly
    shrink_count = (group['Height_Change'] < -0.3).sum()
    grow_count = (group['Height_Change'] > 0.6).sum()
    total = len(group['Height_Change'].dropna())
    
    if total == 0:
        return pd.Series({'Shrink_Pct': 0, 'Grow_Pct': 0, 'Std_Dev': 0})
        
    std_dev = group['Height_Change'].std()
    
    return pd.Series({
        'Shrink_Pct': shrink_count / total, 
        'Grow_Pct': grow_count / total, 
        'Std_Dev': std_dev
    })

# 4. Apply to every plot for every date
chaos_stats = df.groupby(['Compartment', 'Plot', 'Date']).apply(chaotic_signature).reset_index()

# 5. Filter for True Reversals 
# (High Standard Deviation AND a mix of significant shrinkage and growth)
true_reversals = chaos_stats[
    (chaos_stats['Shrink_Pct'] >= 0.10) &  # At least 10% of plot shrunk
    (chaos_stats['Grow_Pct'] >= 0.10) &    # At least 10% of plot grew
    (chaos_stats['Std_Dev'] >= 0.8)        # The variance is biologically impossible
].sort_values('Std_Dev', ascending=False)

# Format as percentages for readability
true_reversals['Shrink_Pct'] = (true_reversals['Shrink_Pct'] * 100).round(1).astype(str) + '%'
true_reversals['Grow_Pct'] = (true_reversals['Grow_Pct'] * 100).round(1).astype(str) + '%'
true_reversals['Std_Dev'] = true_reversals['Std_Dev'].round(3)

print("CONFIRMED GRID REVERSALS (Chaotic Mix of Growth and Shrinkage):")
print(true_reversals[['Compartment', 'Plot', 'Date', 'Std_Dev', 'Shrink_Pct', 'Grow_Pct']].to_string(index=False))