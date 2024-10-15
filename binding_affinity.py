import pandas as pd
import os

# List of CSV file paths
csv_files = ['/Users/mac/Desktop/file-converter/Round-1/HDAC10.csv',
             '/Users/mac/Desktop/file-converter/Round-2/HDAC10_2.csv',
             '/Users/mac/Desktop/file-converter/Round-3/selected/HDAC10_3.csv'
             ]  

# Create an empty dataframe to store binding affinity data from each file
binding_affinity_df = pd.DataFrame()

# Loop through each file
for i, file in enumerate(csv_files):
    # Read the CSV file
    data = pd.read_csv(file)
    
    # Check if the 'binding_affinity' column exists
    if 'Binding Affinity' in data.columns:
        # Extract the 'binding_affinity' column and rename it to indicate the file number
        binding_affinity_df[f'binding_affinity_file_{i+1}'] = data['Binding Affinity']
    else:
        print(f"Column 'binding_affinity' not found in {file}")

# Calculate the average of the three binding affinity columns
binding_affinity_df['average_binding_affinity'] = binding_affinity_df.mean(axis=1)

# Display the resulting dataframe
print(binding_affinity_df)

# Optionally, save the combined data to a new CSV file
binding_affinity_df.to_csv('HBAC10_combined_binding_affinity_with_average.csv', index=False)
