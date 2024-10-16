import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

csv_files = [
    '/Users/mac/Desktop/file-converter/mean_binding_affinity/HBAC1_combined_binding_affinity_with_average.csv',
    '/Users/mac/Desktop/file-converter/mean_binding_affinity/HBAC2_combined_binding_affinity_with_average.csv',
    '/Users/mac/Desktop/file-converter/mean_binding_affinity/HBAC3_combined_binding_affinity_with_average.csv',
    '/Users/mac/Desktop/file-converter/mean_binding_affinity/HBAC4_combined_binding_affinity_with_average.csv',
    '/Users/mac/Desktop/file-converter/mean_binding_affinity/HBAC5_combined_binding_affinity_with_average.csv',
    '/Users/mac/Desktop/file-converter/mean_binding_affinity/HBAC6_combined_binding_affinity_with_average.csv',
    '/Users/mac/Desktop/file-converter/mean_binding_affinity/HBAC7_combined_binding_affinity_with_average.csv',
    '/Users/mac/Desktop/file-converter/mean_binding_affinity/HBAC8_combined_binding_affinity_with_average.csv',
    '/Users/mac/Desktop/file-converter/mean_binding_affinity/HBAC9_combined_binding_affinity_with_average.csv',
    '/Users/mac/Desktop/file-converter/mean_binding_affinity/HBAC10_combined_binding_affinity_with_average.csv',
    '/Users/mac/Desktop/file-converter/mean_binding_affinity/HBAC11_combined_binding_affinity_with_average.csv'
]  

# Create an empty dataframe to store the average binding affinity for each file
average_affinity = pd.DataFrame()

for i, file in enumerate(csv_files):
    data = pd.read_csv(file)
    
    if 'average_binding_affinity' in data.columns:
        average_affinity[f'HBAC_{i+1}'] = data['average_binding_affinity']
    else:
        print(f"Column 'average_binding_affinity' not found in {file}")


average_affinity = average_affinity.T

average_affinity.to_csv('combined_average_binding_affinity_HBAC.csv')

# Plotting the heatmap using seaborn
plt.figure(figsize=(25, 10))
sns.heatmap(average_affinity, annot=True, cmap='coolwarm', cbar=True, annot_kws={"size": 6})
plt.title('Heatmap of Average Binding Affinity (HBAC 1 to 11) to 61 ligands')
plt.ylabel('HBAC')
plt.yticks(rotation=45)  
plt.xlabel('Ligands')  
plt.xticks(rotation=45)  

plt.show()