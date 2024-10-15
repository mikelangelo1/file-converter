import csv
from collections import defaultdict

def find_most_negative_binding_affinity_row(csv_file_path):
    ligand_data = {}
    
    with open(csv_file_path, 'r') as file:
        reader = csv.reader(file)
        
        # Get the header
        header = next(reader)
        
        # Process each row in the CSV
        for row in reader:
            ligand = row[0]
            binding_affinity = float(row[1])
            rmsd_ub = float(row[2])
            rmsd_lb = float(row[3])
            
            # If this ligand hasn't been seen before, or if this binding affinity is more negative
            if ligand not in ligand_data or binding_affinity < ligand_data[ligand]['binding_affinity']:
                ligand_data[ligand] = {
                    'binding_affinity': binding_affinity,
                    'rmsd_ub': rmsd_ub,
                    'rmsd_lb': rmsd_lb,
                    'full_row': row
                }
    
    return header, ligand_data

def main(input_csv_path, output_csv_path):
    # Find the row with the most negative binding affinity for each ligand
    header, result = find_most_negative_binding_affinity_row(input_csv_path)

    # Print the results
    print("Header:", ", ".join(header))
    for ligand, data in result.items():
        print(f"Ligand: {ligand}")
        print(f"  Binding Affinity: {data['binding_affinity']}")
        print(f"  RMSD/UB: {data['rmsd_ub']}")
        print(f"  RMSD/LB: {data['rmsd_lb']}")
        print()

    # Write the results to a new CSV file
    with open(output_csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for data in result.values():
            writer.writerow(data['full_row'])

    print(f"Results have been written to '{output_csv_path}'")

if __name__ == "__main__":
    input_file = "/Users/mac/Desktop/hackbio-internship/stage-4/hdac11/hdac11-3.csv" 
    output_file = "HDAC11_3.csv"
    main(input_file, output_file)