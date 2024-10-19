import os

def read_xvg_file(filename):
    """
    Read an XVG file and extract the result column.
    Skips comment lines starting with @ or #.
    """
    results = []
    with open(filename, 'r') as f:
        for line in f:
            # Skip comment and header lines
            if line.startswith('@') or line.startswith('#'):
                continue
            # Split line into columns and take the last column
            try:
                columns = line.strip().split()
                results.append(float(columns[-1]))
            except (ValueError, IndexError):
                continue
    return results

def combine_xvg_results(input_files, output_file):
    """
    Combine result columns from multiple XVG files into a new XVG file.
    """
    # Read all result columns
    all_results = []
    for filename in input_files:
        results = read_xvg_file(filename)
        all_results.append(results)
    
    # Ensure all files have the same number of data points
    lengths = [len(results) for results in all_results]
    if len(set(lengths)) != 1:
        raise ValueError(f"Files have different numbers of data points: {lengths}")
    
    # Write combined results to output file
    with open(output_file, 'w') as f:
        # Write header
        f.write("# Combined results from files: " + ", ".join(input_files) + "\n")
        f.write("@ title \"Combined Results\"\n")
        f.write("@ xaxis label \"Index\"\n")
        f.write("@ yaxis label \"Values\"\n")
        
        # Write data points
        for i in range(lengths[0]):
            line = f"{i}"
            for results in all_results:
                line += f"\t{results[i]}"
            f.write(line + "\n")

# Example usage
if __name__ == "__main__":
    # List your XVG files here
    input_files = [
        "/Users/mac/Desktop/file-converter/RMSD/Cluster_1_rmsd_4.5.xvg",
        "/Users/mac/Desktop/file-converter/RMSD/Cluster_2_rmsd_4.5.xvg",
        "/Users/mac/Desktop/file-converter/RMSD/mutant_rmsd_4.5.xvg",
        "/Users/mac/Desktop/file-converter/RMSD/wild_rmsd_4.5.xvg"
    ]
    output_file = "combined_rmsd_results.xvg"
    
    try:
        combine_xvg_results(input_files, output_file)
        print(f"Successfully combined results into {output_file}")
    except Exception as e:
        print(f"Error: {str(e)}")