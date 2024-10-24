import sys

def order_rmsf_xvg(input_file, output_file):
    """
    Read RMSF XVG file and renumber residues sequentially from 1.
    Maintains original data order, only changes the numbering.
    
    Args:
        input_file (str): Path to input XVG file
        output_file (str): Path to output ordered XVG file
    """
    # Store header lines and data separately
    header_lines = []
    data_points = []
    
    # Read input file
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                header_lines.append(line)
            else:
                # Split line into residue number and RMSF value
                values = line.strip().split()
                if len(values) >= 2:  # Ensure line has at least 2 values
                    try:
                        rmsf = float(values[1])
                        data_points.append(rmsf)
                    except ValueError:
                        continue

    # Write output file
    with open(output_file, 'w') as f:
        # Write header lines
        for line in header_lines:
            f.write(line)
        
        # Write data with sequential numbering
        for i, rmsf in enumerate(data_points, 1):
            f.write(f"{i:8.3f} {rmsf:10.5f}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.xvg output.xvg")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    order_rmsf_xvg(input_file, output_file)