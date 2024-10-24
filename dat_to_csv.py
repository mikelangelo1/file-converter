import csv

def dat_to_csv(dat_file, csv_file, delimiter=' '):
    """
    Converts a .dat file to a .csv file.
    
    :param dat_file: Path to the .dat file
    :param csv_file: Path to the output .csv file
    :param delimiter: Delimiter used in the .dat file (default is space ' ')
    """
    with open(dat_file, 'r') as infile, open(csv_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter=delimiter)
        writer = csv.writer(outfile)
        
        for row in reader:
            writer.writerow(row)
    print(f'Converted {dat_file} to {csv_file}')

# usage
dat_file_path = '/Users/mac/Desktop/file-converter/FINAL_RESULTS_MMPBSA.dat'  
csv_file_path = 'output.csv'   # Output CSV file path
dat_to_csv(dat_file_path, csv_file_path, delimiter=' ')  # Specify delimiter if needed
