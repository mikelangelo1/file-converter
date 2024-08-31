import re
import argparse
import numpy as np
from typing import List, Tuple
import os

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert multiple XVG files to CSV format.")
    parser.add_argument("-xvg", nargs='+', help="Input XVG files", required=True)
    parser.add_argument("-out_dir", help="Output directory for CSV files", default=".")
    return parser.parse_args()

def extract_metadata(file_path: str) -> Tuple[List[str], List[str]]:
    x_axis_pattern = re.compile(r'xaxis\s+label\s+"(.+?)"')
    legend_pattern = re.compile(r'legend\s+"(.+?)"')
    
    with open(file_path, 'r') as file:
        content = file.read()
        x_axis_labels = x_axis_pattern.findall(content)
        legends = legend_pattern.findall(content)
    
    return x_axis_labels, legends

def create_header(x_axis_labels: List[str], legends: List[str]) -> str:
    all_labels = x_axis_labels + legends
    return ','.join(all_labels)

def read_xvg_data(file_path: str) -> np.ndarray:
    with open(file_path, 'r') as file:
        lines = (line for line in file if not line.startswith(('@', '#')))
        return np.loadtxt(lines)

def write_csv(file_path: str, header: str, data: np.ndarray) -> None:
    with open(file_path, 'w') as file:
        file.write(f"{header}\n")
        np.savetxt(file, data, fmt="%.3f", delimiter=',')

def process_xvg_file(xvg_path: str, out_dir: str) -> None:
    try:
        x_axis_labels, legends = extract_metadata(xvg_path)
        header = create_header(x_axis_labels, legends)
        data = read_xvg_data(xvg_path)
        
        base_name = os.path.splitext(os.path.basename(xvg_path))[0]
        csv_path = os.path.join(out_dir, f"{base_name}.csv")
        
        write_csv(csv_path, header, data)
        print(f"Converted {xvg_path} to {csv_path}")
        print(f"Columns: {header}")
    except Exception as e:
        print(f"Error processing {xvg_path}: {str(e)}")

def main() -> None:
    args = parse_arguments()
    
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    
    for xvg_file in args.xvg:
        process_xvg_file(xvg_file, args.out_dir)
    
    print("Conversion complete.")

if __name__ == "__main__":
    main()