import re
import argparse
import numpy as np
from typing import List, Tuple

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert XVG file to CSV format.")
    parser.add_argument("-xvg", help="Input XVG file", default="test.xvg")
    parser.add_argument("-csv", help="Output CSV file", default="test.csv")
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

def main() -> None:
    args = parse_arguments()
    
    x_axis_labels, legends = extract_metadata(args.xvg)
    header = create_header(x_axis_labels, legends)
    data = read_xvg_data(args.xvg)
    
    write_csv(args.csv, header, data)
    print(f"Conversion complete. Output written to {args.csv}")
    print(f"Columns: {header}")

if __name__ == "__main__":
    main()