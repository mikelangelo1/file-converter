# file-converter
# XVG to CSV Converter

This Python script converts multiple XVG (GROMACS visualization) files to CSV (Comma-Separated Values) format. It's designed to handle batch processing of XVG files, making it useful for researchers and data analysts working with molecular dynamics simulation data.

## Features

- Convert multiple XVG files to CSV in a single run
- Extract metadata (x-axis labels and legends) from XVG files
- Specify custom output directory for CSV files
- Error handling for individual file processing
- Informative output messages

## Requirements

- Python 3.6+
- NumPy library

## Installation

1. Ensure you have Python 3.6 or higher installed on your system.
2. Install the required NumPy library:

   ```
   pip install numpy
   ```

3. Download the `xvg_to_csv_converter.py` script to your local machine.

## Usage

Run the script from the command line with the following syntax:

```
python xvg_to_csv_converter.py -xvg [XVG_FILES] -out_dir [OUTPUT_DIRECTORY]
```

### Arguments:

- `-xvg`: One or more XVG files to convert (required)
- `-out_dir`: Directory to save the output CSV files (optional, default is current directory)

### Example:

```
python xvg_to_csv_converter.py -xvg file1.xvg file2.xvg file3.xvg -out_dir ./output_csvs
```

This command will convert `file1.xvg`, `file2.xvg`, and `file3.xvg` to CSV format and save the resulting files in the `./output_csvs` directory.

## Output

For each successfully processed XVG file, the script will:

1. Create a corresponding CSV file in the specified output directory
2. Print a message confirming the conversion
3. Display the column headers of the created CSV file

If any errors occur during processing, the script will print an error message for the affected file and continue with the next file.

## Note

This script assumes that the XVG files follow the standard GROMACS format. If your XVG files have a different structure, you may need to modify the `extract_metadata` function in the script.

## Contributing

Feel free to fork this repository and submit pull requests with any enhancements. You can also open issues for any bugs found or features you'd like to see added.

## License

This script is released under the MIT License. See the LICENSE file for more details.