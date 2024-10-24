import os
import cairosvg

def convert_all_svg_to_png(folder_path, dpi=300):
    """
    Convert all SVG files in the specified folder to PNG with a given DPI.

    :param folder_path: Path to the folder containing SVG files.
    :param dpi: DPI (dots per inch) resolution for the output PNG (default: 300).
    """
    # Convert DPI to scale factor (default DPI for SVG is 96)
    scale_factor = dpi / 96.0

    # Iterate over all files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith(".svg"):
            svg_file_path = os.path.join(folder_path, filename)
            png_file_name = f"{os.path.splitext(filename)[0]}.png"
            png_file_path = os.path.join(folder_path, png_file_name)
            
            # Convert SVG to PNG
            with open(svg_file_path, "rb") as svg_file:
                svg_data = svg_file.read()
                cairosvg.svg2png(bytestring=svg_data, write_to=png_file_path, scale=scale_factor)
                
            print(f"Converted {filename} to {png_file_name} with DPI {dpi}")

# Example usage:
folder_path = "/Users/mac/Desktop/file-converter/7.4-decomps/wild-decomp" 
desired_dpi = 600

convert_all_svg_to_png(folder_path, dpi=desired_dpi)
