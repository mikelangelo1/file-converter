#!/bin/bash

# Define DPI value
dpi=600

# Loop over all SVG files in the current folder
for svg_file in *.svg; do
    # Extract the base name without extension
    base_name=$(basename "$svg_file" .svg)

    # Define the output file name
    png_file="${base_name}.png"

    cairosvg "$svg_file" -d "$dpi" -o "$png_file"

    echo "Converted $svg_file to $png_file with DPI $dpi"
done