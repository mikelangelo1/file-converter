#!/bin/bash

# Set desired DPI, defaulting to 600 if not specified
DPI=${1:-600}

# Loop through all .ps files in the current directory
for psfile in *.ps; do
    outputfile="${psfile%.ps}.png"
    echo "Converting $psfile to $outputfile at $DPI DPI..."
    gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r$DPI -sOutputFile="$outputfile" "$psfile"
done

echo "Conversion completed!"
