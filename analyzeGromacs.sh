#!/bin/bash

# Function to create a directory if it doesn't exist
create_dir() {
    if [ ! -d "$1" ]; then
        mkdir "$1"
    fi
}

# Function to run Grace and save the result as a PNG
save_grace_result() {
    local input_file="$1"
    local output_file="$2"
    gracebat "$input_file" -hdevice PNG -autoscale xy -printfile "$output_file" -fixed 3840 2160
}

# Recenter and rewrap coordinates
create_dir "Rec_and_Rewrap"
cd Rec_and_Rewrap || exit
gmx trjconv -s ../md.tpr -f ../md.xtc -o md_center.xtc -center -pbc mol -ur compact << EOF
1   # Select "Protein" for centering
0   # Select "System" for output
EOF
cd ..

# RMSD Calculations
create_dir "RMSD"
cd RMSD || exit
gmx rms -s ../md.tpr -f ../Rec_and_Rewrap/md_center.xtc -n ../index.ndx -o rmsd.xvg -tu ns << EOF
1   # Select "Protein" for reference
3   # Select "C-alpha" for comparison
EOF
save_grace_result "rmsd.xvg" "rmsd.png"
cd ..

# RMSF Calculations
create_dir "RMSF"
cd RMSF || exit
gmx rmsf -s ../md.tpr -f ../Rec_and_Rewrap/md_center.xtc -n ../index.ndx -o rmsf.xvg << EOF
3   # Select "C-alpha" for comparison

EOF
save_grace_result "rmsf.xvg" "rmsf.png"
cd ..

# Hydrogen Bonds
create_dir "Hbond"
cd Hbond || exit
gmx hbond -s ../md.tpr -f ../Rec_and_Rewrap/md_center.xtc -n ../index.ndx -num hb.xvg -tu ns << EOF
1   # Select "Protein" for RMSF analysis
EOF
save_grace_result "hb.xvg" "hb.png"
cd ..

# Gyration Radius
create_dir "Gyration"
cd Gyration || exit
gmx gyrate -s ../md.tpr -f ../Rec_and_Rewrap/md_center.xtc -n ../index.ndx -o gyrate1.xvg -tu ns << EOF
1   # Select "Protein" for Gyration Radius calculation
EOF
save_grace_result "gyrate1.xvg" "gyrate1.png"
cd ..

# SASA Calculations
create_dir "SASA"
cd SASA || exit
gmx sasa -s ../md.tpr -f ../Rec_and_Rewrap/md_center.xtc -n ../index.ndx -o sasa.xvg -tu ns << EOF
1   # Select the group for SASA analysis (adjust if needed)
EOF
save_grace_result "sasa.xvg" "sasa.png"
cd ..

echo "All analysis completed and results saved."
