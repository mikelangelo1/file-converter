#!/bin/bash

# Check if GROMACS is installed
if ! command -v gmx &> /dev/null
then
    echo "GROMACS is not installed or not in PATH. Please install GROMACS and try again."
    exit 1
fi

# Set variables
TRAJ="md.xtc"
TPR="md.tpr"
EDR="md.edr"
OUTPUT_DIR="analysis_results"

# Function to check if a file exists
check_file() {
    if [ ! -f "$1" ]; then
        echo "Error: File $1 not found in the current directory."
        exit 1
    fi
}

# Check if all required files exist
check_file "$TRAJ"
check_file "$TPR"
check_file "$EDR"

# Create output directory
mkdir -p $OUTPUT_DIR

# 1. Energy analysis
echo "0" | gmx energy -f $EDR -o $OUTPUT_DIR/potential.xvg
echo "Potential energy analysis completed."

# 2. RMSD analysis
echo "4 4" | gmx rms -s $TPR -f $TRAJ -o $OUTPUT_DIR/rmsd.xvg -tu ns
echo "RMSD analysis completed."

# 3. RMSF analysis
echo "4" | gmx rmsf -s $TPR -f $TRAJ -o $OUTPUT_DIR/rmsf.xvg -res
echo "RMSF analysis completed."

# 4. Radius of gyration analysis
echo "1" | gmx gyrate -s $TPR -f $TRAJ -o $OUTPUT_DIR/gyrate.xvg
echo "Radius of gyration analysis completed."

# 5. Hydrogen bond analysis
echo "1 1" | gmx hbond -s $TPR -f $TRAJ -num $OUTPUT_DIR/hbnum.xvg
echo "Hydrogen bond analysis completed."

# 6. Distance analysis (assuming you want to measure distance between groups 1 and 2)
echo "1 2" | gmx distance -s $TPR -f $TRAJ -select 'com of group 1 plus com of group 2' -oall $OUTPUT_DIR/distance.xvg
echo "Distance analysis completed."

# 7. Generate a PDB file for visualization
echo "1" | gmx trjconv -s $TPR -f $TRAJ -o $OUTPUT_DIR/trajectory.pdb -pbc mol -center
echo "PDB file generated for visualization."

echo "All analyses completed. Results are in the '$OUTPUT_DIR' directory."
echo "You can visualize the .xvg files using Grace or another plotting tool."
echo "The trajectory.pdb file can be visualized using VMD, PyMOL, or similar software."