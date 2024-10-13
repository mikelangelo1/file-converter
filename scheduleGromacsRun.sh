#!/bin/bash

# Check if GROMACS is installed
if ! command -v gmx &> /dev/null
then
    echo "GROMACS is not installed or not in PATH. Please install GROMACS and try again."
    exit 1
fi

# Check if at command is available
if ! command -v at &> /dev/null
then
    echo "The 'at' command is not available. Please install it and try again."
    exit 1
fi

# Function to validate time format
validate_time() {
    if ! date -d "$1" >/dev/null 2>&1
    then
        echo "Invalid time format. Please use HH:MM format (24-hour clock)."
        exit 1
    fi
}

# Get the time to run the simulation
read -p "Enter the time to run the simulation (HH:MM, 24-hour format): " run_time

# Validate the entered time
validate_time "$run_time"

# Set variables for GROMACS files
TOPO="topol.top"
CONF="conf.gro"
INDEX="index.ndx"
MDP_EM="em.mdp"
MDP_NVT="nvt.mdp"
MDP_NPT="npt.mdp"
MDP_MD="md.mdp"

# Function to check if a file exists
check_file() {
    if [ ! -f "$1" ]; then
        echo "Error: File $1 not found in the current directory."
        exit 1
    fi
}

# Check if all required files exist
check_file "$TOPO"
check_file "$CONF"
check_file "$INDEX"
check_file "$MDP_EM"
check_file "$MDP_NVT"
check_file "$MDP_NPT"
check_file "$MDP_MD"

# Create a temporary script to run the GROMACS simulation
cat << EOF > run_gromacs_temp.sh
#!/bin/bash
cd "$(pwd)"

# Energy Minimization
gmx grompp -f $MDP_EM -c $CONF -p $TOPO -o em.tpr
gmx mdrun -v -deffnm em

# NVT Equilibration
gmx grompp -f $MDP_NVT -c em.gro -r em.gro -p $TOPO -n $INDEX -o nvt.tpr
gmx mdrun -v -deffnm nvt

# NPT Equilibration
gmx grompp -f $MDP_NPT -c nvt.gro -r nvt.gro -t nvt.cpt -p $TOPO -n $INDEX -o npt.tpr
gmx mdrun -v -deffnm npt

# Production MD
gmx grompp -f $MDP_MD -c npt.gro -t npt.cpt -p $TOPO -n $INDEX -o md.tpr
gmx mdrun -v -deffnm md

echo "GROMACS simulation completed at $(date)"
EOF

# Make the temporary script executable
chmod +x run_gromacs_temp.sh

# Schedule the job
at "$run_time" < run_gromacs_temp.sh

echo "GROMACS simulation scheduled to run at $run_time"
echo "You can check scheduled jobs with the 'atq' command"
echo "To cancel the job, use 'atrm' followed by the job number from 'atq'"