#!/bin/bash

# Check if GROMACS is installed
if ! command -v gmx &> /dev/null
then
    echo "GROMACS is not installed or not in PATH. Please install GROMACS and try again."
    exit 1
fi

# Set variables
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

echo "All required files are present. Proceeding with the simulation."

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

echo "GROMACS simulation completed successfully!"