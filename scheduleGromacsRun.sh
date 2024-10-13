#!/bin/bash

# Max retry attempts
MAX_RETRIES=3

# Check if GROMACS is installed
check_gromacs() {
    if ! command -v gmx &> /dev/null
    then
        echo "GROMACS is not installed or not in PATH. Please install GROMACS and try again."
        return 1
    fi
    return 0
}

# Function to validate time format
validate_time() {
    if ! date -d "$1" >/dev/null 2>&1
    then
        echo "Invalid time format. Please use HH:MM format (24-hour clock)."
        return 1
    fi
    return 0
}

# Function to check if a file exists
check_file() {
    if [ ! -f "$1" ]; then
        echo "Error: File $1 not found in the current directory."
        return 1
    fi
    return 0
}

# Function to prompt user for retry or exit
prompt_retry_or_exit() {
    while true; do
        read -p "Do you want to retry? (y/n): " choice
        case "$choice" in 
            y|Y ) return 0;;
            n|N ) echo "Exiting script."; exit 0;;
            * ) echo "Please answer y or n.";;
        esac
    done
}

# Set variables for GROMACS files
TOPO="topol.top"
CONF="conf.gro"
INDEX="index.ndx"
MDP_EM="em.mdp"
MDP_NVT="nvt.mdp"
MDP_NPT="npt.mdp"
MDP_MD="md.mdp"

# Function to run GROMACS steps
run_gromacs() {
    # Energy Minimization
    if ! gmx grompp -f $MDP_EM -c $CONF -p $TOPO -o em.tpr || \
       ! gmx mdrun -v -deffnm em; then
        echo "Error during Energy Minimization. Exiting GROMACS steps."
        return 1
    fi

    # NVT Equilibration
    if ! gmx grompp -f $MDP_NVT -c em.gro -r em.gro -p $TOPO -n $INDEX -o nvt.tpr || \
       ! gmx mdrun -v -deffnm nvt; then
        echo "Error during NVT Equilibration. Exiting GROMACS steps."
        return 1
    fi

    # NPT Equilibration
    if ! gmx grompp -f $MDP_NPT -c nvt.gro -r nvt.gro -t nvt.cpt -p $TOPO -n $INDEX -o npt.tpr || \
       ! gmx mdrun -v -deffnm npt; then
        echo "Error during NPT Equilibration. Exiting GROMACS steps."
        return 1
    fi

    # Production MD
    if ! gmx grompp -f $MDP_MD -c npt.gro -t npt.cpt -p $TOPO -n $INDEX -o md.tpr || \
       ! gmx mdrun -v -deffnm md; then
        echo "Error during Production MD. Exiting GROMACS steps."
        return 1
    fi

    return 0
}

retry_count=0

while true; do
    # Check GROMACS installation
    while ! check_gromacs; do
        if [ "$retry_count" -ge "$MAX_RETRIES" ]; then
            echo "Maximum retries reached. Exiting script."
            exit 1
        fi
        retry_count=$((retry_count+1))
        echo "Retry attempt $retry_count of $MAX_RETRIES."
        if [ "$retry_count" -lt "$MAX_RETRIES" ]; then
            prompt_retry_or_exit
        else
            echo "Maximum retries reached. Exiting script."
            exit 1
        fi
    done

    # Reset retry count after successful check
    retry_count=0

    # Get the time to run the simulation
    while true; do
        read -p "Enter the time to run the simulation (HH:MM, 24-hour format): " run_time
        if validate_time "$run_time"; then
            break
        fi
        if [ "$retry_count" -ge "$MAX_RETRIES" ]; then
            echo "Maximum retries reached. Exiting script."
            exit 1
        fi
        retry_count=$((retry_count+1))
        echo "Retry attempt $retry_count of $MAX_RETRIES."
        if [ "$retry_count" -lt "$MAX_RETRIES" ]; then
            prompt_retry_or_exit
        else
            echo "Maximum retries reached. Exiting script."
            exit 1
        fi
    done

    # Reset retry count after successful time validation
    retry_count=0

    # Check if all required files exist
    files_missing=false
    for file in "$TOPO" "$CONF" "$INDEX" "$MDP_EM" "$MDP_NVT" "$MDP_NPT" "$MDP_MD"; do
        if ! check_file "$file"; then
            files_missing=true
        fi
    done

    if [ "$files_missing" = true ]; then
        if [ "$retry_count" -ge "$MAX_RETRIES" ]; then
            echo "Maximum retries reached. Exiting script."
            exit 1
        fi
        retry_count=$((retry_count+1))
        echo "Retry attempt $retry_count of $MAX_RETRIES."
        if [ "$retry_count" -lt "$MAX_RETRIES" ]; then
            prompt_retry_or_exit
        else
            echo "Maximum retries reached. Exiting script."
            exit 1
        fi
        continue
    fi

    # Reset retry count after successful file checks
    retry_count=0

    echo "All required files found. Waiting to start simulation at $run_time..."

    # Wait until the specified time
    current_time=$(date +%s)
    target_time=$(date -d "$run_time" +%s)

    if [ $target_time -lt $current_time ]; then
        target_time=$(date -d "tomorrow $run_time" +%s)
    fi

    sleep_duration=$((target_time - current_time))
    sleep $sleep_duration

    echo "Starting GROMACS simulation at $(date)"

    # Run GROMACS simulation
    if run_gromacs; then
        echo "GROMACS simulation completed successfully at $(date)"
        break  # Exit loop on successful completion
    else
        if [ "$retry_count" -ge "$MAX_RETRIES" ]; then
            echo "Maximum retries reached during GROMACS execution. Exiting script."
            exit 1
        fi
        retry_count=$((retry_count+1))
        echo "Retry attempt $retry_count of $MAX_RETRIES."
        if [ "$retry_count" -lt "$MAX_RETRIES" ]; then
            prompt_retry_or_exit
        else
            echo "Maximum retries reached. Exiting script."
            exit 1
        fi
    fi
done
