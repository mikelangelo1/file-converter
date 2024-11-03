#!/usr/bin/env python3
import os
import subprocess
import argparse
import sys
from pathlib import Path

def run_command(cmd, stdin_text=None, check_error=True):
    """
    Run a shell command with proper error handling
    """
    if stdin_text is not None:
        process = subprocess.Popen(
            cmd,
            shell=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        stdout, stderr = process.communicate(stdin_text)
    else:
        process = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True
        )
        stdout, stderr = process.stdout, process.stderr
    
    if check_error and process.returncode != 0:
        print(f"Error running command: {cmd}")
        print(f"Error message: {stderr}")
        sys.exit(1)
    
    return stdout, stderr

def prepare_trajectory(xtc_file, tpr_file, output_xtc):
    """
    Prepare trajectory by:
    1. Centering on protein
    2. Making molecule whole
    3. Removing periodic boundary conditions
    """
    print("Preparing trajectory...")
    
    # Center on protein and make molecule whole
    cmd = f"gmx trjconv -f {xtc_file} -s {tpr_file} -pbc mol -center"
    stdout, _ = run_command(cmd, stdin_text="1\n1\n")  # Select protein for both centering and output
    
    # Remove periodic boundary conditions
    temp_xtc = "temp_nopbc.xtc"
    cmd = f"gmx trjconv -f {output_xtc} -s {tpr_file} -pbc nojump -o {temp_xtc}"
    stdout, _ = run_command(cmd, stdin_text="1\n")  # Select protein
    
    # Move temporary file to final output
    os.rename(temp_xtc, output_xtc)

def extract_frames(tpr_file, xtc_file, output_dir="pdb_frames", stride=10, center=True, 
                  remove_water=True, align=True, reference_group="Backbone"):
    """
    Extract PDB frames from Gromacs trajectory with enhanced processing
    
    Parameters:
    -----------
    tpr_file : str
        Path to TPR file
    xtc_file : str
        Path to XTC trajectory file
    output_dir : str
        Directory to store output PDB files
    stride : int
        Extract every nth frame
    center : bool
        Center the protein in the box
    remove_water : bool
        Remove water and ions
    align : bool
        Align frames to first frame
    reference_group : str
        Group to use for alignment ("Backbone" or "C-alpha")
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create processed trajectory directory
    proc_dir = output_dir / "processed"
    proc_dir.mkdir(exist_ok=True)
    
    # Process trajectory first
    processed_xtc = proc_dir / "processed.xtc"
    if center or align:
        prepare_trajectory(xtc_file, tpr_file, processed_xtc)
    else:
        processed_xtc = xtc_file
    
    # Create index file with protein and necessary groups
    ndx_file = proc_dir / "index.ndx"
    if not ndx_file.exists():
        cmd = f"gmx make_ndx -f {tpr_file} -o {ndx_file}"
        if reference_group == "C-alpha":
            stdin_text = "1\nq\n"  # Select protein and quit
        else:  # Backbone
            stdin_text = "1\n3\nq\n"  # Select protein, backbone, and quit
        run_command(cmd, stdin_text=stdin_text)
    
    # Get total number of frames
    cmd = f"gmx check -f {processed_xtc}"
    stdout, _ = run_command(cmd)
    
    # Extract frames
    print(f"\nExtracting frames (stride={stride})...")
    frame_num = 0
    
    while True:
        output_file = output_dir / f'frame_{frame_num}.pdb'
        
        # Build trjconv command
        cmd = f'gmx trjconv -f {processed_xtc} -s {tpr_file} -n {ndx_file} -o {output_file} -dump {frame_num}'
        
        if remove_water:
            cmd += ' -pbc mol'
        
        # Select appropriate groups based on settings
        if remove_water:
            stdin_text = "1\n"  # Select protein-only
        else:
            stdin_text = "0\n"  # Select system
            
        stdout, stderr = run_command(cmd, stdin_text=stdin_text, check_error=False)
        
        # Check if we've reached the end of trajectory
        if stderr and 'Frame out of trajectory' in stderr:
            break
            
        print(f"Extracted frame {frame_num}")
        frame_num += stride
    
    print(f"\nExtracted {frame_num//stride} frames to {output_dir}")
    print("\nFrames are ready for POVME analysis!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract and prepare PDB frames from Gromacs trajectory for POVME')
    parser.add_argument('--tpr', required=True, help='TPR file')
    parser.add_argument('--xtc', required=True, help='XTC trajectory file')
    parser.add_argument('--output', default='pdb_frames', help='Output directory')
    parser.add_argument('--stride', type=int, default=10, help='Extract every nth frame')
    parser.add_argument('--no-center', action='store_false', dest='center', 
                      help='Disable protein centering')
    parser.add_argument('--keep-water', action='store_false', dest='remove_water',
                      help='Keep water molecules and ions')
    parser.add_argument('--no-align', action='store_false', dest='align',
                      help='Disable frame alignment')
    parser.add_argument('--reference', choices=['Backbone', 'C-alpha'], default='Backbone',
                      help='Reference group for alignment')
    
    args = parser.parse_args()
    
    # Print settings
    print("\nSettings:")
    print(f"TPR file: {args.tpr}")
    print(f"Trajectory: {args.xtc}")
    print(f"Output directory: {args.output}")
    print(f"Frame stride: {args.stride}")
    print(f"Center protein: {args.center}")
    print(f"Remove water: {args.remove_water}")
    print(f"Align frames: {args.align}")
    print(f"Reference group: {args.reference}")
    print("")
    
    extract_frames(args.tpr, args.xtc, args.output, args.stride, 
                  args.center, args.remove_water, args.align, args.reference)