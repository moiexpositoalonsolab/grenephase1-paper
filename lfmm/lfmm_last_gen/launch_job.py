#!/usr/bin/env python3

"""
LFMM Analysis Job Launcher (Last Generation)

This script launches LFMM analysis jobs for multiple bioclimatic variables using the optimal K value
determined by variance_explained_k.r. It creates and submits SLURM jobs for each bioclimatic variable
for the last generation data.

Required Input Files:
1. lfmm_fullresults_all_k/variance_explained_k.csv
   - Results from variance_explained_k.r containing optimal K value
2. env_train_files/env_train_*.csv
   - Environmental training data files for each bioclimatic variable

Output:
- Creates and submits SLURM job scripts for each bioclimatic variable
- Each job runs run_lfmm_full_genome_ridge_multiplebiovars.r with the appropriate parameters
"""

import os
import pandas as pd
import subprocess

# --- Configuration ---
# Base directories
BASE_DIR = '.'
RESULTS_DIR = os.path.join(BASE_DIR, 'lfmm_fullresults_all_k')
ENV_TRAIN_DIR = os.path.join(BASE_DIR, 'env_train_files')

# List of bioclimatic variables to analyze
BIOVARS = ['bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10',
           'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19']

def get_optimal_k():
    """Read the variance explained results and determine optimal K value."""
    var_explained_file = os.path.join(RESULTS_DIR, 'variance_explained_k.csv')
    df = pd.read_csv(var_explained_file)
    
    # Find K value where variance explained starts to level off
    # This is a simple heuristic - you might want to adjust this logic
    optimal_k = df.loc[df['Variance_Explained'].diff().abs() < 0.01, 'K'].iloc[0]
    return optimal_k

def create_job_script(biovar, k):
    """Create a SLURM job script for a specific bioclimatic variable."""
    script_content = f"""#!/bin/bash
#SBATCH --job-name=lfmm_lastgen_{biovar}
#SBATCH --output=lfmm_lastgen_{biovar}_%j.out
#SBATCH --error=lfmm_lastgen_{biovar}_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

# Load required modules
module load R/4.1.0

# Run the LFMM analysis
Rscript run_lfmm_full_genome_ridge_multiplebiovars.r "{biovar}" {k}
"""
    
    # Write the script to a file
    script_path = os.path.join(BASE_DIR, f'run_lfmm_lastgen_{biovar}.sh')
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    return script_path

def main():
    """Main function to launch jobs for all bioclimatic variables."""
    # Get optimal K value
    optimal_k = get_optimal_k()
    print(f"Using optimal K value: {optimal_k}")
    
    # Create and submit jobs for each bioclimatic variable
    for biovar in BIOVARS:
        print(f"Processing {biovar}...")
        
        # Create job script
        script_path = create_job_script(biovar, optimal_k)
        
        # Submit the job
        try:
            subprocess.run(['sbatch', script_path], check=True)
            print(f"Successfully submitted job for {biovar}")
        except subprocess.CalledProcessError as e:
            print(f"Error submitting job for {biovar}: {e}")

if __name__ == "__main__":
    main()
