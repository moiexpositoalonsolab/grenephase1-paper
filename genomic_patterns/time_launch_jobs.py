#!/usr/bin/env python
# coding: utf-8

# This script is a job launcher designed to automate the execution of the
# `calc_delta_time_allplots_binreg.py` script for multiple sites on a high-performance computing cluster (SLURM).
# It identifies site-plot combinations with sufficient data and generates and submits sbatch scripts
# to perform binomial regression of allele frequency change for each site.

# --- Required Input Files ---
# Ensure this file is present relative to the script's location:
# - ../key_files/merged_hapFIRE_allele_frequency.txt (Used to identify sites/plots with sufficient data)
# The script will generate and submit sbatch files to a SLURM cluster.
# The `calc_delta_time_allplots_binreg.py` script (also in this directory) is required for the launched jobs.

# --- Import Libraries ---

import pandas as pd
import random
import subprocess
import os # Import os for path joining

# --- Identify Sites/Plots with Sufficient Data ---

# Read the header of the allele frequency file to identify sample names and extract site, generation, and plot.
merged_hapFIRE_allele_frequency_header = pd.read_csv('../key_files/merged_hapFIRE_allele_frequency.txt', nrows=1, sep = '\t')
merged_hapFIRE_allele_frequency_header = merged_hapFIRE_allele_frequency_header.T.reset_index()

# Extract site, generation, and plot information from sample names
merged_hapFIRE_allele_frequency_header['site'] = merged_hapFIRE_allele_frequency_header['index'].str.split('_').str[0]
merged_hapFIRE_allele_frequency_header['gen'] = merged_hapFIRE_allele_frequency_header['index'].str.split('_').str[1]
merged_hapFIRE_allele_frequency_header['plot'] = merged_hapFIRE_allele_frequency_header['index'].str.split('_').str[2]

# Group by site and plot and count the number of samples (representing generations)
# Identify site-plot combinations that have data for at least 2 generations
samples_at_least_2_years = merged_hapFIRE_allele_frequency_header.groupby(['site', 'plot']).size().reset_index(name='n_samples')
sites_plots_w_at_least_2_years = samples_at_least_2_years[samples_at_least_2_years['n_samples'] > 1 ]

# Get a dictionary of sites with their corresponding plots that have at least 2 years of data
# The original script seemed to iterate over sites based on this, so we extract unique sites.
unique_sites_with_sufficient_data = sites_plots_w_at_least_2_years['site'].unique()

print(f"Identified {len(unique_sites_with_sufficient_data)} sites with at least 2 years of data for job submission.")

# --- Generate and Submit SLURM sbatch Jobs ---

# Create a directory for sbatch scripts and job outputs if it doesn't exist
sbatch_dir = 'sbatch_scripts'
os.makedirs(sbatch_dir, exist_ok=True)

shfiles = [] # List to store the names of the generated sbatch files

# Iterate through each unique site identified
for site in unique_sites_with_sufficient_data:
    # Generate a random seed (though not used in the sbatch script in the original code)
    # seed = random.randint(1,100000000)
    
    # Define the sbatch filename for the current site
    file = os.path.join(sbatch_dir, f'changep_time_{site}.sh')
    
    # Define the command to be executed on the cluster
    # This command runs the calc_delta_time_allplots_binreg.py script with the current site as argument.
    cmd = f'python calc_delta_time_allplots_binreg.py {site}'
    
    # Define the content of the sbatch script
    text = f'''#!/bin/bash
#SBATCH --job-name=changep_time_{site}    # Job name that will appear in squeue
#SBATCH --time=24:00:00             # Maximum runtime (DD-HH:MM:SS)
#SBATCH --ntasks=1                  # Number of tasks (usually 1 for a single script)
#SBATCH --mem-per-cpu=30gb          # Memory per CPU (adjust as needed)
#SBATCH --output={sbatch_dir}/changep_time_{site}_%j.out    # Standard output and error log (incorporates job ID %j)
#SBATCH --error={sbatch_dir}/changep_time_{site}_%j.err     # Standard error log
#SBATCH --mail-user=tbellagio@carnegiescience.edu # Email address for job notifications
#SBATCH --mail-type=FAIL            # Send email on job failure

# Load the necessary conda environment
source /home/tbellagio/miniforge3/etc/profile.d/conda.sh
conda activate /home/tbellagio/miniforge3/envs/pipeline_snakemake

# Change to the directory where the script and data are located on the cluster
# Assumes the script is run from /carnegie/nobackup/scratch/tbellagio/gea_grene-net/snp_origin
cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/snp_origin

# Execute the python script
{cmd}

'''
    
    # Write the sbatch script content to the file
    with open(file, 'w') as o:
        o.write("%s" % text)
        
    # Add the generated sbatch filename to the list
    shfiles.append(file)

print(f"Generated {len(shfiles)} sbatch files in the {sbatch_dir} directory.")

# --- Submit Jobs to SLURM ---

# Submit each generated sbatch file to the SLURM scheduler
print("Submitting jobs to SLURM...")
for shfile in shfiles:
    try:
        # Use subprocess.run to execute the sbatch command
        result = subprocess.run(["sbatch", shfile], check=True, capture_output=True, text=True)
        print(f"Submitted {shfile}: {result.stdout.strip()}")
    except subprocess.CalledProcessError as e:
        print(f"Error submitting {shfile}: {e.stderr.strip()}")
    except FileNotFoundError:
        print(f"Error: sbatch command not found. Ensure SLURM is in your PATH.")

print("Job submission process completed.")
