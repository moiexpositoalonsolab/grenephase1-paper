#!/usr/bin/env python
# coding: utf-8

"""
Launch WZA Analysis for LFMM Results with Scale Checking

This script performs WZA (Weighted-Z Analysis) on LFMM results with additional scale checking
and Fisher's method for p-value aggregation. It includes statistical validation and visualization
of the results.

Required Files:
1. Input Files:
   - ../kendall_tau/kendall_corr_{biovar}.csv: Kendall Tau correlation results
   - before_filtering_wza_df.csv: Pre-filtered WZA results

2. Script Files:
   - general_WZA_script_mod_polynomial_order7.py: WZA analysis script

Outputs:
1. Analysis Results:
   - wza_results_lfmm_{biovar}_poly7check.csv: WZA analysis results with scale checking
   - fisher_aggregated_results.csv: Results of Fisher's method aggregation

Script Outline:
1. Import required libraries
2. Define Fisher's method for p-value aggregation
3. Process Kendall Tau correlation results
4. Create and submit SLURM job for WZA analysis
5. Perform Fisher's method aggregation
6. Generate statistical summaries
"""

# --- Import Libraries ---
import pandas as pd
import os 
import random
import subprocess
import dask.dataframe as dd
import seaborn as sns
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
import scipy.stats
from scipy.stats import norm, chi2
from pandas.api.types import is_string_dtype, is_numeric_dtype

# --- Define Fisher's Method Function ---
def fishers_method(p_values):
    """
    Applies Fisher's method to combine p-values into a single p-value.
    
    Parameters:
    p_values (list or np.array): A list or array of p-values for SNPs in a window.
    
    Returns:
    float: Aggregated p-value using Fisher's method.
    """
    # Ensure all p-values are above 0 (log of 0 is undefined)
    p_values = np.clip(p_values, 1e-10, 1)  # Avoid taking log(0)
    
    # Calculate the test statistic
    chi_square_stat = -2 * np.sum(np.log(p_values))
    
    # Degrees of freedom (2 * number of p-values)
    df = 2 * len(p_values)
    
    # Compute the p-value from the chi-square distribution
    combined_p_value = 1 - chi2.cdf(chi_square_stat, df)
    
    return combined_p_value

# --- Configuration and File Paths ---
# Path to input files
maf = pd.read_csv('../data/maf_all_samples.csv')
dict_blocks = '../data/blocks_snpsid_dict.pkl'

# Path to working directory
path = './wza'

# --- Set Parameters ---
biovar = 'bio1'  # Can be modified for different bioclimatic variables

# --- Load and Process Data ---
print(f"Loading Kendall Tau correlation results for {biovar}...")
kendall_corr = pd.read_csv(f'../kendall_tau/kendall_corr_{biovar}.csv')

# Get minimum p-value per block
print("Calculating minimum p-values per block...")
kendall_corr_top = kendall_corr.groupby('block')['K_tau_p'].min().reset_index()

# Apply Fisher's method
print("Applying Fisher's method for p-value aggregation...")
aggregated_pvals = kendall_corr.groupby('block')['K_tau_p'].apply(fishers_method)
fisher_results = aggregated_pvals.reset_index()
fisher_results.columns = ['block', 'fisher_agg']

# --- Create and Submit SLURM Job ---
print("Creating SLURM job script...")
shfiles = []

# Create SLURM job script
seed = random.randint(1, 100000000)
file = 'wza.sh'
cmd = f'''python general_WZA_script_mod_polynomial_order7.py \\
        --correlations ../kendall_tau/kendall_corr_{biovar}.csv \\
        --summary_stat K_tau_p \\
        --window block \\
        --output wza_results_lfmm_{biovar}_poly7check.csv \\
        --sep ","'''

text = f'''#!/bin/bash
#SBATCH --job-name=wza
#SBATCH --time=1:00:00  # Time limit set to 4 days
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --output=wza_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

module load python/3.11_conda
conda activate /home/tbellagio/miniforge3/envs/pipeline_snakemake
export LD_LIBRARY_PATH="/home/tbellagio/miniforge3/envs/run_baypass/lib:$LD_LIBRARY_PATH"
cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/wza
{cmd}
'''

with open(file, 'w') as o:
    o.write(text)
shfiles.append(file)

# Submit jobs
print("Submitting SLURM jobs...")
for shfile in shfiles:
    subprocess.run(["sbatch", shfile], check=True)

# --- Load WZA Results ---
print("Loading WZA results...")
wza_df = pd.read_csv('before_filtering_wza_df.csv')

# Change to the appropriate directory
os.chdir(path)

print("Script execution completed.")
