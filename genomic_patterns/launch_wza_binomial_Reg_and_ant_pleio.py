#!/usr/bin/env python
# coding: utf-8

# This script launches the Weighted Z-score Analysis (WZA) on the results of the binomial regression
# performed by the previous script (`calc_delta_time_allplots_binreg.py`).
# It generates `sbatch` files for submitting WZA jobs to a Cedar server, merges the WZA results
# with the binomial regression results and block origin information, and creates a hexbin plot grid
# to visualize the relationship between SNP origin and allele frequency change across different sites.

# --- Required Input Files ---
# Ensure these files are present relative to the script's location:
# - ../key_files/var_pos_grenenet.csv
# - ../key_files/blocks_snpsid_dict.pkl
# - ../key_files/dict_snps.csv
# - ../key_files/block_origin_bio1_1001g.csv
# - ../key_files/bioclimvars_experimental_sites_era5.csv
# - Files in the 'binom_reg/' directory with the pattern 'site_*.csv'
# The script will also write output files to the 'binom_reg/' directory and the current directory.

# --- Import Libraries ---

import pandas as pd
import os
import random
import subprocess
import dask.dataframe as dd # Imported but not used in the provided code
import seaborn as sns
import numpy as np
import statsmodels.api as sm # Imported but not used in the provided code
import matplotlib.pyplot as plt
import pickle
import math

# --- Data Loading and Preparation ---

# Read the list of significant site results from the 'binom_reg/' directory.
sign_results = [i for i in os.listdir('binom_reg/') if 'sign' in i]
sites = pd.Series(sign_results).str.split('_').str[-1].str.replace('.csv', '').astype(int).reset_index()

# Read SNP names and positions.
snps_names = pd.read_csv('../key_files/var_pos_grenenet.csv')

# Read the dictionary mapping SNPs to blocks.
with open('../key_files/blocks_snpsid_dict.pkl', 'rb') as f:
    dict_blocks = pickle.load(f)

# Create a reverse mapping from SNP ID to block ID.
reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}

# Read the dictionary of SNPs.
dict_snps = pd.read_csv('../key_files/dict_snps.csv')

# --- Prepare Binomial Regression Results for WZA ---

# Iterate through the list of significant sites.
for i in sites[0].values:
    # Read the binomial regression results for the current site.
    binomial_reg = pd.read_csv(f'binom_reg/site_{i}.csv')
    # Add a 'block' column based on the reverse mapping.
    binomial_reg['block'] = binomial_reg['id'].map(reverse_mapping)
    # Set a default MAF value.
    binomial_reg['MAF'] = 0.06
    # Save the modified results to a new CSV file.
    binomial_reg.to_csv(f'binom_reg/site_{i}_wmaf.csv', index=None)

# --- Generate SBATCH Files for WZA ---

# Create a list to store the names of the generated `sbatch` files.
shfiles = []

# Iterate through the list of significant sites.
for site in sites[0].values:
    # Generate a random seed for the WZA analysis.
    seed = random.randint(1, 100000000)
    # Define the name of the `sbatch` file.
    file = f'wza_{site}.sh'
    # Define the command to run the WZA analysis.
    cmd = f'python general_WZA_script_mod.py --correlations binom_reg/site_{site}_wmaf.csv --summary_stat p_value --window block --output binom_reg/wza_site_{site}.csv --sep "," --sample_snps 1729 --resample 50'
    # Define the content of the `sbatch` file.
    text = f'''#!/bin/bash
#SBATCH --job-name=wza_{site}
#SBATCH --time=1:00:00  # Time limit set to 1 hour
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --output=wza_{site}_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

module load python/3.11_conda
conda activate /home/tbellagio/miniforge3/envs/pipeline_snakemake
export LD_LIBRARY_PATH="/home/tbellagio/miniforge3/envs/run_baypass/lib:$LD_LIBRARY_PATH"
cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/snp_origin
{cmd}
'''
    # Write the content of the `sbatch` file to the file.
    with open(file, 'w') as o:
        o.write("%s" % text)
    # Add the name of the `sbatch` file to the list.
    shfiles.append(file)

# --- Submit SBATCH Files to SLURM Scheduler ---

# Iterate through the list of `sbatch` files.
for shfile in shfiles:
    # Submit each `sbatch` file to the SLURM scheduler.
    subprocess.run(["sbatch", shfile], check=True)

# --- Merge WZA Results with Binomial Regression Results and Block Origin Information ---

# Read the block origin information from a CSV file.
block_origin_bio1 = pd.read_csv('../key_files/block_origin_bio1_1001g.csv')

# Read the list of significant site results and the corresponding site numbers.
sign_results = [i for i in os.listdir('binom_reg/') if 'sign' in i]
sites = pd.Series(sign_results).str.split('_').str[-1].str.replace('.csv', '').astype(int).reset_index()
sites.columns = ['ign', 'site']

# Read the climate data for the experimental sites and merge it with the site numbers.
clim_sites_during_exp = pd.read_csv('../key_files/bioclimvars_experimental_sites_era5.csv')
clim_sites_during_exp = clim_sites_during_exp[['site', 'bio1']]
clim_sites_during_exp = clim_sites_during_exp.merge(sites['site'])
clim_sites_during_exp_ordered = clim_sites_during_exp.sort_values('bio1')
unique_sites = clim_sites_during_exp_ordered['site'].values

# Iterate through the unique site numbers.
for i in unique_sites:
    print(i)
    # Read the WZA results for the current site.
    wza = pd.read_csv(f'binom_reg/wza_site_{i}.csv')
    # Read the binomial regression results for the current site.
    result_site = pd.read_csv(f'binom_reg/site_{i}.csv')
    # Add a 'block' column based on the reverse mapping.
    result_site['block'] = result_site['id'].map(reverse_mapping)
    # Calculate the median slope for each block.
    slopes = result_site.groupby('block')['slope'].median().reset_index()
    # Merge the WZA results with the median slopes and the block origin information.
    wza = wza.merge(slopes, left_on='gene', right_on='block').merge(block_origin_bio1, left_on='gene', right_on='blocks')
    # Save the merged results to a new CSV file.
    wza.to_csv(f'binom_reg/wza_site_{i}_pr.csv', index=None)

# --- Create Hexbin Plot Grid ---

# Define a function to convert log odds to probabilities.
def log_odds_to_probability(log_odds):
    return 1 / (1 + np.exp(-log_odds))

# Set the font type for PDF output.
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

# Create a grid of hexbin plots, one for each site.
n_plots = 6  # Number of subplots
fig, axes = plt.subplots(nrows=n_plots, ncols=4, figsize=(12, 2 * n_plots))
axes = axes.flatten()

# Iterate through the unique site numbers.
for i, (ax, site_value) in enumerate(zip(axes, unique_sites)):
    # Read the merged WZA and binomial regression results for the current site.
    site = pd.read_csv(f'binom_reg/wza_site_{site_value}_pr.csv')
    print(site_value)

    # Convert log odds to probabilities.
    site['slope'] = log_odds_to_probability(site['slope'])

    # Create a hexbin plot for the current site.
    hb = ax.hexbin(
        site['snp_origin_bio1'],
        site['slope'],
        gridsize=35,
        cmap='Greys_r',
        bins='log',
        mincnt=1
    )

    # Set title for the current subplot.
    ax.set_title(f"Site: {site_value}")

    # Add dotted grey horizontal and vertical lines at y=0.5 and x=0.
    ax.axhline(0.5, color='grey', linestyle='--', linewidth=1)
    ax.axvline(0, color='grey', linestyle='--', linewidth=1)

    # Set axis labels for the current subplot.
    ax.set_xlabel("SNP origin (mean annual temperature)")
    ax.set_ylabel("Allele Frequency")

# Remove any unused subplots.
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

# Configure the last subplot for the colorbar.
cbar_ax = axes[-1]
fig.colorbar(hb, cax=cbar_ax, orientation='vertical', label='Density')
cbar_ax.set_visible(False)

# Adjust layout and save the figure.
plt.tight_layout()
plt.savefig('block_origin_across_time_hexbin.pdf')
plt.show()

