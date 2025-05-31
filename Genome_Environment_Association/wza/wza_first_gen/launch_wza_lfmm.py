#!/usr/bin/env python
# coding: utf-8

"""
Launch WZA Analysis for LFMM Results

This script prepares and launches WZA (Weighted-Z Analysis) for LFMM (Latent Factor Mixed Model) results.
It processes the LFMM results, adds MAF information, and creates visualization plots.

Required Files:
1. Input Files:
   - ../data/blocks_snpsid_dict.pkl: Dictionary mapping SNP IDs to block IDs
   - ../data/maf_all_samples.csv: Minor allele frequencies
   - ../key_files/p0_average_seed_mix.csv: Seed mix data
   - ../key_files/var_pos_grenenet.csv: SNP position information
   - ./lfmm_full/lfmm_fullresults_all_k/lfmm_{biovar}_k25_results.csv: LFMM results

2. Script Files:
   - general_WZA_script_mod_polynomial_order7.py: WZA analysis script

Outputs:
1. Processed Files:
   - lfmm_{biovar}_results_all_samples_k25_wmaf.csv: LFMM results with MAF
   - wza_results_lfmm_{biovar}_poly7.csv: WZA analysis results

2. Visualization Files:
   - qq_plot_{biovar}.png: QQ plot of p-values
   - manhattan_lfmm_{biovar}.png: Manhattan plot of results

Script Outline:
1. Load and prepare SNP and block mapping data
2. Process LFMM results and add MAF information
3. Create and submit SLURM job for WZA analysis
4. Generate QQ plot
5. Create Manhattan plot with significance threshold
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
import pickle

# --- Load SNP and Block Mapping ---
print("Loading SNP and block mapping...")
dict_blocks = '../data/blocks_snpsid_dict.pkl'

with open(dict_blocks, 'rb') as file:
    dict_blocks = pickle.load(file)

reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}

# --- Load MAF Data ---
print("Loading MAF data...")
maf = pd.read_csv('../data/maf_all_samples.csv')
maf.columns = ['MAF']

# --- Set Working Directory and Bioclimatic Variable ---
wd = './lfmm_full/lfmm_fullresults_all_k/'
biovar = 'bio18'  # Can be modified for different bioclimatic variables

# --- Process LFMM Results ---
print(f"Processing LFMM results for {biovar}...")
lfmm = f'./lfmm_full/lfmm_fullresults_all_k/lfmm_{biovar}_k25_results.csv'
lfmm = pd.read_csv(lfmm)

# Add MAF information
lfmm = pd.concat([lfmm, maf['MAF']], axis=1)
lfmm.to_csv(f'lfmm_{biovar}_results_all_samples_k25_wmaf.csv', index=None)

# --- Create and Submit SLURM Job ---
print("Creating SLURM job script...")
path = './wza'
shfiles = []

# Create SLURM job script
seed = random.randint(1, 100000000)
file = 'wza.sh'
cmd = f'''python general_WZA_script_mod_polynomial_order7.py \\
        --correlations lfmm_{biovar}_results_all_samples_k25_wmaf.csv \\
        --summary_stat p_value \\
        --window block \\
        --output wza_results_lfmm_{biovar}_poly7.csv \\
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

# --- Process WZA Results ---
print("Loading and processing WZA results...")
wza_lfmm = pd.read_csv(f'wza_results_lfmm_{biovar}_poly7.csv')
wza_lfmm = wza_lfmm.dropna()

# Load additional data for analysis
p0_average_seed_mix = pd.read_csv('../key_files/p0_average_seed_mix.csv')
snps_dict = pd.read_csv('../key_files/var_pos_grenenet.csv')
p0_average_seed_mix = pd.concat([p0_average_seed_mix, snps_dict], axis=1)

# --- Generate QQ Plot ---
print("Generating QQ plot...")
observed_quantiles = -np.log10(np.sort(wza_lfmm['Z_pVal'].values))
expected_quantiles = -np.log10(np.linspace(1 / len(wza_lfmm), 1, len(wza_lfmm)))

plt.figure(figsize=(10, 10))
sns.scatterplot(x=expected_quantiles, y=observed_quantiles, edgecolor='b', facecolor='none', alpha=0.5)
plt.plot([min(expected_quantiles), max(expected_quantiles)], 
         [min(expected_quantiles), max(expected_quantiles)], 'r--')

plt.xlabel("Expected -log10(p-values)")
plt.ylabel("Observed -log10(p-values)")
plt.title(f'QQ Plot for {biovar} WZA')
plt.savefig(f'qq_plot_{biovar}.png')
plt.close()

# --- Create Manhattan Plot ---
print("Creating Manhattan plot...")
wza_lfmm['chrom'] = wza_lfmm['gene'].str.split('_').str[0].astype(int)
wza_lfmm['pos'] = wza_lfmm['gene'].str.split('_').str[1].astype(int)

# Calculate significance threshold
threshold_value = 0.05 / len(wza_lfmm)
threshold = -np.log10(threshold_value)

# Prepare data for plotting
wza_lfmm['chrom_pos'] = wza_lfmm['chrom'].astype(str) + '_' + wza_lfmm['pos'].astype(str)
df = wza_lfmm[['Z_pVal', 'pos', 'chrom', 'chrom_pos']].copy()

# Create Manhattan plot
plt.figure(figsize=(20, 6))
colors = sns.color_palette("crest", n_colors=5)

# Calculate chromosome offsets
chromosome_offsets = {}
offset = 0
for chrom in sorted(df['chrom'].unique()):
    chromosome_offsets[chrom] = offset
    max_position = df[df['chrom'] == chrom]['pos'].max()
    offset += max_position + 200  # Buffer between chromosomes

# Apply offsets to positions
df['adjusted_position'] = df.apply(
    lambda row: row['pos'] + chromosome_offsets[row['chrom']], 
    axis=1
)

# Plot each chromosome
for chrom in sorted(df['chrom'].unique()):
    subset = df[df['chrom'] == chrom]
    plt.scatter(
        subset['adjusted_position'],
        -np.log10(subset['Z_pVal']),
        alpha=0.7,
        c=colors[chrom % len(colors)],
        label=f'Chr {chrom}',
        s=20
    )

# Add aesthetics
plt.xlabel('Adjusted Position')
plt.ylabel('-log10(p-value)')
plt.title(f'{biovar} LFMM')
plt.grid(axis='y')

# Add significance threshold
plt.axhline(y=threshold, color='grey', linestyle='dashed')

# Save plot
plt.tight_layout()
plt.savefig(f'manhattan_lfmm_{biovar}.png')
plt.close()

print("Script execution completed.")


# In[ ]:




