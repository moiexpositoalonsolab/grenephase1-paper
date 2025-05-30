#!/usr/bin/env python
# coding: utf-8

"""
Launch WZA Analysis for Binomial Regression Results (Bio12)

This script performs WZA (Weighted-Z Analysis) on binomial regression results for bio12.
It processes the regression results, creates visualizations, and submits the analysis to the cluster.

Required Files:
1. Input Files:
   - ../key_files/var_pos_grenenet.csv: SNP position information
   - ../key_files/blocks_snpsid_dict.pkl: Dictionary mapping SNP IDs to block IDs
   - ../key_files/maf_all_samples_last_gen.csv: Minor allele frequencies
   - ../binomial_regression_lastgen/bio12_binomial_reg_results_last_gen.csv: Binomial regression results

2. Script Files:
   - general_WZA_script_mod_polynomial_order7.py: WZA analysis script

Outputs:
1. Analysis Results:
   - bio12_binomial_reg_lastgen_wmaf.csv: Processed regression results with MAF
   - wza_binomial_regression_bio12_poly7.csv: WZA analysis results
   - wza.sh: SLURM job script

2. Visualization Files:
   - QQ plot showing observed vs expected p-values
   - Manhattan plot showing significant regions

Script Outline:
1. Load SNP and block mapping data
2. Process binomial regression results for bio12
3. Create and submit SLURM job for WZA analysis
4. Generate QQ and Manhattan plots
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
snps_names = pd.read_csv('../key_files/var_pos_grenenet.csv')

with open('../key_files/blocks_snpsid_dict.pkl', 'rb') as f:
    dict_blocks = pickle.load(f)

reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}
snps_names = snps_names[snps_names['total_alleles05filter_firstgen'].notna()].reset_index(drop=True)

# --- Process Binomial Regression Results ---
print("Processing binomial regression results...")
binomial_reg = pd.read_csv('../binomial_regression_lastgen/bio12_binomial_reg_results_last_gen.csv')

# Map SNP IDs to blocks
binomial_reg['block'] = binomial_reg['snp_id'].map(reverse_mapping)

# Add MAF information
maf = pd.read_csv('../key_files/maf_all_samples_last_gen.csv')
binomial_reg = pd.concat([binomial_reg, maf], axis=1)
binomial_reg.columns = ['slope', 'pvalue', 'snp_id', 'block', 'MAF']

# Save processed results
binomial_reg.to_csv('../binomial_regression_lastgen/bio12_binomial_reg_lastgen_wmaf.csv', index=None)

# --- Create and Submit SLURM Job ---
print("Creating SLURM job script...")
path = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/wza_last_gen/'
shfiles = []

seed = random.randint(1, 100000000)
file = 'wza.sh'
cmd = f'''python general_WZA_script_mod_polynomial_order7.py \\
    --correlations ../binomial_regression_lastgen/bio12_binomial_reg_lastgen_wmaf.csv \\
    --summary_stat pvalue \\
    --window "block" \\
    --output wza_binomial_regression_bio12_poly7.csv \\
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
cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/wza_last_gen
{cmd}
'''

with open(file, 'w') as o:
    o.write(text)
shfiles.append(file)

# Submit job
print("Submitting SLURM job...")
for shfile in shfiles:
    subprocess.run(["sbatch", shfile], check=True)

# --- Generate Visualizations ---
print("Generating visualizations...")
wza = pd.read_csv('wza_binomial_regression_bio12_poly7.csv')
wza = wza.dropna()

# Process data for plotting
wza['chrom'] = wza['gene'].str.split('_').str[0].astype(int)
wza['pos'] = wza['gene'].str.split('_').str[1].astype(int)

# Create QQ plot
print("Creating QQ plot...")
observed_quantiles = -np.log10(np.sort(wza['Z_pVal'].values))
expected_quantiles = -np.log10(np.linspace(1 / len(wza), 1, len(wza)))

plt.figure(figsize=(10, 10))
sns.scatterplot(x=expected_quantiles, y=observed_quantiles, edgecolor='b', facecolor='none', alpha=0.5)
plt.plot([min(expected_quantiles), max(expected_quantiles)], 
         [min(expected_quantiles), max(expected_quantiles)], 'r--')
plt.xlabel("Expected -log10(p-values)")
plt.ylabel("Observed -log10(p-values)")
plt.title('QQ Plot for Bio12 WZA')
plt.savefig('qq_plot_bio12.png')
plt.close()

# Create Manhattan plot
print("Creating Manhattan plot...")
df = wza[['Z_pVal', 'pos', 'chrom']].copy()
df['chromosome'] = df['chrom']
df['position'] = df['pos']
df['-log10(pvalue)'] = -np.log10(df['Z_pVal'])

colors = sns.color_palette("crest", n_colors=5)
chromosome_offsets = {}
offset = 0

for chrom in sorted(df['chromosome'].unique()):
    chromosome_offsets[chrom] = offset
    max_position = df[df['chromosome'] == chrom]['position'].max()
    offset += max_position + 200

df['adjusted_position'] = df.apply(
    lambda row: row['position'] + chromosome_offsets[row['chromosome']], 
    axis=1
)

plt.figure(figsize=(20, 6))
for chrom in sorted(df['chromosome'].unique()):
    subset = df[df['chromosome'] == chrom]
    plt.scatter(
        subset['adjusted_position'],
        subset['-log10(pvalue)'],
        alpha=0.7,
        c=colors[chrom % len(colors)],
        label=f'Chr {chrom}',
        s=20
    )

plt.xlabel('Position')
plt.ylabel('-log10(p-value)')
plt.title('Bio12 Manhattan Plot')
plt.grid(axis='y')
plt.tight_layout()
plt.savefig('manhattan_bio12.png')
plt.close()

print("Script execution completed.")
