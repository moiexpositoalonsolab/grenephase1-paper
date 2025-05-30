#!/usr/bin/env python
# coding: utf-8

"""
Launch WZA Analysis for Kendall Tau Results

This script performs WZA (Weighted-Z Analysis) on Kendall Tau correlation results.
It processes the correlation results, creates visualizations, and submits the analysis to the cluster.

Required Files:
1. Input Files:
   - ../kendall_tau_last_gen/kendall_corr_{biovar}.csv: Kendall Tau correlation results for each bioclimatic variable

2. Script Files:
   - general_WZA_script_mod_polynomial_order7.py: WZA analysis script

Outputs:
1. Analysis Results:
   - wza_kendalltau_results_{biovar}_poly7.csv: WZA analysis results for each bioclimatic variable
   - wza_{biovar}.sh: SLURM job scripts

2. Visualization Files:
   - QQ plot showing observed vs expected p-values
   - Manhattan plot showing significant regions

Script Outline:
1. Process Kendall Tau correlation results for each bioclimatic variable
2. Create and submit SLURM jobs for WZA analysis
3. Generate QQ and Manhattan plots for each variable
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
from matplotlib import gridspec
from statsmodels.stats.multitest import multipletests

# --- Define Bioclimatic Variables ---
biovars = [
    'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9',
    'bio10', 'bio11', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio19'
]

# --- Create and Submit SLURM Jobs ---
print("Creating SLURM job scripts...")
path = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/wza_last_gen/'
shfiles = []

for biovar in biovars:
    seed = random.randint(1, 100000000)
    file = f'wza_{biovar}.sh'
    cmd = f'''python general_WZA_script_mod_polynomial_order7.py \\
        --correlations ../kendall_tau_last_gen/kendall_corr_{biovar}.csv \\
        --summary_stat K_tau_p \\
        --window "block" \\
        --output wza_kendalltau_results_{biovar}_poly7.csv \\
        --sep ","'''

    text = f'''#!/bin/bash
#SBATCH --job-name=wza_{biovar}
#SBATCH --time=1:00:00  # Time limit set to 4 days
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --output=wza_%j_{biovar}.out
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

# Submit jobs
print("Submitting SLURM jobs...")
for shfile in shfiles:
    subprocess.run(["sbatch", shfile], check=True)

# --- Generate Visualizations ---
print("Generating visualizations...")
sns.set_context("talk")

for biovar in biovars:
    print(f"Processing {biovar}...")
    wza = pd.read_csv(f'wza_kendalltau_results_{biovar}_poly7.csv')
    
    # Process data for plotting
    wza['chrom'] = wza['gene'].str.split('_').str[0].astype(int)
    wza['pos'] = wza['gene'].str.split('_').str[1].astype(int)
    
    # Calculate QQ plot data
    observed_quantiles = -np.log10(np.sort(wza['Z_pVal'].values))
    expected_quantiles = -np.log10(np.linspace(1 / len(wza), 1, len(wza)))
    
    # Apply Benjamini-Hochberg correction
    _, adjusted_pvals, _, _ = multipletests(wza['Z_pVal'], alpha=0.05, method='fdr_bh')
    wza['adjusted_pval'] = adjusted_pvals
    
    # Find critical p-value
    critical_pvalue = wza.loc[wza['adjusted_pval'] <= 0.05, 'Z_pVal'].max()
    
    # Prepare data for Manhattan plot
    df = wza[['Z_pVal', 'pos', 'chrom']].copy()
    df['chromosome'] = df['chrom']
    df['position'] = df['pos']
    df['-log10(pvalue)'] = -np.log10(df['Z_pVal'])
    
    # Calculate chromosome offsets
    colors = sns.color_palette("crest", n_colors=5)
    chromosome_offsets = {}
    offset = 0
    chrom_ends = {}
    
    for chrom in sorted(wza['chrom'].unique()):
        chromosome_offsets[chrom] = offset
        max_position = wza[wza['chrom'] == chrom]['pos'].max()
        offset += max_position + 200
        chrom_ends[offset] = (chrom, max_position)
    
    df['adjusted_position'] = df.apply(
        lambda row: row['position'] + chromosome_offsets[row['chromosome']], 
        axis=1
    )
    
    # Create figure with QQ and Manhattan plots
    fig = plt.figure(figsize=(20, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1])
    
    # Manhattan plot
    ax1 = plt.subplot(gs[0])
    for chrom in sorted(df['chromosome'].unique()):
        subset = df[df['chromosome'] == chrom]
        ax1.scatter(
            subset['adjusted_position'],
            subset['-log10(pvalue)'],
            alpha=0.7,
            c=colors[chrom % len(colors)],
            s=50
        )
    
    ax1.set_xlabel('Chromosomes')
    ax1.set_ylabel('-log10(pvalue)')
    ax1.set_title(f'{biovar} Kendall Tau + WZA')
    ax1.axhline(y=-np.log10(critical_pvalue), color='grey', linestyle='dashed')
    
    # QQ plot
    ax2 = plt.subplot(gs[1])
    sns.scatterplot(
        x=expected_quantiles,
        y=observed_quantiles,
        edgecolor='b',
        facecolor='none',
        alpha=0.5,
        ax=ax2
    )
    ax2.plot(
        [min(expected_quantiles), max(expected_quantiles)],
        [min(expected_quantiles), max(expected_quantiles)],
        'r--'
    )
    ax2.set_xlabel("Expected -log10(p-values)")
    ax2.set_ylabel("Observed -log10(p-values)")
    ax2.set_title('QQ Plot')
    
    plt.tight_layout()
    plt.savefig(f'manhattan_qq_{biovar}.png')
    plt.close()

print("Script execution completed.")


# In[ ]:




