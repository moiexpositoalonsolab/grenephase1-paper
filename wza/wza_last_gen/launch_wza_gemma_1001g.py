#!/usr/bin/env python
# coding: utf-8

"""
Launch WZA Analysis for GEMMA 1001G Results

This script performs WZA (Weighted-Z Analysis) on GEMMA GWAS results from the 1001G dataset.
It processes the GEMMA results, creates Manhattan plots, and submits the analysis to the cluster.

Required Files:
1. Input Files:
   - ../key_files/var_pos_grenenet.csv: SNP position information
   - ../key_files/blocks_snpsid_dict.pkl: Dictionary mapping SNP IDs to block IDs
   - /carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/allele_assoc_runs/lmm_gemma/{biovar}/output/results_lmm.csv: GEMMA results

2. Script Files:
   - general_WZA_script_mod_polynomial_order7.py: WZA analysis script

Outputs:
1. Analysis Results:
   - gemma_1001_wza_{biovar}.csv: WZA analysis results for each bioclimatic variable
   - wza_{biovar}.sh: SLURM job scripts for each bioclimatic variable

2. Visualization Files:
   - Manhattan plots for each bioclimatic variable showing significant regions

Script Outline:
1. Load SNP and block mapping data
2. Process GEMMA results for each bioclimatic variable
3. Create and submit SLURM jobs for WZA analysis
4. Generate Manhattan plots with significance thresholds
5. Apply multiple testing corrections
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
import matplotlib.gridspec as gridspec
from statsmodels.stats.multitest import multipletests
import pickle

# --- Load SNP and Block Mapping ---
print("Loading SNP and block mapping...")
snps_names = pd.read_csv('../key_files/var_pos_grenenet.csv')

with open('../key_files/blocks_snpsid_dict.pkl', 'rb') as f:
    dict_blocks = pickle.load(f)

reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}
snps_names = snps_names[snps_names['total_alleles05filter_lastgen'].notna()].reset_index(drop=True)

# --- Define Bioclimatic Variables ---
biovars = [
    'bio1', 'bio12', 'bio18', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6',
    'bio7', 'bio8', 'bio9', 'bio10', 'bio11', 'bio13', 'bio14', 'bio15',
    'bio16', 'bio17'
]

# --- Process GEMMA Results ---
print("Processing GEMMA results...")
for biovar in biovars:
    path_gemma = f'/carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/allele_assoc_runs/lmm_gemma/{biovar}/output/results_lmm.csv'
    gemma_result = pd.read_csv(path_gemma)
    gemma_result.columns = ['rs', 'p_wald', 'beta', 'significant', 'blocks', 'MAF']
    path_gemma_n = f'/carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/allele_assoc_runs/lmm_gemma/{biovar}/output/results_lmm_maf.csv'
    gemma_result.to_csv(path_gemma_n, index=None)

# --- Create and Submit SLURM Jobs ---
print("Creating SLURM job scripts...")
path = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/wza_last_gen/'
shfiles = []

for biovar in biovars:
    seed = random.randint(1, 100000000)
    file = f'wza_{biovar}.sh'
    cmd = f'''python general_WZA_script_mod_polynomial_order7.py \\
        --correlations /carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/allele_assoc_runs/lmm_gemma/{biovar}/output/results_lmm_maf.csv \\
        --summary_stat p_wald \\
        --window blocks \\
        --output gemma_1001_wza_{biovar}.csv \\
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

# Submit jobs
print("Submitting SLURM jobs...")
for shfile in shfiles:
    subprocess.run(["sbatch", shfile], check=True)

# --- Generate Manhattan Plots ---
print("Generating Manhattan plots...")
sns.set_context("talk")

for biovar in biovars:
    wza = pd.read_csv(f'gemma_231_wza_{biovar}.csv')
    
    # Process data for plotting
    wza['chrom'] = wza['gene'].str.split('_').str[0].astype(int)
    wza['pos'] = wza['gene'].str.split('_').str[1].astype(int)
    
    # Calculate quantiles for QQ plot
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
    
    # Set up colors and chromosome offsets
    colors = sns.color_palette("crest", n_colors=5)
    chromosome_offsets = {}
    offset = 0
    chrom_ends = {}
    
    for chrom in sorted(wza['chrom'].unique()):
        chromosome_offsets[chrom] = offset
        max_position = wza[wza['chrom'] == chrom]['pos'].max()
        offset += max_position + 200
        chrom_ends[offset] = (chrom, max_position)
    
    # Create Manhattan plot
    df['adjusted_position'] = df.apply(
        lambda row: row['position'] + chromosome_offsets[row['chromosome']], 
        axis=1
    )
    
    # Create figure with custom subplot widths
    fig = plt.figure(figsize=(20, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1])
    
    # Plot Manhattan plot
    ax1 = plt.subplot(gs[0])
    for chrom in sorted(df['chromosome'].unique()):
        subset = df[df['chromosome'] == chrom]
        ax1.scatter(
            subset['adjusted_position'],
            subset['-log10(pvalue)'],
            alpha=0.7,
            c=colors[chrom % len(colors)],
            label=f'Chr {chrom}',
            s=20
        )
    
    # Add significance threshold
    ax1.axhline(y=-np.log10(critical_pvalue), color='grey', linestyle='dashed')
    
    # Add aesthetics
    ax1.set_xlabel('Position')
    ax1.set_ylabel('-log10(p-value)')
    ax1.set_title(f'{biovar} Manhattan Plot')
    ax1.grid(axis='y')
    
    # Save plot
    plt.tight_layout()
    plt.savefig(f'manhattan_{biovar}.png')
    plt.close()

print("Script execution completed.")


# In[ ]:




