#!/usr/bin/env python
# coding: utf-8

"""
Process BSLMM Results for Climate GWAS

This script processes the output files from BSLMM analysis for each bioclimatic variable.
It identifies significant SNPs based on gamma values, maps them to genomic blocks,
and generates Manhattan plots showing the distribution of gamma values.

Required Files:
- blocks_snpsid_dict.pkl: Dictionary mapping SNP IDs to block IDs
- BSLMM output files:
  - /bslmm/{biovar}/output/significant_gammas.csv
  - /bslmm/{biovar}/output/{biovar}.param.txt
- Clumping results: /lfmm_full/clumping/output/output_clumping_{biovar}.clumped

Outputs:
- For each bioclimatic variable:
  - Processed results in CSV format with significance status and block information
  - Manhattan plot showing gamma values and clumped SNPs
"""

# --- Import Libraries ---
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import allel
import seaborn as sns
import os
import pickle
import dask.dataframe as dd

# --- Load Block Mapping Dictionary ---
print("Loading block mapping dictionary...")
with open('/home/tbellagio/HapFM/blocks_snpsid_dict.pkl', 'rb') as file:
    dict_blocks = pickle.load(file)

# Create inverse dictionary mapping SNP ID to Block ID
dict_blocks_inv = {}
for key, values in dict_blocks.items():
    for value in values:
        dict_blocks_inv[value] = key

# --- Set Up Paths ---
path = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/allele_assoc_runs'
path_bslmm = path + '/bslmm'

# Get list of bioclimatic variables
biovars = [i for i in os.listdir(path_bslmm) if 'bio' in i]

# --- Process BSLMM Results for Each Bioclimatic Variable ---
print("Processing BSLMM results for each bioclimatic variable...")
for biovar in biovars:
    print(f"Processing {biovar}...")
    
    # Read significant gammas and parameter files
    significant_gammas = pd.read_csv(path_bslmm + f'/{biovar}/output/significant_gammas.csv')
    params = pd.read_csv(path_bslmm + f'/{biovar}/output/{biovar}.param.txt', sep='\t')
    
    # Process significant gammas
    significant_gammas.columns = ['significant']
    params = params.merge(significant_gammas, how='left', left_on='gamma', right_on='significant')
    params['significant'] = params['significant'].notna()
    
    # Map SNPs to blocks
    params['blocks'] = params['rs'].map(dict_blocks_inv)
    
    # Save processed results
    params[['rs', 'gamma', 'beta', 'significant', 'blocks']].to_csv(
        path_bslmm + f'/{biovar}/output/results_bslmm.csv',
        index=None
    )

# --- Generate Manhattan Plots ---
print("Generating Manhattan plots...")
for biovar in biovars:
    print(f"Generating plot for {biovar}...")
    
    # Read clumping results
    clump = pd.read_csv(
        f'/carnegie/nobackup/scratch/tbellagio/gea_grene-net/lfmm_full/clumping/output/output_clumping_{biovar}.clumped',
        sep='\s+'
    )
    
    # Read processed BSLMM results
    pvalues_file = path_bslmm + f'/{biovar}/output/results_bslmm.csv'
    pvalues = dd.read_csv(pvalues_file).compute()
    
    # Prepare data for plotting
    df = pvalues.copy()
    colors = sns.color_palette("crest", n_colors=5)
    
    # Parse chromosome and position information
    df['chromosome'] = df['rs'].str.split('_').str[0].astype(int)
    df['position'] = df['rs'].str.split('_').str[1].astype(int)
    
    # Calculate chromosome offsets for adjusted positions
    chromosome_offsets = {}
    offset = 0
    for chrom in sorted(df['chromosome'].unique()):
        chromosome_offsets[chrom] = offset
        max_position = df[df['chromosome'] == chrom]['position'].max()
        offset += max_position + 1000000  # Buffer between chromosomes
    
    # Apply offsets to positions
    df['adjusted_position'] = df.apply(
        lambda row: row['position'] + chromosome_offsets[row['chromosome']],
        axis=1
    )
    
    # Create Manhattan plot
    plt.figure(figsize=(20, 6))
    
    # Plot all SNPs
    for chrom in sorted(df['chromosome'].unique()):
        subset = df[df['chromosome'] == chrom]
        plt.scatter(
            subset['adjusted_position'],
            subset['gamma'],
            c=colors[chrom % len(colors)],
            label=f'Chr {chrom}',
            s=10
        )
    
    # Highlight clumped SNPs
    clumped_subset = df[df['significant'] == True]
    plt.scatter(
        clumped_subset['adjusted_position'],
        clumped_subset['gamma'],
        s=50,
        facecolors='none',
        edgecolors='grey',
        linewidths=2,
        label='Clumped SNPs'
    )
    
    # Add aesthetics
    plt.xlabel('Adjusted Position')
    plt.ylabel('gamma')
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.title(f'{biovar}')
    
    # Save plot
    plt.tight_layout()
    plt.savefig(path_bslmm + f'/{biovar}/output/manhattan.png')
    plt.close()

print("Script execution completed.")


# In[ ]:




