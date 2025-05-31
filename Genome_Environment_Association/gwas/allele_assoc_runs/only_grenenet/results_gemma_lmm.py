#!/usr/bin/env python
# coding: utf-8

"""
Process GEMMA LMM Results for Climate GWAS (Only Grenenet Data)

This script processes the output files (.assoc.txt) from GEMMA LMM analysis for the
grenenet data subset. It identifies statistically significant SNPs, maps them to
genomic blocks, and generates Manhattan plots.

Required Files:
- var_pos_grenenet.csv: SNP position information
- blocks_snpsid_dict.pkl: Dictionary mapping SNP IDs to block IDs
- GEMMA output files: /only_grenenet/lmm_gemma/{biovar}/output/{biovar}.assoc.txt

Outputs:
- For each bioclimatic variable:
  - Processed results in CSV format with significance status and block information
  - Manhattan plot showing significant associations
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

# --- Load and Prepare Mapping Dictionaries ---
print("Loading and preparing mapping dictionaries...")
# Load SNP position information
dict_snps = pd.read_csv('../../data/var_pos_grenenet.csv')

with open('/home/tbellagio/HapFM/blocks_snpsid_dict.pkl', 'rb') as file:
    dict_blocks = pickle.load(file)

# Create inverse dictionary mapping SNP ID to Block ID
dict_blocks_inv = {}
for key, values in dict_blocks.items():
    for value in values:
        dict_blocks_inv[value] = key

# --- Set Up Paths ---
path = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/allele_assoc_runs/only_grenenet'
path_gemma = path + '/lmm_gemma'

# Get list of bioclimatic variables
biovars = [i for i in os.listdir(path_gemma) if 'bio' in i]

# --- Process GEMMA Results for Each Bioclimatic Variable ---
print("Processing GEMMA results for each bioclimatic variable...")
for biovar in biovars:
    print(f"Processing {biovar}...")
    
    # Read GEMMA association results
    assoc = pd.read_csv(path_gemma + f'/{biovar}/output/{biovar}.assoc.txt', sep='\t')
    
    # Calculate Bonferroni significance threshold
    th = 0.05 / len(assoc)
    
    # Add significance status
    assoc['significant'] = assoc['p_wald'] < th
    print(f"Significant SNPs in {biovar}:")
    print(assoc['significant'].value_counts())
    
    # Map SNPs to blocks
    assoc['blocks'] = assoc['rs'].map(dict_blocks_inv)
    
    # Save processed results
    assoc[['rs', 'p_wald', 'beta', 'significant', 'blocks', 'af']].to_csv(
        path_gemma + f'/{biovar}/output/results_lmm.csv', 
        index=None
    )

# --- Generate Manhattan Plots ---
print("Generating Manhattan plots...")
for biovar in biovars:
    print(f"Generating plot for {biovar}...")
    
    # Read processed results
    pvalues_file = path_gemma + f'/{biovar}/output/results_lmm.csv'
    pvalues = dd.read_csv(pvalues_file).compute()
    
    # Rename columns for clarity
    pvalues.columns = ['id', 'pvalue', 'beta', 'significant', 'blocks', 'af']
    
    # Calculate significance threshold
    threshold_value = 0.05 / len(pvalues)
    
    # Prepare data for plotting
    df = pvalues.copy()
    colors = sns.color_palette("crest", n_colors=5)
    
    # Parse chromosome and position information
    df['chromosome'] = df['id'].str.split('_').str[0].astype(int)
    df['position'] = df['id'].str.split('_').str[1].astype(int)
    df['-log10(pvalue)'] = -np.log10(df['pvalue'])
    
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
    
    for chrom in sorted(df['chromosome'].unique()):
        subset = df[df['chromosome'] == chrom]
        plt.scatter(
            subset['adjusted_position'], 
            subset['-log10(pvalue)'], 
            c=colors[chrom % len(colors)], 
            label=f'Chr {chrom}', 
            s=10
        )
    
    # Add aesthetics
    plt.xlabel('Adjusted Position')
    plt.ylabel('-log10(pvalue)')
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add significance threshold line
    threshold = -np.log10(threshold_value)
    plt.axhline(y=threshold, color='grey', linestyle='dashed')
    plt.title(f'{biovar}')
    
    # Save plot
    plt.tight_layout()
    plt.savefig(path_gemma + f'/{biovar}/output/manhattan.png')
    plt.close()

print("Script execution completed.")
