#!/usr/bin/env python
# coding: utf-8

# This script processes the results from binomial regression analysis of allele frequency changes
# across different sites. It reads JSON results files, converts them to CSV format, and identifies
# statistically significant results using a Bonferroni-corrected p-value threshold.

# --- Required Input Files ---
# Ensure these files are present relative to the script's location:
# - ../key_files/var_pos_grenenet.csv (SNP position information)
# - ../key_files/merged_hapFIRE_allele_frequency.txt (Sample information)
# - binom_reg/*.jsonl (Results from binomial regression analysis)
# The script will write output to:
# - binom_reg/site_{site}.csv (Processed results for each site)

# --- Import Libraries ---
import pandas as pd
import os
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# --- Load and Prepare SNP Position Data ---
# Read SNP positions and filter based on MAF threshold
var_pos = pd.read_csv('../key_files/var_pos_grenenet.csv')
var_pos = var_pos[var_pos['maf05filter'].notna()].reset_index(drop=True)

# --- Process Sample Information ---
# Read sample information to identify sites and plots with sufficient data
merged_hapFIRE_allele_frequency = pd.read_csv('../key_files/merged_hapFIRE_allele_frequency.txt', nrows=1, sep='\t')
merged_hapFIRE_allele_frequency = merged_hapFIRE_allele_frequency.T.reset_index()

# Extract site, generation, and plot information from sample names
merged_hapFIRE_allele_frequency['site'] = merged_hapFIRE_allele_frequency['index'].str.split('_').str[0]
merged_hapFIRE_allele_frequency['gen'] = merged_hapFIRE_allele_frequency['index'].str.split('_').str[1]
merged_hapFIRE_allele_frequency['plot'] = merged_hapFIRE_allele_frequency['index'].str.split('_').str[2]

# Identify sites and plots with at least 2 years of data
samples = merged_hapFIRE_allele_frequency.drop([0, 'index'], axis=1).drop_duplicates()
samples_at_least_2_years = samples[['site', 'plot']].value_counts()[samples[['site', 'plot']].value_counts() > 1].reset_index()
sites_plots_w_at_least_2_years = samples_at_least_2_years[['site', 'plot']].groupby('site')['plot'].unique().to_dict()

# --- Process Regression Results ---
# Find all JSON result files
result_files = [i for i in os.listdir('binom_reg/') if 'json' in i]

# Process each result file
for site_file in result_files:
    # Extract site identifier from filename
    site = site_file.split('_')[-2]
    print(f"Processing results for site {site}")
    
    # Read and parse JSON results
    results = []
    with open('binom_reg/' + site_file, 'r') as file:
        for line in file:
            data = json.loads(line)
            results.append(data)
    
    # Convert results to DataFrame and merge with SNP positions
    df = pd.DataFrame(results)
    df = pd.concat([df[['slope', 'p_value']], var_pos['id']], axis=1)
    
    # Save processed results
    df.to_csv(f'binom_reg/site_{site}.csv', index=None)

# --- Identify Significant Results ---
# Find all CSV result files
result_files = [i for i in os.listdir('binom_reg') if 'csv' in i]

# Process each result file to identify significant results
for i in result_files:
    site = pd.read_csv(f'binom_reg/{i}')
    
    # Calculate Bonferroni-corrected threshold
    th = 0.05/len(site)
    
    # Identify significant results
    sign_site = site[site['p_value'] < th]
    
    # Print proportion of significant results
    print(f"Site {i}: {len(sign_site)/len(site):.2%} of SNPs are significant")
