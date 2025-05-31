#!/usr/bin/env python3
# coding: utf-8

"""
Process GEMMA LMM Results, Identify Significant Loci, and Visualize

This script processes the output files (.assoc.txt) from GEMMA LMM analysis for fitness GWAS.
It identifies statistically significant SNPs (using Bonferroni correction), maps them to
genomic blocks, and performs downstream analyses and visualizations.

Key steps include:
- Reading GEMMA association results for each site.
- Identifying significant SNPs based on a multiple testing corrected p-value threshold.
- Mapping significant SNPs to predefined genomic blocks.
- Saving processed results including significance status and block information.
- Generating Manhattan plots to visualize GWAS results per site.
- Analyzing and visualizing the relationship between the number of significant loci and
  ecological factors (climate, parallelism, diversity).
- Identifying and quantifying shared significant SNPs and blocks between pairs of sites
  and relating this to climatic distance.

Inputs:
- ../../data/var_pos_grenenet.csv: SNP position information.
- ../../data/blocks_snpsid_dict.pkl: Pickle file containing a dictionary mapping SNP IDs to block IDs.
- site_*/output/*.assoc.txt: GEMMA association results files.
- ../../data/bioclimvars_experimental_sites_era5.csv: Climate data for experimental sites.
- ../../data/generation_1_parallelism.txt: Parallelism data.
- ../../data/unique_ecotypes.csv: Unique ecotype counts.
- ../../data/shannon_div.csv: Shannon diversity data.

Outputs:
- For each site: site_*/output/results_lmm.csv:
  Processed GEMMA results including significance status and block information.
- Interactive Manhattan plots for each site.
- Interactive scatter plots showing the relationship between the number of significant loci and
  ecological factors.
- Interactive scatter plot showing the relationship between shared significant blocks and
  climatic distance.
- A DataFrame (block_count_df) summarizing the frequency of each block being significant across site pairs.
"""

# --- Import Libraries ---
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr, kendalltau
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import allel
import seaborn as sns
import os
import pickle
import dask.dataframe as dd
from itertools import combinations
import statsmodels.api as sm

# --- Configuration and Data Paths ---
# Path to the base directory containing site-specific GEMMA results
gemma_results_base_path = '.'

# Path to SNP position information file
snp_pos_file = '../../data/var_pos_grenenet.csv'

# Path to the pickle file containing block-SNP mapping dictionary
blocks_dict_file = '../../data/blocks_snpsid_dict.pkl'

# Paths to ecological data files
climate_file = '../../data/bioclimvars_experimental_sites_era5.csv'
parallelism_file = '../../data/generation_1_parallelism.txt'
unique_ecotypes_file = '../../data/unique_ecotypes.csv'
shannon_div_file = '../../data/shannon_div.csv'

# --- Load and Prepare Mapping Dictionaries ---
print("Loading and preparing mapping dictionaries...")
# Load SNP position information
dict_snps = pd.read_csv(snp_pos_file)

# Load block-SNP mapping dictionary
with open(blocks_dict_file, 'rb') as file:
    dict_blocks = pickle.load(file)

# Create inverse dictionary mapping SNP ID to Block ID
dict_blocks_inv = {}
for block_id, snp_ids in dict_blocks.items():
    for snp_id in snp_ids:
        dict_blocks_inv[snp_id] = block_id

# --- Process GEMMA Results per Site ---
print("Processing GEMMA results per site...")

def list_only_directories(path):
    """
    Lists names of directories within a given path that contain 'site'.
    """
    all_entries = os.listdir(path)
    site_directories = [entry for entry in all_entries if 'site' in entry and os.path.isdir(os.path.join(path, entry))]
    return site_directories

# Get list of site directories
site_directories = list_only_directories(gemma_results_base_path)

# Iterate through each site directory to process GEMMA results
for site in site_directories:
    just_number_site = site.replace('site_', '')
    gemma_assoc_file = os.path.join(gemma_results_base_path, site, 'output', f'{just_number_site}.assoc.txt')
    output_csv_file = os.path.join(gemma_results_base_path, site, 'output', 'results_lmm.csv')
    
    if not os.path.exists(gemma_assoc_file):
        print(f"Warning: GEMMA association file not found for site {site}: {gemma_assoc_file}")
        continue
        
    print(f"Processing GEMMA results for site: {site}...")
    # Read GEMMA association results
    assoc_df = pd.read_csv(gemma_assoc_file, sep = '\t')
    
    # Calculate Bonferroni significance threshold
    bonferroni_threshold = 0.05 / len(assoc_df)
    
    # Determine significance and map to blocks
    assoc_df['significant'] = assoc_df['p_wald'] < bonferroni_threshold
    
    # Map SNP rs ID to block ID using the inverse dictionary
    assoc_df['blocks'] = assoc_df['rs'].map(dict_blocks_inv)
    
    # Rename frequency column for clarity and select relevant columns
    assoc_df = assoc_df.rename(columns = {'af': 'MAF'})[['rs','p_wald', 'beta', 'significant', 'blocks','MAF']]
    
    # Save processed results to CSV
    assoc_df.to_csv(output_csv_file, index=None)
    print(f"Processed results saved to: {output_csv_file}")

# --- Load Processed Results and Ecological Data for Downstream Analysis ---
print("Loading processed results and ecological data...")

# Load climate data
climate = pd.read_csv(climate_file)[['site', 'bio1']]
climate['site'] = climate['site'].astype(int) # Ensure site column is integer for merging

# Load parallelism data and filter for SNP source
parallelism = pd.read_csv(parallelism_file, sep = '\t')
parallelism_snps = parallelism[parallelism['source'] == 'snp'][['site', 'mean']]
parallelism_snps.columns = ['site', 'parallelism_snp'] # Rename column
parallelism_snps['site'] = parallelism_snps['site'].astype(int) # Ensure site column is integer

# Load unique ecotypes data and aggregate by site
unique_ecotypes = pd.read_csv(unique_ecotypes_file)
unique_ecotypes_persite = unique_ecotypes.groupby('site')['unique_ecotypes'].mean().reset_index()
unique_ecotypes_persite['site'] = unique_ecotypes_persite['site'].astype(int) # Ensure site column is integer

# Load Shannon diversity data and aggregate by site
shannon_div = pd.read_csv(shannon_div_file)
shannon_div = shannon_div.groupby('site')['0'].mean().reset_index() # Assuming column '0' contains diversity values
shannon_div.columns = ['site', 'shannon_diversity'] # Rename column
shannon_div['site'] = shannon_div['site'].astype(int) # Ensure site column is integer

# Aggregate significant SNP and block counts per site
sign_snps_number_dict = {}
sign_snps_names_dict = {}
sign_blocks_names_dict = {}

for site in site_directories:
    site_number = int(site.replace('site_', ''))
    processed_results_file = os.path.join(gemma_results_base_path, site, 'output', 'results_lmm.csv')
    
    if not os.path.exists(processed_results_file):
        print(f"Warning: Processed results file not found for site {site}: {processed_results_file}")
        continue
        
    # Read processed results using dask for potentially large files, then compute
    pvalues_df = dd.read_csv(processed_results_file).compute()
    
    # Count significant SNPs and collect their names and corresponding block names
    significant_snps_df = pvalues_df[pvalues_df['significant']==True]
    
    sign_snps_number_dict[site_number] = len(significant_snps_df)
    sign_snps_names_dict[site_number] = significant_snps_df['rs'].tolist() # Store as list
    sign_blocks_names_dict[site_number] = significant_snps_df['blocks'].dropna().tolist() # Store as list, drop NaNs

# Create a DataFrame for significant SNP counts and merge with ecological data
sign_snps_number_df = pd.DataFrame.from_dict(sign_snps_number_dict, orient='index', columns=['significant_snps_count']).reset_index()
sign_snps_number_df.columns = ['site', 'significant_snps_count']
sign_snps_number_df['site'] = sign_snps_number_df['site'].astype(int) # Ensure site column is integer

# Merge with ecological dataframes
sign_snps_analysis_df = sign_snps_number_df.merge(climate, on='site', how='left')
sign_snps_analysis_df = sign_snps_analysis_df.merge(parallelism_snps, on='site', how='left')
sign_snps_analysis_df = sign_snps_analysis_df.merge(unique_ecotypes_persite, on='site', how='left')
sign_snps_analysis_df = sign_snps_analysis_df.merge(shannon_div, on='site', how='left')

# --- Visualize Significant Loci vs. Ecological Factors ---
print("Generating plots for significant loci vs. ecological factors...")

# Plot Significant SNPs vs. bio1
plt.figure(figsize=(8, 6))
sns.scatterplot(data=sign_snps_analysis_df, x='bio1', y='significant_snps_count')
for i in range(sign_snps_analysis_df.shape[0]):
    plt.text(
        sign_snps_analysis_df['bio1'].iloc[i],
        sign_snps_analysis_df['significant_snps_count'].iloc[i],
        sign_snps_analysis_df['site'].iloc[i],
        fontsize=9, ha='right'
    )
plt.xlabel('Bio1 (Annual Mean Temperature)')
plt.ylabel('Number of Significant SNPs')
plt.title('Significant SNPs vs. Bio1')
plt.tight_layout()
plt.show()

# Plot Significant SNPs vs. Parallelism (SNP source)
plt.figure(figsize=(8, 6))
sns.scatterplot(data=sign_snps_analysis_df, x='parallelism_snp', y='significant_snps_count')
for i in range(sign_snps_analysis_df.shape[0]):
    plt.text(
        sign_snps_analysis_df['parallelism_snp'].iloc[i],
        sign_snps_analysis_df['significant_snps_count'].iloc[i],
        sign_snps_analysis_df['site'].iloc[i],
        fontsize=9, ha='right'
    )
plt.xlabel('Parallelism (SNP source)')
plt.ylabel('Number of Significant SNPs')
plt.title('Significant SNPs vs. Parallelism')
plt.tight_layout()
plt.show()

# Plot Significant SNPs vs. Unique Ecotype Count
plt.figure(figsize=(8, 6))
sns.scatterplot(data=sign_snps_analysis_df, x='unique_ecotypes', y='significant_snps_count')
for i in range(sign_snps_analysis_df.shape[0]):
    plt.text(
        sign_snps_analysis_df['unique_ecotypes'].iloc[i],
        sign_snps_analysis_df['significant_snps_count'].iloc[i],
        sign_snps_analysis_df['site'].iloc[i],
        fontsize=9, ha='right'
    )
plt.xlabel('Unique Ecotype Count')
plt.ylabel('Number of Significant SNPs')
plt.title('Significant SNPs vs. Unique Ecotype Count')
plt.tight_layout()
plt.show()

# Plot Significant SNPs vs. Shannon Diversity
plt.figure(figsize=(8, 6))
sns.scatterplot(data=sign_snps_analysis_df, x='shannon_diversity', y='significant_snps_count')
for i in range(sign_snps_analysis_df.shape[0]):
    plt.text(
        sign_snps_analysis_df['shannon_diversity'].iloc[i],
        sign_snps_analysis_df['significant_snps_count'].iloc[i],
        sign_snps_analysis_df['site'].iloc[i],
        fontsize=9, ha='right'
    )
plt.xlabel('Shannon Diversity')
plt.ylabel('Number of Significant SNPs')
plt.title('Significant SNPs vs. Shannon Diversity')
plt.tight_layout()
plt.show()

# --- Analyze Shared Significant Loci Between Sites ---
print("Analyzing shared significant loci between sites...")

# Find all pairwise combinations of sites for shared SNP analysis
pairwise_intersections_snps = {}
# Use keys from sign_snps_names_dict as they are site numbers (int)
for (site1, site2) in combinations(sign_snps_names_dict.keys(), 2):
    # Find the intersection of significant SNP names (rs IDs) between sites
    intersection = set(sign_snps_names_dict[site1]).intersection(set(sign_snps_names_dict[site2]))
    # Store the intersection (as list) with the site pair as key
    pairwise_intersections_snps[(site1, site2)] = list(intersection)

# Filter to keep only pairs with shared significant SNPs
filtered_intersections_snps = {pair: snps for pair, snps in pairwise_intersections_snps.items() if snps}

# Find all pairwise combinations of sites for shared block analysis
pairwise_intersections_blocks = {}
# Use keys from sign_blocks_names_dict as they are site numbers (int)
for (site1, site2) in combinations(sign_blocks_names_dict.keys(), 2):
    # Find the intersection of significant block names between sites
    intersection = set(sign_blocks_names_dict[site1]).intersection(set(sign_blocks_names_dict[site2]))
    # Store the intersection (as list) with the site pair as key
    pairwise_intersections_blocks[(site1, site2)] = list(intersection)

# Filter to keep only pairs with shared significant blocks
filtered_intersections_blocks = {pair: blocks for pair, blocks in pairwise_intersections_blocks.items() if blocks}

# Calculate climatic distances between all pairwise combinations of sites
climatic_distances = []
# Use site numbers from the keys of the intersection dictionaries
sites_for_distance = sorted(list(sign_snps_names_dict.keys())) # Get all site numbers
climate_filtered = climate[climate['site'].isin(sites_for_distance)] # Filter climate data for these sites

for (site1, site2) in combinations(climate_filtered['site'], 2):
    bio1_site1 = climate_filtered.loc[climate_filtered['site'] == site1, 'bio1'].values[0]
    bio1_site2 = climate_filtered.loc[climate_filtered['site'] == site2, 'bio1'].values[0]
    distance = (bio1_site1 - bio1_site2) ** 2
    climatic_distances.append({'site_pair': (site1, site2), 'climatic_distance': distance})

climatic_distances_df = pd.DataFrame(climatic_distances)

# Add columns for number and names of shared SNPs to the climatic_distances_df
climatic_distances_df['shared_snps_count'] = climatic_distances_df['site_pair'].apply(
    lambda pair: len(pairwise_intersections_snps.get(pair, [])) # Use pairwise_intersections_snps before filtering
)
climatic_distances_df['shared_snps_names'] = climatic_distances_df['site_pair'].apply(
    lambda pair: pairwise_intersections_snps.get(pair, []) # Use pairwise_intersections_snps before filtering
)

# Add columns for number and names of shared blocks to the climatic_distances_df
climatic_distances_df['shared_blocks_count'] = climatic_distances_df['site_pair'].apply(
    lambda pair: len(pairwise_intersections_blocks.get(pair, [])) # Use pairwise_intersections_blocks before filtering
)
climatic_distances_df['shared_blocks_names'] = climatic_distances_df['site_pair'].apply(
    lambda pair: pairwise_intersections_blocks.get(pair, []) # Use pairwise_intersections_blocks before filtering
)

# Display pairs of sites with shared significant SNPs
print("Pairs of sites with shared significant SNPs (count > 0):")
print(climatic_distances_df[climatic_distances_df['shared_snps_count'] > 0][['site_pair', 'climatic_distance', 'shared_snps_count']])

# Display pairs of sites with shared significant blocks
print("Pairs of sites with shared significant blocks (count > 0):")
print(climatic_distances_df[climatic_distances_df['shared_blocks_count'] > 0][['site_pair', 'climatic_distance', 'shared_blocks_count']])

# Plot shared blocks vs. climatic distance
plt.figure(figsize=(8, 6))
sns.scatterplot(data=climatic_distances_df, x = 'climatic_distance', y='shared_blocks_count')
plt.xlabel('Climatic Distance (Squared Difference in Bio1)')
plt.ylabel('Number of Shared Significant Blocks')
plt.title('Shared Significant Blocks vs. Climatic Distance')
plt.tight_layout()
plt.show()

# Count the occurrences of each block being significant across site pairs
print("Counting occurrences of each block being significant across site pairs...")
block_count = {}
# Iterate through the lists of shared block names for each site pair
for blocks in pairwise_intersections_blocks.values(): # Use pairwise_intersections_blocks before filtering
    for block in blocks:
        if block in block_count:
            block_count[block] += 1
        else:
            block_count[block] = 1

# Convert the count dictionary into a DataFrame and sort by frequency
block_count_df = pd.DataFrame(list(block_count.items()), columns=['block', 'count']).sort_values(by='count', ascending=False)

print("Occurrences of each block being significant across site pairs:")
print(block_count_df)

print("Script execution completed.")
