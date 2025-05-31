#!/usr/bin/env python
# coding: utf-8

"""
Launch Weighted Z-score Analysis (WZA) Jobs

This script prepares and submits SBATCH job scripts to run the Weighted Z-score Analysis (WZA)
script (`general_WZA_script_mod_polynomial_order7.py`) for the processed GEMMA results
from each experimental site.

It iterates through site-specific directories containing GEMMA results, generates a unique
SBATCH script for each site with the necessary command-line arguments for the WZA script,
and submits these scripts to a Slurm cluster scheduler using `sbatch`.

Inputs:
- Processed GEMMA results files: `/carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/fitness_gwas/site_*/output/results_lmm.csv`
- The WZA script: `general_WZA_script_mod_polynomial_order7.py` (assumed to be in the same directory or accessible in the PATH)
- ../../data/bioclimvars_experimental_sites_era5.csv: Climate data for experimental sites (used for downstream analysis/plotting, not for launching WZA).
- ../../data/blocks_snpsid_dict.pkl: Pickle file containing a dictionary mapping SNP IDs to block IDs (used for downstream analysis/plotting, not for launching WZA).

Outputs:
- For each site: `/carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/fitness_gwas/wza_site_<site_id>.sh`: SBATCH script for running WZA.
- Submitted WZA jobs on the cluster via sbatch.
- (Potential) Interactive plots visualizing WZA results and shared significant blocks vs. climatic distance.
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
from itertools import combinations
import pickle

# --- Configuration and File Paths ---
# Path to the base directory containing site-specific GEMMA results
gemma_results_base_path = '.'

# Name of the WZA script
wza_script_name = 'general_WZA_script_mod_polynomial_order7.py'

# Path to climate data (for downstream analysis/plotting)
climate_file = '../../data/bioclimvars_experimental_sites_era5.csv'

# Path to blocks dictionary (for downstream analysis/plotting)
blocks_dict_file = '../../data/blocks_snpsid_dict.pkl' # Using relative path based on common file location

# --- Helper Function to List Site Directories ---
def list_only_directories(path):
    """
    Lists names of directories within a given path that contain 'site'.
    """
    all_entries = os.listdir(path)
    site_directories = [entry for entry in all_entries if 'site' in entry and os.path.isdir(os.path.join(path, entry))]
    return site_directories

# Get list of site directories
site_directories = list_only_directories(gemma_results_base_path)

# --- Create and Submit WZA Job Scripts per Site ---
print("Creating and submitting WZA job scripts...")
shfiles = []

# Iterate through each site directory to create and submit WZA job scripts
for site in site_directories:
    site_number = site.replace('site_', '')
    # Construct the command to run the WZA script
    # Ensure paths within the command are correct for the cluster environment
    wza_command = f'python {wza_script_name} --correlations {site}/output/results_lmm.csv --summary_stat p_wald --window blocks --output {site}/wza_fitness_gwas_poly7.csv --sep ","'
    
    # Define SBATCH script content
    # Note: Adjust resource requests (--time, --mem-per-cpu) and module/conda activation as needed for the cluster.
    file_name = f'wza_{site}.sh'
    file_path = os.path.join(gemma_results_base_path, file_name) # Save script in the base fitness_gwas directory
    
    text = f'''#!/bin/bash
#SBATCH --job-name=wza_{site_number}
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --output=wza_{site_number}_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

# Load necessary modules or activate conda environment
source /home/tbellagio/miniforge3/etc/profile.d/conda.sh
conda activate /home/tbellagio/miniforge3/envs/pipeline_snakemake

# Set LD_LIBRARY_PATH if required by the WZA script or its dependencies
export LD_LIBRARY_PATH="/home/tbellagio/miniforge3/envs/run_baypass/lib:$LD_LIBRARY_PATH"

# Change to the fitness_gwas directory where input files and WZA script are located
cd {gemma_results_base_path}

# Execute the WZA command
{wza_command}

'''
    
    # Write the SBATCH script to file
    with open(file_path, 'w') as o:
        o.write(text)
    shfiles.append(file_path)
    print(f"Created SBATCH script: {file_path}")

# Submit SBATCH jobs
print("Submitting SBATCH jobs...")
for i, file in enumerate(shfiles):
    try:
        # Use check=True to raise an exception if the command fails
        subprocess.run(['sbatch', file], check=True)
        print(f"Submitted job {i+1}/{len(shfiles)}: {file}")
    except subprocess.CalledProcessError as e:
        print(f"Error submitting job {i+1}/{len(shfiles)} {file}: {e}")
    except FileNotFoundError:
        print(f"Error: sbatch command not found. Make sure Slurm is in your PATH.")

print("SBATCH job submission process completed.")

# --- Downstream Analysis and Visualization (Optional) ---
# The following code appears to perform analysis and visualization of WZA results.
# It is kept here for completeness but is commented out as it might be run separately
# or depends on the WZA jobs completing successfully.

# print("Performing downstream analysis and visualization of WZA results...")

# # Load climate data (for shared block analysis vs. climatic distance)
# try:
#     climate = pd.read_csv(climate_file)[['site', 'bio1']]
#     climate['site'] = climate['site'].astype(int)
# except FileNotFoundError:
#     print(f"Warning: Climate file not found at {climate_file}. Skipping shared block analysis vs. climatic distance.")
#     climate = pd.DataFrame()

# # Load blocks dictionary (for mapping block IDs)
# try:
#     with open(blocks_dict_file, 'rb') as f:
#         dict_blocks = pickle.load(f)
#     reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}
# except FileNotFoundError:
#     print(f"Warning: Blocks dictionary file not found at {blocks_dict_file}. Skipping some downstream analysis.")
#     dict_blocks = {}
#     reverse_mapping = {}

# # Aggregate significant WZA blocks per site
# sign_blocks_names_dict_wza = {}

# for site in site_directories:
#     site_number = int(site.replace('site_', ''))
#     wza_results_file = os.path.join(gemma_results_base_path, site, 'wza_fitness_gwas_poly7.csv')
    
#     if not os.path.exists(wza_results_file):
#         print(f"Warning: WZA results file not found for site {site}: {wza_results_file}. Skipping.")
#         continue
        
#     try:
#         wza_df = pd.read_csv(wza_results_file)
        
#         # Assuming significance is based on 'Z_pVal' and a Bonferroni correction
#         # Need to determine the correct denominator for the Bonferroni correction here.
#         # A common approach for WZA is correction by the number of windows/blocks tested.
#         # Let's assume len(wza_df) is the number of blocks tested for this site.
#         wza_threshold = 0.05 / len(wza_df)
        
#         significant_blocks = wza_df[wza_df['Z_pVal'] < wza_threshold]['gene'].unique().tolist()
#         sign_blocks_names_dict_wza[site_number] = significant_blocks
#     except pd.errors.EmptyDataError:
#          print(f"Warning: WZA results file is empty for site {site}: {wza_results_file}. Skipping.")
#     except KeyError as e:
#          print(f"Warning: Expected column {e} not found in WZA results for site {site}: {wza_results_file}. Skipping.")

# # Analysis of shared significant blocks between sites based on WZA results
# if sign_blocks_names_dict_wza and not climate.empty:
#     print("Analyzing shared significant WZA blocks between sites...")
    
#     # Find all pairwise combinations of sites
#     pairwise_intersections_wza_blocks = {}
#     for (site1, site2) in combinations(sign_blocks_names_dict_wza.keys(), 2):
#         intersection = set(sign_blocks_names_dict_wza[site1]).intersection(set(sign_blocks_names_dict_wza[site2]))
#         pairwise_intersections_wza_blocks[(site1, site2)] = list(intersection)

#     # Filter to keep only pairs with shared significant WZA blocks
#     filtered_intersections_wza_blocks = {pair: blocks for pair, blocks in pairwise_intersections_wza_blocks.items() if blocks}

#     # Calculate climatic distances between all pairwise combinations of sites with WZA results
#     climatic_distances_wza_sites = []
#     sites_for_distance_wza = sorted(list(sign_blocks_names_dict_wza.keys()))
#     climate_filtered_wza = climate[climate['site'].isin(sites_for_distance_wza)]

#     for (site1, site2) in combinations(climate_filtered_wza['site'], 2):
#         bio1_site1 = climate_filtered_wza.loc[climate_filtered_wza['site'] == site1, 'bio1'].values[0]
#         bio1_site2 = climate_filtered_wza.loc[climate_filtered_wza['site'] == site2, 'bio1'].values[0]
#         distance = (bio1_site1 - bio1_site2) ** 2
#         climatic_distances_wza_sites.append({'site_pair': (site1, site2), 'climatic_distance': distance})

#     climatic_distances_wza_df = pd.DataFrame(climatic_distances_wza_sites)

#     # Add columns for number and names of shared blocks to the climatic_distances_wza_df
#     climatic_distances_wza_df['shared_blocks_count_wza'] = climatic_distances_wza_df['site_pair'].apply(
#         lambda pair: len(pairwise_intersections_wza_blocks.get(pair, []))
#     )
#     climatic_distances_wza_df['shared_blocks_names_wza'] = climatic_distances_wza_df['site_pair'].apply(
#         lambda pair: pairwise_intersections_wza_blocks.get(pair, [])
#     )

#     # Display pairs of sites with shared significant WZA blocks
#     print("Pairs of sites with shared significant WZA blocks (count > 0):")
#     print(climatic_distances_wza_df[climatic_distances_wza_df['shared_blocks_count_wza'] > 0][['site_pair', 'climatic_distance', 'shared_blocks_count_wza']])

#     # Plot shared WZA blocks vs. climatic distance
#     plt.figure(figsize=(8, 6))
#     sns.scatterplot(data=climatic_distances_wza_df, x = 'climatic_distance', y='shared_blocks_count_wza')
#     plt.xlabel('Climatic Distance (Squared Difference in Bio1)')
#     plt.ylabel('Number of Shared Significant WZA Blocks')
#     plt.title('Shared Significant WZA Blocks vs. Climatic Distance')
#     plt.tight_layout()
#     plt.show()

#     # Count the occurrences of each block being significant across site pairs (based on WZA)
#     print("Counting occurrences of each block being significant across site pairs (WZA)...")
#     block_count_wza = {}
#     for blocks in pairwise_intersections_wza_blocks.values():
#         for block in blocks:
#             if block in block_count_wza:
#                 block_count_wza[block] += 1
#             else:
#                 block_count_wza[block] = 1

#     block_count_wza_df = pd.DataFrame(list(block_count_wza.items()), columns=['block', 'count']).sort_values(by='count', ascending=False)

#     print("Occurrences of each block being significant across site pairs (WZA):")
#     print(block_count_wza_df)

# # WZA Results Visualization (Manhattan-like plots)
# # This section plots the Z-score p-values from WZA results.
# print("Generating WZA results plots...")

# for site in site_directories:
#     site_number = site.replace('site_', '')
#     wza_results_file = os.path.join(gemma_results_base_path, site, 'wza_fitness_gwas_poly7.csv')
    
#     if not os.path.exists(wza_results_file):
#         print(f"Warning: WZA results file not found for site {site}: {wza_results_file}. Skipping plot generation.")
#         continue
        
#     try:
#         wza_df = pd.read_csv(wza_results_file)

#         if wza_df.empty:
#              print(f"Warning: WZA results file is empty for site {site}: {wza_results_file}. Skipping plot generation.")
#              continue
        
#         # Prepare data for plotting
#         # Assuming 'gene' column contains block IDs like 'chrom_pos'
#         wza_df['chrom'] = wza_df['gene'].str.split('_').str[0].astype(int)
#         wza_df['pos'] = wza_df['gene'].str.split('_').str[1].astype(int)
        
#         # Determine significance threshold for plotting (using WZA Z_pVal)
#         wza_threshold_plot = 0.05 / len(wza_df)
        
#         plot_df = wza_df[['Z_pVal', 'pos', 'chrom']].copy()
#         plot_df.columns = ['pvalue', 'position', 'chromosome']
#         plot_df['-log10(pvalue)'] = -np.log10(plot_df['pvalue'])
        
#         colors = sns.color_palette("crest", n_colors=len(plot_df['chromosome'].unique())) # Adjust colors based on number of chromosomes
        
#         # Calculate chromosome offsets for adjusted position
#         chromosome_offsets = {}
#         offset = 0
#         for chrom in sorted(plot_df['chromosome'].unique()):
#             chromosome_offsets[chrom] = offset
#             max_position = plot_df[plot_df['chromosome'] == chrom]['position'].max()
#             offset += max_position + 200 # Buffer between chromosomes
        
#         plot_df['adjusted_position'] = plot_df.apply(lambda row: row['position'] + chromosome_offsets[row['chromosome']], axis=1)
        
#         # Creating the Manhattan-like plot for WZA results
#         plt.figure(figsize=(20, 6))
        
#         for i, chrom in enumerate(sorted(plot_df['chromosome'].unique())):
#             subset = plot_df[plot_df['chromosome'] == chrom]
#             plt.scatter(
#                 subset['adjusted_position'],
#                 subset['-log10(pvalue)'],
#                 alpha=0.7,
#                 c=colors[i % len(colors)],
#                 label=f'Chr {chrom}',
#                 s=20
#             )
        
#         # Aesthetics
#         plt.xlabel('Adjusted Position (Blocks)')
#         plt.ylabel('-log10(WZA Z-score p-value)')
#         plt.title(f'WZA Results Manhattan Plot - Site {site_number}')
#         plt.grid(axis='y')
        
#         # Add circles around specific genes/blocks (example)
#         # genes_to_highlight = ['2_1265', '4_801'] # Example block IDs
#         # for gene in genes_to_highlight:
#         #     if '_ ' in gene:
#         #         chrom_str, pos_str = gene.split('_')
#         #         try:
#         #             chrom = int(chrom_str)
#         #             pos = int(pos_str)
#         #         except ValueError:
#         #             print(f"Warning: Could not parse chromosome or position for gene {gene}. Skipping highlighting.")
#         #             continue
        
#         #         subset_gene = plot_df[(plot_df['chromosome'] == chrom) & (plot_df['position'] == pos)]
#         #         if not subset_gene.empty:
#         #             plt.scatter(
#         #                 subset_gene['adjusted_position'],
#         #                 subset_gene['-log10(pvalue)'],
#         #                 edgecolor='red',
#         #                 linewidth=2,
#         #                 facecolor='none',
#         #                 s=100,
#         #                 label=f'Block {gene}'
#         #             )
                
#         # Threshold line
#         threshold_plot_val = -np.log10(wza_threshold_plot)
#         plt.axhline(y=threshold_plot_val, color='grey', linestyle='dashed', label=f'Bonferroni Threshold ({wza_threshold_plot:.2e})')
        
#         # plt.legend(title="Chromosome", bbox_to_anchor=(1.05, 1), loc='upper left') # Optional: Add legend
        
#         plt.tight_layout()
#         plt.show()
        
#     except pd.errors.EmptyDataError:
#          print(f"Warning: WZA results file is empty for site {site}: {wza_results_file}. Skipping plot generation.")
#     except KeyError as e:
#          print(f"Warning: Expected column {e} not found in WZA results for site {site}: {wza_results_file}. Skipping plot generation.")

# # WZA Results Visualization (Manhattan-like plots using top_candidate_p)
# # This section plots the top candidate p-values from WZA results.
# print("Generating WZA top candidate p-value plots...")

# for site in site_directories:
#     site_number = site.replace('site_', '')
#     wza_results_file = os.path.join(gemma_results_base_path, site, 'wza_fitness_gwas_poly7.csv')
    
#     if not os.path.exists(wza_results_file):
#         print(f"Warning: WZA results file not found for site {site}: {wza_results_file}. Skipping plot generation.")
#         continue
        
#     try:
#         wza_df = pd.read_csv(wza_results_file)

#         if wza_df.empty:
#              print(f"Warning: WZA results file is empty for site {site}: {wza_results_file}. Skipping plot generation.")
#              continue

#         # Prepare data for plotting (using top_candidate_p)
#         wza_df['chrom'] = wza_df['gene'].str.split('_').str[0].astype(int)
#         wza_df['pos'] = wza_df['gene'].str.split('_').str[1].astype(int)

#         # Determine significance threshold for plotting (using top_candidate_p)
#         # Need to determine the correct denominator for the Bonferroni correction here,
#         # potentially based on the number of top candidates or total blocks.
#         # Let's use the same threshold as for Z_pVal for consistency in visualization,
#         # although the interpretation might differ.
#         wza_threshold_plot = 0.05 / len(wza_df) # Using number of blocks as denominator
        
#         plot_df = wza_df[['top_candidate_p', 'pos', 'chrom']].copy()
#         plot_df.columns = ['pvalue', 'position', 'chromosome']
#         plot_df['-log10(pvalue)'] = -np.log10(plot_df['pvalue'])

#         colors = sns.color_palette("crest", n_colors=len(plot_df['chromosome'].unique()))

#         # Calculate chromosome offsets for adjusted position (using same logic)
#         chromosome_offsets = {}
#         offset = 0
#         for chrom in sorted(plot_df['chromosome'].unique()):
#             chromosome_offsets[chrom] = offset
#             max_position = plot_df[plot_df['chromosome'] == chrom]['position'].max()
#             offset += max_position + 200

#         plot_df['adjusted_position'] = plot_df.apply(lambda row: row['position'] + chromosome_offsets[row['chromosome']], axis=1)

#         # Creating the Manhattan-like plot for WZA top candidate p-values
#         plt.figure(figsize=(20, 6))

#         for i, chrom in enumerate(sorted(plot_df['chromosome'].unique())):
#             subset = plot_df[plot_df['chromosome'] == chrom]
#             plt.scatter(
#                 subset['adjusted_position'],
#                 subset['-log10(pvalue)'],
#                 alpha=0.7,
#                 c=colors[i % len(colors)],
#                 label=f'Chr {chrom}',
#                 s=20
#             )

#         # Aesthetics
#         plt.xlabel('Adjusted Position (Blocks)')
#         plt.ylabel('-log10(WZA Top Candidate p-value)')
#         plt.title(f'WZA Top Candidate p-value Plot - Site {site_number}')
#         plt.grid(axis='y')

#         # Add circles around specific genes/blocks (example)
#         # genes_to_highlight = ['2_1265', '4_801'] # Example block IDs
#         # for gene in genes_to_highlight:
#         #     if '_ ' in gene:
#         #         chrom_str, pos_str = gene.split('_')
#         #         try:
#         #             chrom = int(chrom_str)
#         #             pos = int(pos_str)
#         #         except ValueError:
#         #             print(f"Warning: Could not parse chromosome or position for gene {gene}. Skipping highlighting.")
#         #             continue

#         #         subset_gene = plot_df[(plot_df['chromosome'] == chrom) & (plot_df['position'] == pos)]
#         #         if not subset_gene.empty:
#         #             plt.scatter(
#         #                 subset_gene['adjusted_position'],
#         #                 subset_gene['-log10(pvalue)'],
#         #                 edgecolor='red',
#         #                 linewidth=2,
#         #                 facecolor='none',
#         #                 s=100,
#         #                 label=f'Block {gene}'
#         #             )

#         # Threshold line
#         threshold_plot_val = -np.log10(wza_threshold_plot)
#         plt.axhline(y=threshold_plot_val, color='grey', linestyle='dashed', label=f'Bonferroni Threshold ({wza_threshold_plot:.2e})')

#         # plt.legend(title="Chromosome", bbox_to_anchor=(1.05, 1), loc='upper left') # Optional: Add legend

#         plt.tight_layout()
#         plt.show()

#     except pd.errors.EmptyDataError:
#          print(f"Warning: WZA results file is empty for site {site}: {wza_results_file}. Skipping plot generation.")
#     except KeyError as e:
#          print(f"Warning: Expected column {e} not found in WZA results for site {site}: {wza_results_file}. Skipping plot generation.")


# print("Downstream analysis and visualization completed.")
