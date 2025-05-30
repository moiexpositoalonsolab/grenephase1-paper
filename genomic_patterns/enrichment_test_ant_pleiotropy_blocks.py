#!/usr/bin/env python
# coding: utf-8

# This script performs enrichment tests for antagonistic pleiotropy by analyzing block-level WZA results
# across different experimental sites. It tests two main hypotheses:
# 1. Pleiotropy vs Conditional Neutrality: Whether blocks show selection in multiple sites
# 2. Synergistic vs Antagonistic Pleiotropy: Whether selected blocks show consistent or opposing effects
# across sites.
# The analysis is performed using Chi-squared tests, and results are visualized against climatic and geographic distance.

# --- Required Input Files ---
# Ensure these files are present relative to the script's location:
# - ../key_files/snp_origin_bio1_1001gvcf.csv
# - ../key_files/bioclimvars_experimental_sites_era5.csv
# - ../key_files/greneNet_final_v1.1_LDpruned.recode.vcf
# - ../key_files/climatic_distance_sites.csv
# - Files in the 'binom_reg/' directory with the pattern 'wza_site_*_pr.csv' (WZA results for blocks)
# The script will write output files to the current directory, including:
# - odds_ap_cn_blocks.csv
# - chi_square_pvcn_vs_climatic_distance_blocks.pdf
# - chi_square_sva_vs_climatic_distance_blocks.pdf

# --- Import Libraries ---

import pandas as pd
import numpy as np
import os
from scipy.stats import chi2_contingency # Using chi2_contingency for enrichment tests
from itertools import combinations
import allel
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform

# --- Data Loading and Preparation ---

# Read SNP origin information (used for general SNP info, though not directly for block tests)
snp_origin_bio1 = pd.read_csv('../key_files/snp_origin_bio1_1001gvcf.csv')
snp_origin_bio1['id'] = snp_origin_bio1['chrom'].astype(str) + '_' + snp_origin_bio1['pos'].astype(str)

# Read LD-pruned SNPs (used for filtering, though applied at block level here)
snp_pruned = allel.read_vcf('../key_files/greneNet_final_v1.1_LDpruned.recode.vcf')
snp_pruned = snp_pruned['variants/ID'] # This might need adjustment if block IDs are used instead of SNP IDs

# Read and prepare site information
sign_results = [i for i in os.listdir('binom_reg/') if 'wza_site' in i and '_pr.csv' in i]
sites = pd.Series(sign_results).str.split('_').str[-2].astype(int).reset_index()
sites.columns = ['ign', 'site']

# Read climate data for experimental sites
clim_sites_during_exp = pd.read_csv('../key_files/bioclimvars_experimental_sites_era5.csv')
clim_sites_during_exp = clim_sites_during_exp[['site', 'bio1']]
clim_sites_during_exp = clim_sites_during_exp.merge(sites['site'])
clim_sites_during_exp_ordered = clim_sites_during_exp.sort_values('bio1')
unique_sites = clim_sites_during_exp_ordered['site'].values

# Classify sites as hot or cold based on bio1 temperature
clim_sites_during_exp['siteis'] = 'hot'
clim_sites_during_exp.loc[clim_sites_during_exp['bio1'] < 12, 'siteis'] = 'cold'

# Calculate pairwise climatic distances between sites
bio_columns = [col for col in clim_sites_during_exp if col.startswith('bio')]
bio_data = clim_sites_during_exp[bio_columns]
distance_matrix = squareform(pdist(bio_data, metric='euclidean'))
distance_df = pd.DataFrame(distance_matrix, index=clim_sites_during_exp['site'], columns=clim_sites_during_exp['site'])
melted_distance_df = distance_df.reset_index().rename(columns={'site': 'site1'}).melt(id_vars=['site1'], var_name='site2', value_name='climatic_distance')
melted_distance_df = melted_distance_df[melted_distance_df['site1'] != melted_distance_df['site2']]

# Read geographic distance information (assuming this file contains pre-calculated distances)
# Note: The original notebook seemed to regenerate and save this, but listing it as an input
# here for clarity based on typical workflow.
geo_distance_df = pd.read_csv('../key_files/climatic_distance_sites.csv')


# --- Perform Enrichment Tests on Blocks ---

# Dictionary to store results
results = {}

# Loop through all unique pairwise combinations of sites
for site_value_1, site_value_2 in combinations(unique_sites, 2):
    print(f"Analyzing sites {site_value_1} and {site_value_2}")
    
    # Load and prepare WZA data for first site (blocks)
    try:
        site1_wza = pd.read_csv(f'binom_reg/wza_site_{site_value_1}_pr.csv')
        # Ensure 'gene' column exists and is used as the block ID
        if 'gene' not in site1_wza.columns:
             print(f"Warning: 'gene' column not found in wza_site_{site_value_1}_pr.csv. Skipping site pair.")
             continue
        site1_wza.rename(columns={'gene': 'block_id'}, inplace=True)
        th1 = 0.05 / len(site1_wza)  # Bonferroni correction based on number of blocks
        site1_wza['sign'] = site1_wza['Z_pVal'] < th1
        site1_processed = site1_wza[['slope', 'Z_pVal', 'block_id', 'sign']].copy()
        site1_processed.columns = [f'slope_site{site_value_1}', f'p_value{site_value_1}', 'block_id', f'sign_{site_value_1}']
    except FileNotFoundError:
        print(f"Warning: File binom_reg/wza_site_{site_value_1}_pr.csv not found. Skipping site pair.")
        continue

    # Load and prepare WZA data for second site (blocks)
    try:
        site2_wza = pd.read_csv(f'binom_reg/wza_site_{site_value_2}_pr.csv')
        # Ensure 'gene' column exists and is used as the block ID
        if 'gene' not in site2_wza.columns:
             print(f"Warning: 'gene' column not found in wza_site_{site_value_2}_pr.csv. Skipping site pair.")
             continue
        site2_wza.rename(columns={'gene': 'block_id'}, inplace=True)
        th2 = 0.05 / len(site2_wza)  # Bonferroni correction based on number of blocks
        site2_wza['sign'] = site2_wza['Z_pVal'] < th2
        site2_processed = site2_wza[['slope', 'Z_pVal', 'block_id', 'sign']].copy()
        site2_processed.columns = [f'slope_site{site_value_2}', f'p_value{site_value_2}', 'block_id', f'sign_{site_value_2}']
    except FileNotFoundError:
        print(f"Warning: File binom_reg/wza_site_{site_value_2}_pr.csv not found. Skipping site pair.")
        continue
    
    # Merge datasets on 'block_id'
    df = site1_processed.merge(site2_processed, on='block_id', how='inner')
    
    # Skip if no common blocks after merge
    if df.empty:
        print(f"No common blocks between sites {site_value_1} and {site_value_2}. Skipping.")
        continue

    # Test 1: Pleiotropy vs Conditional Neutrality (using Chi-squared)
    # Contingency table: [[sign_site1=False, sign_site2=False], [sign_site1=False, sign_site2=True], [sign_site1=True, sign_site2=False], [sign_site1=True, sign_site2=True]]
    table1 = pd.crosstab(df[f'sign_{site_value_1}'], df[f'sign_{site_value_2}'])
    # Ensure table has both True/False for both sites, fill with 0 if not present
    table1 = table1.reindex(index=[False, True], columns=[False, True], fill_value=0)
    
    # Perform Chi-squared test
    chi1, p_value1, _, _ = chi2_contingency(table1)
    
    # Test 2: Synergistic vs Antagonistic Pleiotropy (using Chi-squared on selected blocks)
    df_both_selected = df[(df[f'sign_{site_value_1}'] == True) & (df[f'sign_{site_value_2}'] == True)].copy()
    
    if df_both_selected.empty:
        results[(site_value_1, site_value_2)] = {
            'Chi_square_PvC': chi1,
            'P_value_PvC': p_value1,
            'Chi_square_SvA': np.nan,
            'P_value_SvA': np.nan
        }
    else:
        # Classify direction of selection based on slope
        df_both_selected['direction_site1'] = np.where(df_both_selected[f'slope_site{site_value_1}'] > 0, 'Positive', 'Negative')
        df_both_selected['direction_site2'] = np.where(df_both_selected[f'slope_site{site_value_2}'] > 0, 'Positive', 'Negative')
        
        # Create contingency table for direction of selection
        # Table: [[direction_site1=Negative, direction_site2=Negative], [direction_site1=Negative, direction_site2=Positive], [direction_site1=Positive, direction_site2=Negative], [direction_site1=Positive, direction_site2=Positive]]
        table2 = pd.crosstab(df_both_selected['direction_site1'], df_both_selected['direction_site2'], dropna=False)
        # Ensure table has both Positive/Negative for both sites, fill with 0 if not present
        categories = ['Negative', 'Positive']
        table2 = table2.reindex(index=categories, columns=categories, fill_value=0)
        
        # Perform Chi-squared test
        chi2, p_value2, _, _ = chi2_contingency(table2)
        
        results[(site_value_1, site_value_2)] = {
            'Chi_square_PvC': chi1,
            'P_value_PvC': p_value1,
            'Chi_square_SvA': chi2,
            'P_value_SvA': p_value2
        }

# --- Process and Save Results ---

# Convert results to DataFrame
results_chi = pd.DataFrame(results).T.reset_index()
results_chi.columns = ['site1', 'site2', 'Chi_square_PvN', 'P_value_PvN', 'Chi_square_SvA', 'P_value_SvA']

# Merge with climatic and geographic distance information
# Merge with climatic distance
results_chi = results_chi.merge(melted_distance_df, on=['site1', 'site2'], how='left')

# Merge with geographic distance
# Assuming 'geo_distance_df' has columns 'site1', 'site2', and 'distance'
results_chi = results_chi.merge(geo_distance_df[['site1', 'site2', 'distance']], on=['site1', 'site2'], how='left')

# Save the results to a CSV file
results_chi.to_csv('odds_ap_cn_blocks.csv', index=None)


# --- Visualization ---

# Scatter plot of Chi-square for Pleiotropy vs Conditional Neutrality vs Climatic Distance
plt.figure(figsize=(10, 6))
sns.scatterplot(data=results_chi, x='climatic_distance', y='Chi_square_PvN')
plt.title('Chi-square (Pleiotropy vs Conditional Neutrality) vs Climatic Distance (Blocks)')
plt.xlabel('Euclidean Climatic Distance')
plt.ylabel('Chi-square Statistic')
plt.savefig('chi_square_pvcn_vs_climatic_distance_blocks.pdf')
plt.close()

# Scatter plot of Chi-square for Synergistic vs Antagonistic Pleiotropy vs Climatic Distance
plt.figure(figsize=(10, 6))
sns.scatterplot(data=results_chi.dropna(subset=['Chi_square_SvA']), x='climatic_distance', y='Chi_square_SvA')
plt.title('Chi-square (Synergistic vs Antagonistic Pleiotropy) vs Climatic Distance (Blocks)')
plt.xlabel('Euclidean Climatic Distance')
plt.ylabel('Chi-square Statistic')
plt.savefig('chi_square_sva_vs_climatic_distance_blocks.pdf')
plt.close()

# Scatter plot of Chi-square for Pleiotropy vs Conditional Neutrality vs Geographic Distance
plt.figure(figsize=(10, 6))
sns.scatterplot(data=results_chi, x='distance', y='Chi_square_PvN')
plt.title('Chi-square (Pleiotropy vs Conditional Neutrality) vs Geographic Distance (Blocks)')
plt.xlabel('Geographic Distance')
plt.ylabel('Chi-square Statistic')
plt.savefig('chi_square_pvcn_vs_geographic_distance_blocks.pdf')
plt.close()

# Scatter plot of Chi-square for Synergistic vs Antagonistic Pleiotropy vs Geographic Distance
plt.figure(figsize=(10, 6))
sns.scatterplot(data=results_chi.dropna(subset=['Chi_square_SvA']), x='distance', y='Chi_square_SvA')
plt.title('Chi-square (Synergistic vs Antagonistic Pleiotropy) vs Geographic Distance (Blocks)')
plt.xlabel('Geographic Distance')
plt.ylabel('Chi-square Statistic')
plt.savefig('chi_square_sva_vs_geographic_distance_blocks.pdf')
plt.close()

# Note: The original script also included analysis based on 'hot' and 'cold' site classifications
# and printed contingency tables. This streamlined version focuses on the primary outputs (CSV and plots).
# The 'hot'/'cold' classification logic is kept in the data loading section for potential future use.
