#!/usr/bin/env python
# coding: utf-8

# This script performs enrichment tests for antagonistic pleiotropy by analyzing SNP selection patterns
# across different experimental sites. It tests two main hypotheses:
# 1. Pleiotropy vs Conditional Neutrality: Whether SNPs show selection in multiple sites
# 2. Synergistic vs Antagonistic Pleiotropy: Whether selected SNPs show consistent or opposing effects
# across sites.

# --- Required Input Files ---
# Ensure these files are present relative to the script's location:
# - ../key_files/snp_origin_bio1_1001gvcf.csv
# - ../key_files/bioclimvars_experimental_sites_era5.csv
# - ../key_files/greneNet_final_v1.1_LDpruned.recode.vcf
# - ../key_files/climatic_distance_sites.csv
# - Files in the 'binom_reg/' directory with the pattern 'site_*.csv'
# The script will write output files to the current directory.

# --- Import Libraries ---

import pandas as pd
import numpy as np
import os
from scipy.stats import fisher_exact
from itertools import combinations
import allel
import seaborn as sns
import matplotlib.pyplot as plt

# --- Data Loading and Preparation ---

# Read SNP origin information
snp_origin_bio1 = pd.read_csv('../key_files/snp_origin_bio1_1001gvcf.csv')
snp_origin_bio1['id'] = snp_origin_bio1['chrom'].astype(str) + '_' + snp_origin_bio1['pos'].astype(str)

# Read LD-pruned SNPs
snp_pruned = allel.read_vcf('../key_files/greneNet_final_v1.1_LDpruned.recode.vcf')
snp_pruned = snp_pruned['variants/ID']

# Read and prepare site information
sign_results = [i for i in os.listdir('binom_reg/') if 'sign' in i]
sites = pd.Series(sign_results).str.split('_').str[-1].str.replace('.csv', '').astype(int).reset_index()
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

# --- Perform Enrichment Tests ---

# Dictionary to store results
results = {}

# Loop through all unique pairwise combinations of sites
for site_value_1, site_value_2 in combinations(unique_sites, 2):
    print(f"Analyzing sites {site_value_1} and {site_value_2}")
    
    # Load and prepare data for first site
    site1 = pd.read_csv(f'binom_reg/site_{site_value_1}.csv')
    site1 = site1[site1['id'].isin(snp_pruned)].reset_index(drop=True)
    th1 = 0.05 / len(site1)  # Bonferroni correction
    site1['sign'] = site1['p_value'] < th1
    site1.columns = [f'slope_site{site_value_1}', f'p_value{site_value_1}', 'id', f'sign_{site_value_1}']
    
    # Load and prepare data for second site
    site2 = pd.read_csv(f'binom_reg/site_{site_value_2}.csv')
    site2 = site2[site2['id'].isin(snp_pruned)].reset_index(drop=True)
    th2 = 0.05 / len(site2)  # Bonferroni correction
    site2['sign'] = site2['p_value'] < th2
    site2.columns = [f'slope_site{site_value_2}', f'p_value{site_value_2}', 'id', f'sign_{site_value_2}']
    
    # Merge datasets
    df = site1.merge(site2, on='id')
    
    # Test 1: Pleiotropy vs Conditional Neutrality
    table1 = pd.crosstab(df[f'sign_{site_value_1}'], df[f'sign_{site_value_2}'])
    table1 = table1.replace(0, 0.1)  # Add small value to avoid division by zero
    odds_ratio1, p_value1 = fisher_exact(table1)
    
    # Test 2: Synergistic vs Antagonistic Pleiotropy
    df_both_selected = df[(df[f'sign_{site_value_1}'] == True) & (df[f'sign_{site_value_2}'] == True)]
    
    if df_both_selected.empty:
        results[(site_value_1, site_value_2)] = {
            'odds_ratio_PvC': odds_ratio1,
            'P_value_PvC': p_value1,
            'odds_ratio_SvA': np.nan,
            'P_value_SvA': np.nan
        }
    else:
        # Classify direction of selection
        df_both_selected['direction_site1'] = np.where(df_both_selected[f'slope_site{site_value_1}'] > 0, 'Positive', 'Negative')
        df_both_selected['direction_site2'] = np.where(df_both_selected[f'slope_site{site_value_2}'] > 0, 'Positive', 'Negative')
        
        # Create contingency table for direction of selection
        table2 = pd.crosstab(df_both_selected['direction_site1'], df_both_selected['direction_site2'], dropna=False)
        categories = ['Negative', 'Positive']
        table2 = table2.reindex(index=categories, columns=categories, fill_value=0)
        
        odds_ratio2, p_value2 = fisher_exact(table2)
        results[(site_value_1, site_value_2)] = {
            'odds_ratio_PvC': odds_ratio1,
            'P_value_PvC': p_value1,
            'odds_ratio_SvA': odds_ratio2,
            'P_value_SvA': p_value2
        }

# --- Process and Save Results ---

# Convert results to DataFrame
results_odds = pd.DataFrame(results).T.reset_index()
results_odds.columns = ['site1', 'site2', 'odds_PvN', 'P_value_PvN', 'odds_SvA', 'P_value_SvA']

# Add climatic information
results_odds = results_odds.merge(clim_sites_during_exp[['site', 'bio1']], left_on='site1', right_on='site', how='left')
results_odds.rename(columns={'bio1': 'bio1_site1'}, inplace=True)
results_odds.drop(columns=['site'], inplace=True)

results_odds = results_odds.merge(clim_sites_during_exp[['site', 'bio1']], left_on='site2', right_on='site', how='left')
results_odds.rename(columns={'bio1': 'bio1_site2'}, inplace=True)
results_odds.drop(columns=['site'], inplace=True)

# Calculate climatic distance
results_odds['climatic_distance'] = (results_odds['bio1_site1'] - results_odds['bio1_site2']) ** 2

# Add geographic distance information
melted_distance_df = pd.read_csv('../key_files/climatic_distance_sites.csv')
results_odds = results_odds.merge(melted_distance_df, on=['site1', 'site2'])

# Save results
results_odds.to_csv('odds_ap_cn_ldprun_fisher.csv', index=None)

# --- Visualization ---

# Create scatter plot of odds ratios vs climatic distance
plt.figure(figsize=(10, 6))
sns.scatterplot(data=results_odds, x='climatic_distance', y='odds_PvN')
plt.title('Relationship between Climatic Distance and Pleiotropy')
plt.xlabel('Squared Climatic Distance')
plt.ylabel('Odds Ratio (Pleiotropy vs Neutrality)')
plt.savefig('pleiotropy_vs_climatic_distance.pdf')
plt.close()
