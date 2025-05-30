#!/usr/bin/env python3
# coding: utf-8

"""
Allele Frequency Change vs Origin Analysis

This script analyzes the relationship between allele frequency changes and their origin
in the genome.

Required Input Files:
1. ../data/var_pos_grenenet.csv
2. ../data/blocks_snpsid_dict.pkl
3. ../data/1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv
4. ../data/founder_ecotype_names.csv

Output Files:
1. allele_freq_change_vs_origin.csv
   - Results of the analysis
2. allele_freq_change_vs_origin.pdf
   - Plot of the results
"""

# Import necessary libraries
import pandas as pd
import allel
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
from adjustText import adjust_text # For adjusting text labels in the plot
import statsmodels.api as sm # For linear regression
import pickle

# Set matplotlib font type for PDF output
matplotlib.rcParams['pdf.fonttype'] = 42

# --- Load Data ---
# Load SNP positions
snps_names = pd.read_csv('../data/var_pos_grenenet.csv')

# Load blocks dictionary
with open('../data/blocks_snpsid_dict.pkl', 'rb') as f:
    blocks_dict = pickle.load(f)

# Load climate data
clim1001 = pd.read_csv('../data/1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv')

# Load founder ecotypes
founder_ecotypes = pd.read_csv('../data/founder_ecotype_names.csv')['0'].values

# --- Process Data ---
# Create a reverse mapping from SNP ID to block identifier
reverse_mapping = {item: key for key, values in blocks_dict.items() for item in values}

# Filter SNPs with sufficient allele data in the first generation
snps_names = snps_names[snps_names['total_alleles05filter_firstgen'].notna()].reset_index(drop=True)

# Convert ecotype ID to string for consistent merging/indexing
clim1001['ecotypeid'] = clim1001['ecotypeid'].astype(str)

# Read VCF file containing genotype data
# VCF (Variant Call Format) files store genetic variation data (SNPs) for individuals.
vcfgrene = allel.read_vcf('../gwas/greneNet_final_v1.1.recode.vcf.gz')

# Extract genotype data and convert to a pandas DataFrame
# The genotype data is converted to represent the number of non-reference alleles (0, 1, or 2 for diploids).
geno = vcfgrene['calldata/GT'].sum(axis=2)
geno = pd.DataFrame(geno)
geno.columns = vcfgrene['samples'] # Set columns to sample (ecotype) names

# Filter climate data to include only ecotypes present in the genotype data
clim1001 = clim1001[clim1001['ecotypeid'].isin(geno.columns)]

# Further filter climate data to include only founder ecotypes
clim1001 = clim1001[clim1001['ecotypeid'].isin(founder_ecotypes.astype(str))]

# Select only ecotype ID and the specific climate variable 'bio1' (Annual Mean Temperature)
clim1001 = clim1001[['ecotypeid', 'bio1']]

# Convert ecotype ID to string again (addresses potential SettingWithCopyWarning)
clim1001['ecotypeid'] = clim1001['ecotypeid'].astype(str)

# Set ecotype ID as the index and sort both genotype and climate dataframes by index
# This aligns the two dataframes based on ecotype.
clim1001 = clim1001.set_index('ecotypeid').sort_index()
geno = geno.T.sort_index() # Transpose genotype data so ecotypes are rows and sort

# Convert climate data to a numpy array and center it around the mean
# Centering is a common step before calculating associations.
clim1001_np = np.array(clim1001['bio1'])
temperature_centered = clim1001_np - np.mean(clim1001_np)

# Calculate SNP-temperature association
# This is done by matrix multiplication of transposed genotype data and centered temperature.
snp_temperature_association = np.array(geno).T.dot(temperature_centered)

# Normalize the association by the number of ecotypes
snp_temperature_association /= len(clim1001_np)

# Convert the association results to a pandas Series
snp_temperature_association = pd.Series(snp_temperature_association.flatten())

# Create a DataFrame for SNP origin (temperature association)
snp_origin_bio1 = pd.DataFrame({
    'chrom': vcfgrene['variants/CHROM'],
    'pos': vcfgrene['variants/POS'],
    'snp_origin_bio1': snp_temperature_association
})

# Create a unique SNP ID (chromosome_position)
snp_origin_bio1['id'] = snp_origin_bio1['chrom'].astype(str) + '_' + snp_origin_bio1['pos'].astype(str)

# Map SNP IDs to their corresponding block identifiers
snp_origin_bio1['blocks'] = snp_origin_bio1['id'].map(reverse_mapping)

# Calculate the mean SNP-temperature association for each block
block_origin_bio1 = snp_origin_bio1.groupby('blocks')['snp_origin_bio1'].mean().reset_index()

# Read results from binomial regression (allele frequency change vs. temperature)
# This file contains the 'slope' which represents the association of temperature with allele frequency change.
binomial_reg = pd.read_csv('../binomial_regression_lastgen/binomial_reg_lastgen_wmaf_bio1.csv')

# Filter SNPs with sufficient allele data in the last generation and get their IDs
dict_snps_lastgen = pd.read_csv('../data/var_pos_grenenet.csv')  # Re-read if needed, or use snps_names if it contains last gen info
binomial_reg_id = dict_snps_lastgen[dict_snps_lastgen['total_alleles05filter_lastgen'].notna()]['id'].reset_index(drop=True)

# Add SNP IDs to the binomial regression results
binomial_reg = pd.concat([binomial_reg, binomial_reg_id], axis=1)

# Map SNP IDs in binomial regression results to block IDs and calculate mean slope per block
binomial_reg['block'] = binomial_reg['id'].map(reverse_mapping)
block_binomial_slope = binomial_reg.groupby('block')['slope'].mean().reset_index()

# Read results from WZA analysis
# WZA (Weighted Z-score Analysis) is another method to assess association.
wza_binomial_regression_bio1 = pd.read_csv('../wza_last_gen/wza_binomial_regression_bio1_poly7.csv')

# Merge WZA results with the mean block slopes
# 'gene' in WZA results is assumed to correspond to 'block'.
wza_binomial_regression_bio1 = wza_binomial_regression_bio1.merge(block_binomial_slope, left_on='gene', right_on='block', how='left')

# Merge with the block origin data
wza_binomial_regression_bio1 = wza_binomial_regression_bio1.merge(block_origin_bio1, left_on='block', right_on='blocks', how='left')

# Read information about significant blocks
# This file identifies blocks that are statistically significant based on other analyses.
sign_blocks_info = pd.read_csv('../signficant_intersection/sign_blocks_union_first_last_gen_BH_final.csv')
sign_blocks_info = sign_blocks_info[sign_blocks_info['gen'].isin(['first_gen,last_gen', 'last_gen'])] # Filter for significance in last gen or both
sign_blocks_info = sign_blocks_info.drop_duplicates('block_id') # Remove duplicates

# Add a 'sign' column to indicate if a block is significant
wza_binomial_regression_bio1['sign'] = wza_binomial_regression_bio1['gene'].isin(sign_blocks_info['block_id'])

# Read gene information with significance (likely includes gene names)
# This file is used to get gene names for significant blocks to add as labels on the plot.
sign_genes_info = pd.read_csv('../signficant_intersection/genes_info_BH_tair10.csv')
sign_genes_info = sign_genes_info[sign_genes_info['gen']!='first_gen'] # Filter out entries from the first generation.
sign_genes_info = sign_genes_info.drop_duplicates('block_id') # Remove duplicate entries based on 'block_id'.

# Merge with the significant gene information to get details like gene name for significant blocks
wza_binomial_regression_bio1 = wza_binomial_regression_bio1.merge(sign_genes_info[['block_id', 'gene_name']], left_on='gene', right_on='block_id', how='left')

# --- Linear Regression Analysis ---
# Fit a linear regression model to assess the overall relationship
X = wza_binomial_regression_bio1['snp_origin_bio1'].dropna() # Independent variable: Block origin (temperature association)
y = wza_binomial_regression_bio1.loc[X.index, 'slope'].dropna() # Dependent variable: Block binomial slope (allele frequency change association)
# Ensure X and y have matching indices after dropping NaNs
common_index = X.index.intersection(y.index)
X = X.loc[common_index]
y = y.loc[common_index]


if not X.empty and not y.empty:
    X_with_intercept = sm.add_constant(X) # Add an intercept
    model = sm.OLS(y, X_with_intercept).fit() # Fit the OLS model

    # Extract and print regression parameters
    slope = model.params[1]
    intercept = model.params[0]
    p_value = model.pvalues[1]
    r_squared = model.rsquared
    print(f"Linear Regression Results:\n Slope: {slope:.4f}\n Intercept: {intercept:.4f}\n P-value: {p_value:.4f}\n R-squared: {r_squared:.4f}")
else:
    print("Not enough data to perform linear regression after dropping NaNs.")


# --- Plotting ---

# Create figure and axes for the hexbin plot
fig, ax = plt.subplots(figsize=(9, 6))

# Create a hexbin plot
# This visualizes the density of blocks based on their origin and allele frequency change association.
hb = ax.hexbin(
    wza_binomial_regression_bio1['snp_origin_bio1'], # X-axis: Block origin (temperature association)
    wza_binomial_regression_bio1['slope'], # Y-axis: Block binomial slope (allele frequency change association)
    gridsize=50,  # Size of the hexagons
    cmap='Greys_r',  # Color map (reversed grayscale)
    bins='log',  # Color based on log10 of count in each hexagon
    mincnt=1  # Minimum count to display a hexagon
)

# Add a color bar for the hexbin plot
cb = fig.colorbar(hb, ax=ax)
cb.set_label('log10(count)') # Label for the color bar

# Filter data for significant blocks
significant_blocks = wza_binomial_regression_bio1[wza_binomial_regression_bio1['sign'] == True].dropna(subset=['snp_origin_bio1', 'slope', 'gene_name'])

# Overlay scatter plot for significant blocks
# These points are highlighted in green.
ax.scatter(
    significant_blocks['snp_origin_bio1'],
    significant_blocks['slope'],
    color='Green',  # Color for significant points
    alpha=1,  # Opacity
    edgecolor='none',  # No edge color
    label='Significant' # Label for legend
)

# Add text labels for significant blocks
texts = []
for i, row in significant_blocks.iterrows():
    # Ensure gene_name is not NaN before adding text
    if pd.notna(row['gene_name']):
        texts.append(ax.text(
            row['snp_origin_bio1'],
            row['slope'],
            row['gene_name'],
            ha='center',
            va='bottom',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.8) # White background for text
        ))

# Adjust text labels to prevent overlap
if texts: # Only attempt to adjust if there are text objects
    adjust_text(texts,
                expand_points=(4, 4),
                only_move={'text':'y+'},
                arrowprops=dict(arrowstyle='->', color='gray', lw=1, alpha=0.8))

# Add grid lines and hide spines for cleaner look
plt.grid(True, color='lightgrey', alpha=0.7, zorder=0)
for spine in ['left', 'bottom', 'top', 'right']:
    ax.spines[spine].set_visible(False)

# Set plot labels
ax.set_xlabel('Block Origin (Temperature Association)')
ax.set_ylabel('Allele Frequency Change (Binomial Slope)')
ax.set_title('Block Origin vs. Allele Frequency Change') # Add a title for clarity

# Add a legend
ax.legend(title='Block Significance')

# Save the plot to a PDF file
plt.savefig('blocks_change_across_space_bio1.pdf')

# Display the plot
plt.show()

