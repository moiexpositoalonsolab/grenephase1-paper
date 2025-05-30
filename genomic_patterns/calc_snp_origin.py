#!/usr/bin/env python
# coding: utf-8

# This script calculates a measure of "SNP origin" based on the association between SNP genotypes
# in the 1001 Genomes reference population and environmental data (Mean Annual Temperature - bio1).
# It then approximates the "origin" for predefined blocks of SNPs.

# --- Required Input Files ---
# Ensure these files are present relative to the script's location:
# - ../key_files/var_pos_grenenet.csv
# - ../key_files/blocks_snpsid_dict.pkl
# - ../key_files/1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv
# - ../key_files/founder_ecotype_names.csv
# - ../key_files/greneNet_final_v1.1_LDpruned.recode.vcf
# The script will write output to:
# - ../key_files/snp_origin_bio1_1001gvcf.csv (SNP-level origin)
# - ../key_files/block_origin_bio1_1001g.csv (Block-level origin)

# --- Import Libraries ---

import pandas as pd
import allel
import numpy as np
import pickle
import os # Import os for directory creation if needed

# --- Data Loading and Preparation ---

# Read SNP position information
snps_names = pd.read_csv('../key_files/var_pos_grenenet.csv')
# Filter SNPs based on total alleles in the first generation (assuming this indicates usable data)
snps_names = snps_names[snps_names['total_alleles05filter_firstgen'].notna()].reset_index(drop=True)

# Load the dictionary mapping SNPs to blocks
with open('../key_files/blocks_snpsid_dict.pkl', 'rb') as f:
    dict_blocks = pickle.load(f)
# Create a reverse mapping from SNP ID to block ID
reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}

# Read climate data for the 1001 Genomes accessions
clim1001 = pd.read_csv('../key_files/1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv')
clim1001['ecotypeid'] = clim1001['ecotypeid'].astype(str) # Ensure ecotypeid is string for merging

# Read the list of founder ecotypes (although not directly used in the final output files, it was in the original notebook)
founder_ecotypes = pd.read_csv('../key_files/founder_ecotype_names.csv')['0'].astype(str)

# Read genotype data from the VCF file
vcfgrene = allel.read_vcf('../key_files/greneNet_final_v1.1_LDpruned.recode.vcf')

# Extract genotype data (sum of alleles per individual per SNP)
geno = vcfgrene['calldata/GT'].sum(axis=2)
geno = pd.DataFrame(geno)
geno.columns = vcfgrene['samples'] # Assign sample names as column headers

# Filter climate data to include only ecotypes present in the genotype data
clim1001 = clim1001[clim1001['ecotypeid'].isin(geno.columns)].reset_index(drop=True)
clim1001 = clim1001[['ecotypeid', 'bio1']].set_index('ecotypeid')

# Transpose genotype data and align with climate data by sorting index
geno = geno.T.sort_index().reset_index(drop=True)
clim1001 = clim1001.sort_index().reset_index(drop=True)

# --- Calculate SNP Origin (Association with Temperature) ---

# Extract bio1 temperature data and mean center it
temperature_centered = clim1001['bio1'] - np.mean(clim1001['bio1'])

# Calculate the association between SNP genotypes and centered temperature
# This is done using a dot product of the transposed genotype matrix and the centered temperature vector.
snp_temperature_association = np.array(geno).T.dot(np.array(temperature_centered))

# Normalize the association by the number of ecotypes
snp_temperature_association /= len(clim1001)

# Convert the result to a pandas Series
snp_temperature_association = pd.Series(snp_temperature_association.flatten())

# --- Save SNP-level Origin Data ---

# Create a DataFrame with SNP information and the calculated origin
snp_origin_bio1_df = pd.DataFrame({'chrom': vcfgrene['variants/CHROM'], 'pos': vcfgrene['variants/POS'], 'snp_origin_bio1': snp_temperature_association})

# Define the output path for SNP-level origin data
snp_origin_output_path = '../key_files/snp_origin_bio1_1001gvcf.csv'

# Ensure the output directory exists (assuming ../key_files already exists based on inputs)
# os.makedirs(os.path.dirname(snp_origin_output_path), exist_ok=True)

# Save the SNP-level origin data to a CSV file
snp_origin_bio1_df.to_csv(snp_origin_output_path, index=None)

print(f"SNP-level origin data saved to {snp_origin_output_path}")

# --- Approximate Block Origin ---

# Add SNP IDs to the SNP-level origin DataFrame
snp_origin_bio1_df['id'] = snp_origin_bio1_df['chrom'].astype(str) + '_' + snp_origin_bio1_df['pos'].astype(str)

# Map SNPs to blocks using the reverse mapping dictionary
snp_origin_bio1_df['blocks'] = snp_origin_bio1_df['id'].map(reverse_mapping)

# Calculate the mean SNP origin for each block
# Drop blocks with no mapped SNPs (NaN) before grouping if necessary, although groupby handles NaNs by default.
block_origin_bio1_df = snp_origin_bio1_df.groupby('blocks')['snp_origin_bio1'].mean().reset_index()

# Define the output path for block-level origin data
block_origin_output_path = '../key_files/block_origin_bio1_1001g.csv'

# Ensure the output directory exists (assuming ../key_files already exists)
# os.makedirs(os.path.dirname(block_origin_output_path), exist_ok=True)

# Save the block-level origin data to a CSV file
block_origin_bio1_df.to_csv(block_origin_output_path, index=None)

print(f"Block-level origin data saved to {block_origin_output_path}")
