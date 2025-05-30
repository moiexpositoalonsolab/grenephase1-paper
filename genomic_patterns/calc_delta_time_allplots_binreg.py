#!/usr/bin/env python
# coding: utf-8

# This script calculates the change in allele frequency over generations for each SNP within a specified site
# using binomial regression and saves the results to a JSON Lines file.
# It processes data across all plots within the site.
# It takes the site identifier as a command-line argument.

# --- Required Input Files ---
# Ensure these files are present relative to the script's location:
# - ../key_files/merged_sample_table.csv
# - ../key_files/var_pos_grenenet.csv
# - ../key_files/merged_hapFIRE_allele_frequency.txt
# The script will write a log file and results file to the 'binom_reg/' directory.

# --- Import Libraries ---

import pandas as pd
import statsmodels.api as sm
import json
import argparse
import logging
import os # Import os for directory creation

# --- Argument Parsing ---

# Set up the argument parser to read the site identifier from the command line.
parser = argparse.ArgumentParser(description='Calculate binomial regression of allele frequency change over generations for a site (all plots).')

# Add argument for the site identifier.
parser.add_argument('site', type=str, help='The site identifier to process')

# Parse the command line arguments.
args = parser.parse_args()

# Extract the site identifier.
site = args.site
print(f"Processing site: {site}") # Print the site being processed

# --- Logging Setup ---

# Create the binomial regression results/log directory if it doesn't exist
output_dir = 'binom_reg'
os.makedirs(output_dir, exist_ok=True)

# Set up logging to a file specific to the site being processed within the binom_reg directory.
# Log messages at DEBUG level and higher will be recorded.
log_file_path = os.path.join(output_dir, f'log_site_{site}.log')
logging.basicConfig(filename=log_file_path, level=logging.DEBUG, format='%(asctime)s:%(levelname)s:%(message)s')

logging.info(f'Starting script for site: {site}')

# --- Data Loading and Preparation ---

# Read the merged sample table to get flower counts for each sample.
flowers = pd.read_csv('../key_files/merged_sample_table.csv')
# Create a dictionary mapping sample names to total flower counts.
num_flowers_map = flowers.set_index('sample_name')['total_flower_counts'].to_dict()
logging.debug('Sample table loaded and flower count map created.')

# Read SNP position information and get the mask for filtering based on 0.05 MAF.
var_pos = pd.read_csv('../key_files/var_pos_grenenet.csv')
mask = var_pos['maf05filter'].notna()
logging.debug('SNP position and MAF filter loaded.')

# Read only the header of the allele frequency file to identify columns (samples).
# This is done to find samples belonging to the current site and plots.
merged_hapFIRE_allele_frequency_header = pd.read_csv('../key_files/merged_hapFIRE_allele_frequency.txt', nrows=1, sep = '\t')
merged_hapFIRE_allele_frequency_header = merged_hapFIRE_allele_frequency_header.T.reset_index()

# Identify sample names that start with the current site identifier followed by an underscore.
site_plot_samples = [i for i in merged_hapFIRE_allele_frequency_header['index'] if i.startswith(str(site) + '_')]

logging.debug(f'Identified samples for site {site}: {site_plot_samples}')

# Read the allele frequency data, using only the columns (samples) identified for the current site.
merged_hapFIRE_allele_frequency = pd.read_csv('../key_files/merged_hapFIRE_allele_frequency.txt', sep = '\t', usecols = site_plot_samples)

# Apply the MAF filter mask to the allele frequency data.
merged_hapFIRE_allele_frequency = merged_hapFIRE_allele_frequency[mask].reset_index(drop=True)

logging.debug(f'Data loaded and filtered for site {site}, shape: {merged_hapFIRE_allele_frequency.shape}')

# --- Allele Count Calculation ---

# Calculate alternative and reference allele counts for each sample and SNP.
# Allele counts are derived from allele frequencies and total flower counts (assuming diploidy, hence * 2).
allele_counts_alt = {}
allele_counts_ref = {}
for col in merged_hapFIRE_allele_frequency.columns:
    sample_name = col
    num_flowers = num_flowers_map.get(sample_name, 0) # Use get with default 0 if sample not in map
    # Calculate alternative allele counts: frequency * total alleles (2 * flowers)
    allele_counts_alt[sample_name] = merged_hapFIRE_allele_frequency[col] * num_flowers * 2
    # Calculate reference allele counts: (1 - frequency) * total alleles (2 * flowers)
    maj_allele_freq = 1 - merged_hapFIRE_allele_frequency[col]
    allele_counts_ref[sample_name] = maj_allele_freq * num_flowers * 2

# Convert allele count dictionaries to DataFrames and transpose
allele_counts_alt_df = pd.concat(allele_counts_alt, axis=1).T.reset_index()
allele_counts_ref_df = pd.concat(allele_counts_ref, axis=1).T.reset_index()

# Reshape DataFrames and rename columns for merging
allele_counts_alt_melted = allele_counts_alt_df.melt(id_vars=['index'])
allele_counts_alt_melted.columns = ['sample', 'snp', 'count_alt']

allele_counts_ref_melted = allele_counts_ref_df.melt(id_vars=['index'])
allele_counts_ref_melted.columns = ['sample', 'snp', 'count_ref']

# Extract generation and plot information from sample names and add as columns
allele_counts_alt_melted['gen'] = allele_counts_alt_melted['sample'].str.split('_').str[1].astype(int)
# 'plot' is extracted but not used in the regression grouped by SNP across all plots
allele_counts_alt_melted['plot'] = allele_counts_alt_melted['sample'].str.split('_').str[2].astype(int)

allele_counts_ref_melted['gen'] = allele_counts_ref_melted['sample'].str.split('_').str[1].astype(int)
allele_counts_ref_melted['plot'] = allele_counts_ref_melted['sample'].str.split('_').str[2].astype(int)

# Merge the alternative and reference allele counts DataFrames on sample and snp
allele_counts = pd.merge(allele_counts_alt_melted, allele_counts_ref_melted[['sample', 'snp', 'count_ref']], on=['sample', 'snp'])

# Drop the original sample name column as it's no longer needed for analysis
allele_counts = allele_counts.drop('sample', axis=1)

logging.debug(f'Allele counts calculated and merged, shape: {allele_counts.shape}')

# --- Binomial Regression Analysis and Output ---

# Define the output file path for binomial regression results within the binom_reg directory.
output_jsonl_path = os.path.join(output_dir, f'results_site_{site}_br.jsonl')

# Open the output file in append mode to add to existing results without overwriting.
with open(output_jsonl_path, 'a') as file:
    buffer = [] # Buffer to hold results before writing in chunks
    buffer_chunk_size = 100 # Size of the buffer chunk

    # Group data by SNP to perform binomial regression for each SNP across all samples in the site.
    for snp_id, group in allele_counts.groupby('snp'):
        # Ensure there is enough data points (generations/plots) for regression (> 1 data point)
        # and that there is variation in generation.
        if len(group) > 1 and group['gen'].nunique() > 1:
            # Prepare data for binomial regression.
            # The independent variable X is the generation.
            X = sm.add_constant(group['gen'])  # Add a constant for the intercept term.
            # The dependent variable y is the allele counts as a two-column array [count_alt, count_ref].
            y = group.loc[:, ['count_alt', 'count_ref']].values # Select count columns and convert to numpy array

            # Fit the Generalized Linear Model (GLM) with a Binomial family.
            # This models the log-odds of the alternative allele count relative to the total count across generations.
            try:
                model = sm.GLM(y, X, family=sm.families.Binomial())
                result = model.fit()

                # Extract the slope (coefficient for generation) and its p-value.
                # The slope represents the estimated change in the log-odds of allele frequency per generation.
                # Handle cases where 'gen' is not in results (e.g., if only one generation is present after filtering)
                slope = result.params.get('gen', None) # Use .get to safely retrieve parameter
                p_value = result.pvalues.get('gen', None) # Use .get to safely retrieve p-value

                # Prepare the result as a dictionary.
                regression_result = {
                    "snp": str(snp_id),
                    "slope": round(slope, 5) if slope is not None else None, # Round if not None
                    "p_value": round(p_value, 5) if p_value is not None else None, # Round if not None
                }

                # Append result to buffer.
                buffer.append(regression_result)

                # Write buffered items to the file in chunks to improve efficiency.
                if len(buffer) >= buffer_chunk_size:
                    for item in buffer:
                        file.write(json.dumps(item) + '\n')
                    buffer.clear() # Clear the buffer after writing

            except Exception as e:
                logging.error(f"Error fitting binomial regression for SNP {snp_id}: {e}")
                # Continue to the next SNP if regression fails
                pass

    # Write any remaining buffered items to the file after the loop finishes.
    if buffer:
        for item in buffer:
            file.write(json.dumps(item) + '\n')
        buffer.clear()

logging.info(f'Script finished for site: {site}. Results saved to {output_jsonl_path}')
print(f'Script finished for site: {site}. Results saved to {output_jsonl_path}')




