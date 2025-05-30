#!/usr/bin/env python
# coding: utf-8

# This script calculates the change in allele counts over generations for individual SNPs
# across all plots within a specific site using Ordinary Least Squares (OLS) regression.
# It takes a site identifier as a command-line argument and outputs regression results
# (slope and p-value) for each SNP to a JSON Lines file.

# --- Required Input Files ---
# Ensure these files are present relative to the script's location:
# - ../key_files/merged_sample_table.csv
# - ../key_files/var_pos_grenenet.csv
# - ../key_files/merged_hapFIRE_allele_frequency.txt
# The script will write output to:
# - results_site_<site>.jsonl
# It will also create a log file:
# - log_site_<site>.log

# --- Import Libraries ---

import pandas as pd
import statsmodels.api as sm
import json
from statsmodels.formula.api import ols
import argparse
import logging
import os # Import os for directory creation

# --- Command Line Argument Parsing ---

# Set up the argument parser
parser = argparse.ArgumentParser(description='Calculate allele count change over time for all plots in a site.')

# Add argument for the site identifier
parser.add_argument('site', type=str, help='The site identifier to process')

# Parse the command line arguments
args = parser.parse_args()

# Extract the site value
site = args.site
print(f"Processing site: {site}")

# --- Logging Setup ---

# Create the log directory if it doesn't exist
log_dir = 'logs_delta_time'
os.makedirs(log_dir, exist_ok=True)

# Configure logging to a file specific to the site
log_file_path = os.path.join(log_dir, f'log_site_{site}.log')
logging.basicConfig(filename=log_file_path, level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

logging.info(f'Starting script for site: {site}')

# --- Data Loading and Preparation ---

# Read sample table information (for flower counts)
flowers = pd.read_csv('../key_files/merged_sample_table.csv')
# Create a dictionary mapping sample names to total flower counts
num_flowers_map = flowers.set_index('sample_name')['total_flower_counts'].to_dict()
logging.debug('Sample table loaded and flower count map created.')

# Read SNP position and filtering information
var_pos = pd.read_csv('../key_files/var_pos_grenenet.csv')
# Get a mask for filtering based on 0.05 MAF (Minor Allele Frequency)
mask = var_pos['maf05filter'].notna()
logging.debug('SNP position and MAF filter loaded.')

# Read allele frequency data and filter for relevant samples
# Read only the header initially to find columns for the specific site (all plots)
merged_hapFIRE_allele_frequency_header = pd.read_csv('../key_files/merged_hapFIRE_allele_frequency.txt', nrows=1, sep = '\t')
merged_hapFIRE_allele_frequency_header = merged_hapFIRE_allele_frequency_header.T.reset_index()
# Identify columns (samples) belonging to the specified site
site_samples = [i for i in merged_hapFIRE_allele_frequency_header['index'] if i.startswith(str(site) + '_')]
logging.debug(f'Identified {len(site_samples)} samples for site {site}.')

# Read the full allele frequency data, using only the selected columns
merged_hapFIRE_allele_frequency = pd.read_csv('../key_files/merged_hapFIRE_allele_frequency.txt', sep = '\t', usecols = site_samples)
logging.debug(f'Allele frequency data loaded for site {site}, shape: {merged_hapFIRE_allele_frequency.shape}')

# Apply MAF filter to the allele frequency data
merged_hapFIRE_allele_frequency = merged_hapFIRE_allele_frequency[mask].reset_index(drop=True)
logging.debug(f'MAF filter applied, shape: {merged_hapFIRE_allele_frequency.shape}')

# --- Allele Count Calculation ---

# Calculate allele counts for each sample
allele_counts_list = []
for col in merged_hapFIRE_allele_frequency.columns:
    # Calculate allele count: frequency * number of diploid individuals * 2
    # Assuming number of individuals is represented by flower counts / 2 (adjust if necessary)
    num_individuals = num_flowers_map.get(col, 0) # Use get with default 0 if sample not in map
    allele_count = merged_hapFIRE_allele_frequency[col] * num_individuals * 2
    allele_counts_list.append(allele_count.rename(col))

allele_counts_df = pd.concat(allele_counts_list, axis=1)

# Reshape DataFrame for regression analysis
allele_counts_df = allele_counts_df.T.reset_index()
allele_counts_df.columns = ['sample_name'] + allele_counts_df.columns[1:].tolist() # Rename the first column
allele_counts_melted = allele_counts_df.melt(id_vars=['sample_name'], var_name='snp', value_name='count')

# Extract generation and plot from sample names
allele_counts_melted['gen'] = allele_counts_melted['sample_name'].str.split('_').str[1].astype(int)
allele_counts_melted['plot'] = allele_counts_melted['sample_name'].str.split('_').str[2].astype(int) # Keep plot for grouping if needed later, though regression is per SNP across all plots
allele_counts_melted = allele_counts_melted.drop('sample_name', axis=1) # Drop sample name column

logging.debug(f"Prepared {len(allele_counts_melted)} allele counts for regression.")

# --- Perform OLS Regression and Save Results ---

# Create the results directory if it doesn't exist
results_dir = 'results_delta_time'
os.makedirs(results_dir, exist_ok=True)

# Define the output file path
output_file_path = os.path.join(results_dir, f'results_site_{site}.jsonl')

# Open the output file in append mode to add to existing results without overwriting
with open(output_file_path, 'a') as file:
    buffer = []
    # Group data by SNP to perform regression for each SNP across all samples in the site
    for snp, group in allele_counts_melted.groupby('snp'):
        # Ensure there's enough data (at least 2 data points with different generations) for regression
        if len(group) > 1 and group['gen'].nunique() > 1:
            # Perform OLS regression: count ~ gen
            try:
                model = ols('count ~ gen', data=group).fit()
                # Extract slope and p-value for the 'gen' (generation) coefficient
                slope = model.params['gen'].item()  # Convert numpy float64 to Python float
                p_value = model.pvalues['gen'].item()  # Convert numpy float64 to Python float
                
                # Prepare the result as a JSON object
                result = {
                    "snp": str(snp),
                    "slope": round(slope, 5), # Round slope for smaller file size
                    "p_value": round(p_value, 5) # Round p-value
                }
                
                # Append result to buffer
                buffer.append(result)
                
                # Write buffered items to the file in chunks
                if len(buffer) == 100:
                    for item in buffer:
                        file.write(json.dumps(item) + '\n')
                    buffer.clear()
            except Exception as e:
                logging.error(f"Error performing OLS for snp {snp}: {e}")
                # Optionally, log the error or handle it differently

    # Write any remaining items in the buffer
    if buffer:
        for item in buffer:
            file.write(json.dumps(item) + '\n')

logging.info(f"OLS regression results saved to {output_file_path}")
print(f"OLS regression results saved to {output_file_path}")




