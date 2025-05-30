#!/usr/bin/env python
# coding: utf-8

"""
Required input files:
- 'env_{biovar}.csv': Standardized environmental variable data (created by binomial_reg_last_gen.py)
- '../baypass_lastgen/individual_gfiles_last_gen/partition_{partition}.txt': Partition data
- '../baypass_lastgen/individual_gfiles_last_gen/column_names_partition_{partition}': Column names for the partition

High-level outline:
This script performs binomial regression analysis for a specific partition and bioclimatic variable.
It is called by the main script (binomial_reg_last_gen.py) for each combination of partition and bioclimatic variable.
The script processes allele frequency data and performs binomial regression to test for associations
between allele frequencies and environmental variables.

Bioclimatic Variables:
This script can be run with any of the 19 bioclimatic variables (bio1-bio19):
"""

import pandas as pd
import argparse
import statsmodels.api as sm
from statsmodels.formula.api import glm
import pickle
import os

# Set up argument parser
parser = argparse.ArgumentParser(description='Process binomial regression for a specific partition and bioclimatic variable.')
parser.add_argument('partition', type=str, help='The partition identifier to process')
parser.add_argument('biovar', type=str, help='The bioclimatic variable to use (e.g., bio1, bio12, bio18)')

# Parse arguments
args = parser.parse_args()
biovar = args.biovar  # Store the biovar argument
partition = int(args.partition)

# Read standardized environmental variable
env_variable = pd.read_csv(f'env_{biovar}.csv')

# Read partition data
partition_0 = pd.read_csv(f'../baypass_lastgen/individual_gfiles_last_gen/partition_{partition}.txt', header=None)

# Load column names for the partition
pickle_file_path = f'../baypass_lastgen/individual_gfiles_last_gen/column_names_partition_{partition}'
with open(pickle_file_path, 'rb') as file:
    data0 = pickle.load(file)

partition_0.columns = data0

# Remove site 33 data
site_33_delete = [i for i in partition_0.columns if i.startswith('33_')]
partition_0 = partition_0.drop(site_33_delete, axis=1)

# Process allele counts
minor_columns = partition_0.filter(like='minor')
major_columns = partition_0.filter(like='major')

# Calculate total allele counts
total_allele_counts = minor_columns.values + major_columns.values
total_allele_counts_df = pd.DataFrame(total_allele_counts, 
                                    columns=minor_columns.columns.str.replace('_minor', '_total'))

# Prepare data for regression
minor_alleles = minor_columns.T.reset_index(drop=True)
major_columns = major_columns.T.reset_index(drop=True)
total_alleles = total_allele_counts_df.T.reset_index(drop=True)

# Print shapes for debugging
print(env_variable.shape)
print(minor_alleles.shape)

# Run binomial regression for each SNP
coefficients_df = {}

for i in range(minor_alleles.shape[1]):
    # Prepare response and predictors
    successes = minor_alleles.iloc[:,i]
    failures = major_columns.iloc[:,i]
    
    # Set up the model
    X = sm.add_constant(env_variable)  # Adding constant for intercept
    y = pd.concat([successes, failures], axis=1)
    
    # Fit the model
    model = sm.GLM(y, X, family=sm.families.Binomial())
    result = model.fit()

    # Extract coefficients
    slope = result.params[1]  # Coefficient for the environmental variable
    p_value = result.pvalues[1]  # P-value for the environmental variable
    
    # Store results
    coefficients_df[i] = [slope, p_value]

# Convert results to DataFrame
coefficients_df = pd.DataFrame(coefficients_df).T
coefficients_df.columns = ['slope', 'pvalue']

# Create output directory if it doesn't exist
output_dir = f'results_{biovar}'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Save results
coefficients_df.to_csv(f'{output_dir}/partition{partition}.csv', index=None)


