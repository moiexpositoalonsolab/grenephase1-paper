#!/usr/bin/env python
# coding: utf-8

"""
Required input files:
- 'env_{biovar}.csv': Standardized environmental variable data
- '../baypass_first_gen/individual_gfiles/partition_{partition}.txt': Partition data
- '../baypass_first_gen/individual_gfiles/column_names_partition_{partition}': Column names for the partition

High-level outline:
This script performs binomial regression analysis for a specific partition and bioclimatic variable.
It processes genetic data to analyze the relationship between allele frequencies and environmental variables.
The script can be run with any bioclimatic variable (bio1-bio19) by specifying the biovar parameter.

The script:
1. Reads the standardized environmental variable data
2. Processes genetic data from a specific partition
3. Performs binomial regression for each SNP
4. Saves the results (slopes and p-values) to a CSV file

Output:
- Creates a directory 'results_{biovar}' if it doesn't exist
- Saves results as 'results_{biovar}/partition{partition}.csv'
"""

import pandas as pd
import argparse
import statsmodels.api as sm
from statsmodels.formula.api import glm
import pickle
import os

# Set up argument parser
parser = argparse.ArgumentParser(description='Perform binomial regression analysis for a specific partition and bioclimatic variable.')
parser.add_argument('partition', type=str, help='The partition identifier to process')
parser.add_argument('biovar', type=str, help='The environmental variable to use for the filename')

# Parse arguments
args = parser.parse_args()
biovar = args.biovar  # Store the biovar argument
partition = int(args.partition)

# Read the standardized environmental variable data
env_variable = pd.read_csv(f'env_{biovar}.csv')

# Read partition data and column names
partition_0 = pd.read_csv(f'../baypass_first_gen/individual_gfiles/partition_{partition}.txt',header=None) 

pickle_file_path = f'../baypass_first_gen/individual_gfiles/column_names_partition_{partition}'
with open(pickle_file_path, 'rb') as file:
    data0 = pickle.load(file)

partition_0.columns = data0

# Extract minor and major allele columns
minor_columns = partition_0.filter(like='minor')
major_columns = partition_0.filter(like='major')

# Calculate total allele counts
total_allele_counts = minor_columns.values + major_columns.values
total_allele_counts_df = pd.DataFrame(total_allele_counts, columns=minor_columns.columns.str.replace('_minor', '_total'))

# Prepare data for regression
minor_alleles = minor_columns.T.reset_index(drop=True)
major_columns = major_columns.T.reset_index(drop=True)
total_alleles = total_allele_counts_df.T.reset_index(drop=True)

# Running binomial regression for each SNP separately
coefficients_df = {}

for i in range(minor_alleles.shape[1]):
    # Prepare response (successes, failures) and predictors (environmental variable)
    successes = minor_alleles.iloc[:,i]
    failures = major_columns.iloc[:,i]
    # Set up the binomial regression model
    X = sm.add_constant(env_variable)  # Adding constant for intercept
    y = pd.concat([successes,failures],axis=1)
    # Fit the model
    model = sm.GLM(y, X, family=sm.families.Binomial())
    result = model.fit()

    # Extract slope (coefficient for environmental variable) and p-value
    slope = result.params[1]  # Coefficient for the environmental variable
    p_value = result.pvalues[1]  # P-value for the environmental variable
    
    # Append the results to the DataFrame
    coefficients_df[i] = [slope, p_value]

# Convert results to DataFrame and save
coefficients_df = pd.DataFrame(coefficients_df).T
coefficients_df.columns = ['slope', 'pvalue']

# Create output directory if it doesn't exist
output_dir = f'results_{biovar}'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Save results
coefficients_df.to_csv(f'{output_dir}/partition{partition}.csv', index=None)

