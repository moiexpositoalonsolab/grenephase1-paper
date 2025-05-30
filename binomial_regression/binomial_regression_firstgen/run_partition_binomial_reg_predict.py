#!/usr/bin/env python
# coding: utf-8

"""
Required input files:
- 'env_{biovar}.csv': Standardized environmental variable data for training sites
- 'env_scaled_test.csv': Standardized environmental variable data for prediction sites
- '../baypass_first_gen/individual_gfiles/partition_{partition}.txt': Partition data
- '../baypass_first_gen/individual_gfiles/column_names_partition_{partition}': Column names for the partition

High-level outline:
This script performs binomial regression analysis and makes predictions for new sites.
It extends the basic binomial regression by:
1. Training the model on the original sites
2. Using the trained model to predict allele frequencies for new sites
3. The script can be run with any bioclimatic variable (bio1-bio19) by specifying the biovar parameter

Key differences from run_partition_binomial_reg.py:
- This script includes prediction functionality for new sites
- It uses additional environmental data from 'env_scaled_test.csv'
- Output includes predicted probabilities for each SNP at each test site
- Results are organized by SNP and site for easier analysis

Output:
- Creates a directory 'results_{biovar}/split_{partition}' if it doesn't exist
- Saves predictions as 'results_{biovar}/split_{partition}/prediction_sites.csv'
- The output includes:
  * SNP identifier
  * Slope and p-value from the regression
  * Intercept term
  * Predicted probabilities for each test site
"""

import pandas as pd
import argparse
import statsmodels.api as sm
from statsmodels.formula.api import glm
import pickle
import os
import numpy as np

# Set up argument parser
parser = argparse.ArgumentParser(description='Predict allele frequencies for new sites using binomial regression.')
parser.add_argument('partition', type=str, help='The partition identifier to process')
parser.add_argument('biovar', type=str, help='The environmental variable to use for the filename')

# Parse arguments
args = parser.parse_args()
biovar = args.biovar  # Store the biovar argument
partition = int(args.partition)

# Read the standardized environmental variable data for training
env_variable = pd.read_csv(f'env_{biovar}.csv')

# Read the standardized environmental variable data for prediction sites
env_site_test_scaled = pd.read_csv('env_scaled_test.csv')
test_sites = env_site_test_scaled['site'].unique()

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

# Running binomial regression and making predictions for each SNP
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

    # Extract model parameters
    slope = result.params.iloc[1]  # Coefficient for the environmental variable
    p_value = result.pvalues.iloc[1]  # P-value for the environmental variable
    intercept = result.params.iloc[0]  # Intercept term
    
    # Make predictions for each test site
    for index, row in env_site_test_scaled.iterrows():
        site = row['site']
        X_new = np.column_stack((np.ones(1), np.array([row['env_scaled']])))
        predicted_probability = result.predict(X_new)  # Get the predicted probability
        
        # Store predictions
        pred = [site, predicted_probability[0], slope, p_value, intercept]
        coefficients_df[str(i) + '_' + str(site)] = pred

# Process and organize results
coefficients_df = pd.DataFrame(coefficients_df).T
coefficients_df = coefficients_df.reset_index()
coefficients_df.columns = ['snp','site', 'pred_proba', 'slope', 'pvalue','intercept']
coefficients_df['snp'] = coefficients_df['snp'].str.split('_').str[0]

# Pivot the data to organize by SNP and site
coefficients_df = coefficients_df.pivot_table(
    columns=['site'], 
    values=['pred_proba'], 
    index=['snp','slope', 'pvalue', 'intercept']
).reset_index()

# Flatten column names
flat_columns = [
    f"{col[0]}_{int(col[1])}" if isinstance(col[1], float) else col[0]
    for col in coefficients_df.columns
]
coefficients_df.columns = flat_columns

# Sort results by SNP
coefficients_df['snp'] = coefficients_df['snp'].astype(int)
coefficients_df = coefficients_df.sort_values('snp')

# Create output directory and save results
output_dir = f'results_{biovar}/split_{partition}'
os.makedirs(output_dir, exist_ok=True)
coefficients_df.to_csv(f'results_{biovar}/split_{partition}/prediction_sites.csv', index=False)




