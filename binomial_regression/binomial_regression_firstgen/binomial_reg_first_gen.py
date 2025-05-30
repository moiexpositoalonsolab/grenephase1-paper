#!/usr/bin/env python
# coding: utf-8

"""
Required input files:
- '../key_files/generation_1_sample_names.txt': Contains sample names for the first generation
- '../key_files/bioclimvars_sites_era5_year_2018.csv': Contains bioclimatic variables for experimental sites
- '../baypass_first_gen/individual_gfiles/partition_0.txt': Contains partition data
- '../baypass_first_gen/individual_gfiles/column_names_partition_0': Contains column names for the partition

High-level outline:
This script performs binomial regression analysis for the first generation data.
It processes genetic data to analyze the relationship between allele frequencies and environmental variables.
The script can be run with any bioclimatic variable (bio1-bio19) by modifying the biovar variable.

The script:
1. Reads and processes sample names and environmental data
2. Standardizes the environmental variable
3. Processes genetic data from partition 0
4. Performs binomial regression for each SNP
5. Generates QQ plots for p-value analysis
6. Saves results and intermediate files

Output files:
- 'env_{biovar}.csv': Standardized environmental variable
- 'sample_allele_counts.csv': Sample allele counts
- QQ plot visualization
- Results DataFrame with slopes and p-values
"""

import pandas as pd
import pickle
import statsmodels.api as sm
from statsmodels.formula.api import glm
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import os

# Set the bioclimatic variable to analyze
biovar = 'bio1'  # Can be modified to any bioclimatic variable (bio1-bio19)

# Read and process sample names
first_gen_samples = pd.read_csv('../key_files/generation_1_sample_names.txt', header=None)[0]
samples = first_gen_samples.to_list()

# Read and process environmental data
clim_sites_during_exp = pd.read_csv('../key_files/bioclimvars_sites_era5_year_2018.csv')
sites_af = pd.Series(samples).str.split('_').str[0].astype(int)
sites_af.name = 'site'
env = sites_af.reset_index().merge(clim_sites_during_exp).drop(['index'], axis=1)

# Get and standardize the environmental variable
env_variable = env[biovar]
scaler = StandardScaler()
env_variable_scaled = scaler.fit_transform(env_variable.values.reshape(-1, 1))

# Save standardized environmental variable
pd.DataFrame(env_variable_scaled).to_csv(f'env_{biovar}.csv', index=None)

# Read the standardized environmental variable
env_variable = pd.read_csv(f'env_{biovar}.csv')

# Read partition data and column names
partition_0 = pd.read_csv('../baypass_first_gen/individual_gfiles/partition_0.txt', header=None)
pickle_file_path = '../baypass_first_gen/individual_gfiles/column_names_partition_0'
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

# Running binomial regression for each SNP
coefficients_df = {}

for i in range(minor_alleles.shape[1]):
    # Prepare response (successes, failures) and predictors (environmental variable)
    successes = minor_alleles.iloc[:,i]
    failures = major_columns.iloc[:,i]
    # Set up the binomial regression model
    X = sm.add_constant(env_variable)  # Adding constant for intercept
    y = pd.concat([successes,failures], axis=1)
    # Fit the model
    model = sm.GLM(y, X, family=sm.families.Binomial())
    result = model.fit()

    # Extract slope (coefficient for environmental variable) and p-value
    slope = result.params[1]  # Coefficient for the environmental variable
    p_value = result.pvalues[1]  # P-value for the environmental variable
    
    # Append the results to the DataFrame
    coefficients_df[i] = [slope, p_value]

# Convert results to DataFrame
coefficients_df = pd.DataFrame(coefficients_df).T
coefficients_df.columns = ['slope', 'pvalue']

# Generate QQ plot for p-values
observed_quantiles = -np.log10(np.sort(coefficients_df['pvalue'].values))
expected_quantiles = -np.log10(np.linspace(1 / len(coefficients_df['pvalue']), 1, len(coefficients_df['pvalue'])))

plt.figure(figsize=(10, 6))
sns.scatterplot(x=expected_quantiles, y=observed_quantiles, edgecolor='b', facecolor='none', alpha=0.5)
plt.plot([min(expected_quantiles), max(expected_quantiles)], [min(expected_quantiles), max(expected_quantiles)], 'r--')
plt.xlabel("Expected -log10(p-values)")
plt.ylabel("Observed -log10(p-values)")
plt.title(f'QQ Plot for {biovar}')
plt.savefig(f'qq_plot_{biovar}.png')
plt.close()

# Save sample allele counts
partition_0.to_csv('sample_allele_counts.csv', index=None)
