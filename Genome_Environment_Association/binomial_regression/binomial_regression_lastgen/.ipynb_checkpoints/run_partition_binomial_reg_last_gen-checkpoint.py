import pandas as pd
import argparse
import statsmodels.api as sm
from statsmodels.formula.api import glm
import pickle
import os

# Set up argument parser
parser = argparse.ArgumentParser(description='Process Kendall tau correlations for a specific run.')
parser.add_argument('partition', type=str, help='The partition identifier to process')
parser.add_argument('biovar', type=str, help='The environmental variable to use for the filename')

# Parse arguments
args = parser.parse_args()
biovar = args.biovar  # Store the biovar argument
partition = int(args.partition)

env_variable = pd.read_csv(f'env_{biovar}.csv')

partition_0 = pd.read_csv(f'../baypass_lastgen/individual_gfiles_last_gen/partition_{partition}.txt',header=None) 

pickle_file_path = f'../baypass_lastgen/individual_gfiles_last_gen/column_names_partition_{partition}'
with open(pickle_file_path, 'rb') as file:
    data0 = pickle.load(file)

partition_0.columns = data0

## eliminate site 33 

site_33_delete = [i for i in partition_0.columns if i.startswith('33_')]

partition_0 = partition_0.drop(site_33_delete,axis=1)

## eliminate site 33 

minor_columns = partition_0.filter(like='minor')
major_columns = partition_0.filter(like='major')

# Properly calculate total allele counts by summing corresponding minor and major allele counts
total_allele_counts = minor_columns.values + major_columns.values

# Create a new DataFrame for total counts with proper column names
total_allele_counts_df = pd.DataFrame(total_allele_counts, columns=minor_columns.columns.str.replace('_minor', '_total'))

minor_alleles = minor_columns.T.reset_index(drop=True)
major_columns = major_columns.T.reset_index(drop=True)

# Total allele counts as the number of trials (total attempts)
total_alleles = total_allele_counts_df.T.reset_index(drop=True)

successes = minor_alleles.iloc[:,1]
failures = major_columns.iloc[:,1]

# Running binomial regression for each SNP separately
coefficients_df = {}
print(env_variable.shape)
print(minor_alleles.shape)

for i in range(minor_alleles.shape[1]):
    # Prepare response (successes, failures) and predictors (environmental variable)
    successes = minor_alleles.iloc[:,i]
    failures = major_columns.iloc[:,i]
    # Set up the binomial regression model
    X = sm.add_constant(env_variable)  # Adding constant for intercept
    y = pd.concat([successes,failures],axis=1)
    # Fit the modelb
    model = sm.GLM(y, X, family=sm.families.Binomial())
    result = model.fit()

    # Extract slope (coefficient for environmental variable) and p-value
    slope = result.params[1]  # Coefficient for the environmental variable
    p_value = result.pvalues[1]  # P-value for the environmental variable
    
    # Append the results to the DataFrame
    coefficients_df[i] = [slope, p_value]

coefficients_df = pd.DataFrame(coefficients_df).T

coefficients_df.columns = ['slope', 'pvalue']


output_dir = f'results_{biovar}'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

coefficients_df.to_csv(f'{output_dir}/partition{partition}.csv', index=None)


