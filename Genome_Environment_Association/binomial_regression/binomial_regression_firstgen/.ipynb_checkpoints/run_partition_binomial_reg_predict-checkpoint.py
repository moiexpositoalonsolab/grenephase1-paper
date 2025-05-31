import pandas as pd
import argparse
import statsmodels.api as sm
from statsmodels.formula.api import glm
import pickle
import os
import numpy as np
# Set up argument parser
# Set up argument parser
parser = argparse.ArgumentParser(description='predcit binom')
parser.add_argument('partition', type=str, help='The partition identifier to process')
parser.add_argument('biovar', type=str, help='The environmental variable to use for the filename')

# Parse arguments
args = parser.parse_args()
biovar = args.biovar  # Store the biovar argument
partition = int(args.partition)

env_variable = pd.read_csv(f'env_{biovar}.csv')


## get the predict sites 
env_site_test_scaled = pd.read_csv('env_scaled_test.csv')
test_sites = env_site_test_scaled['site'].unique()

partition_0 = pd.read_csv(f'../baypass_first_gen/individual_gfiles/partition_{partition}.txt',header=None) 

pickle_file_path = f'../baypass_first_gen/individual_gfiles/column_names_partition_{partition}'
with open(pickle_file_path, 'rb') as file:
    data0 = pickle.load(file)

partition_0.columns = data0

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
    slope = result.params.iloc[1]  # Coefficient for the environmental variable
    p_value = result.pvalues.iloc[1]  # P-value for the environmental variable
    intercept = result.params.iloc[0]  # Intercept term
    ### predict
    
    for index, row in env_site_test_scaled.iterrows():
        # Predict the probability using the logistic model
        site = row['site']
        X_new = np.column_stack((np.ones(1), np.array([row['env_scaled']])))
        predicted_probability = result.predict(X_new)  # Get the predicted probability
        # Store predictions as a dictionary
        #sites_list.append(site)
        #all_sites_pred_proba.append(predicted_probability[0])
        
        pred = [site,predicted_probability[0],slope,p_value,intercept]
        coefficients_df[str(i) + '_' + str(site)] = pred

coefficients_df = pd.DataFrame(coefficients_df).T

coefficients_df = coefficients_df.reset_index()

coefficients_df.columns = ['snp','site', 'pred_proba', 'slope', 'pvalue','intercept']

coefficients_df['snp'] = coefficients_df['snp'].str.split('_').str[0]

coefficients_df = coefficients_df.pivot_table(columns = ['site'], values = ['pred_proba'], index= ['snp','slope', 'pvalue', 'intercept']).reset_index()

flat_columns = [
    f"{col[0]}_{int(col[1])}" if isinstance(col[1], float) else col[0]
    for col in coefficients_df.columns
]

coefficients_df.columns = flat_columns

coefficients_df['snp'] = coefficients_df['snp'].astype(int)

coefficients_df = coefficients_df.sort_values('snp')

# Ensure the output directory exists
output_dir = f'results_{biovar}/split_{partition}'
os.makedirs(output_dir, exist_ok=True)


# Save the DataFrame to a CSV file
coefficients_df.to_csv(f'results_{biovar}/split_{partition}/prediction_sites.csv', index=False)




