import pandas as pd
import argparse
import statsmodels.api as sm
from statsmodels.formula.api import glm
import pickle
import numpy as np
# Set up argument parser
# Create the argument parser
parser = argparse.ArgumentParser(description='')

# Add 'partition' argument
parser.add_argument('partition', type=str, help='')

# Add the new 'split' argument
parser.add_argument('split', type=int, help='')

# Parse arguments
args = parser.parse_args()

# Access the arguments
partition = int(args.partition)
split = int(args.split)

file_path = '../jacknife_lfmm/splits_samples_last_gen_no33.pkl'

# Open and load the .pkl file
with open(file_path, 'rb') as file:
    samples = pickle.load(file)

train = samples[split][0]
test = samples[split][1]

##get the sites to train/test 
sites_u_train = pd.Series(train).str.split('_').str[0].unique()
sites_u_test = pd.Series(test).str.split('_').str[0].unique()
# conver tthem to int for easiest manipulations
sites_u_train_int =  [int(i) for i in sites_u_train ]
sites_u_test_int =  [int(i) for i in sites_u_test ]

## setting up the train env 
env_site_scaled = pd.read_csv('env_site_scaled.csv')
env_variable_train = env_site_scaled[env_site_scaled['site'].isin(sites_u_train_int)]['env_scaled'].reset_index(drop=True)

## setting up the test environemtn 
env_variable_test = env_site_scaled[env_site_scaled['site'].isin(sites_u_test_int)]
env_site_scaled_unique = env_variable_test.drop_duplicates().set_index('site')

partition_0 = pd.read_csv(f'../baypass_lastgen/individual_gfiles_last_gen/partition_{partition}.txt',header=None) 

pickle_file_path = f'../baypass_lastgen/individual_gfiles_last_gen/column_names_partition_{partition}'
with open(pickle_file_path, 'rb') as file:
    data0 = pickle.load(file)

partition_0.columns = data0

## eliminate site 33 

site_33_delete = [i for i in partition_0.columns if i.startswith('33_')]

partition_0 = partition_0.drop(site_33_delete,axis=1)

## eliminate site 33 

## filter only thre training environments 
partition_0 = partition_0[[i for i in partition_0.columns if any(i.startswith(site + '_') for site  in sites_u_train)]]
## filter only thre training environments 

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
    X = sm.add_constant(env_variable_train)  # Adding constant for intercept
    y = pd.concat([successes,failures],axis=1)
    # Fit the modelb
    model = sm.GLM(y, X, family=sm.families.Binomial())
    result = model.fit()

    # Extract slope (coefficient for environmental variable) and p-value
    slope = result.params[1]  # Coefficient for the environmental variable
    p_value = result.pvalues[1]  # P-value for the environmental variable
    intercept = result.params[0] 
    
    all_sites_pred_proba = []
    for site, row in env_site_scaled_unique.iterrows():
        # Predict the probability using the logistic model
        X_new = np.column_stack((np.ones(1), np.array([row['env_scaled']])))
        predicted_probability = result.predict(X_new)  # This gives the log-odds
        all_sites_pred_proba.append(predicted_probability[0])
    # Append the results to the DataFrame

    all_sites_pred_proba.append(slope)
    all_sites_pred_proba.append(p_value)
    all_sites_pred_proba.append(intercept)
    coefficients_df[i] = all_sites_pred_proba

coefficients_df = pd.DataFrame(coefficients_df).T
coefficients_df.columns = [sites_u_test[0], 'slope', 'pvalue','intercept']

coefficients_df.to_csv(f'results/split_{split}/partition{partition}.csv',index=None)