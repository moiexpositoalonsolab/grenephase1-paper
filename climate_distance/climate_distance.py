#!/usr/bin/env python
# coding: utf-8

"""
Analysis of Ecotype Frequency in Relation to Climatic Distance

This script analyzes the relationship between climatic distance and ecotype frequencies using two approaches:
1.  Leave-one-out cross-validation: Trains a linear regression model to predict ecotype frequency based on climatic distance
    and evaluates the predictions by leaving out one site at a time.
2.  Naive correlation: Directly correlates ecotype frequency with climatic distance without a predictive model.

Both approaches calculate various correlation metrics (Spearman, Pearson, Kendall's tau, R-squared) and assess the overlap
of top ecotypes between observed frequencies and model predictions/climatic distance.

Inputs:
- ../gwas/worldclim_ecotypesdata_sorted_20240517.csv: Worldclim data for ecotypes.
- ../key_files/bioclimvars_experimental_sites_era5.csv: ERA5 climate data for experimental sites.
- ../key_files/final_gen.csv: List of final generation samples.
- ../key_files/delta_ecotype_freq.txt: Delta ecotype frequencies.
- ../key_files/merged_ecotype_frequency.txt: Merged ecotype frequencies.
- ../jacknife_lfmm/splits_samples_last_gen_no33.pkl: Pickle file containing sample splits for cross-validation.

Outputs:
- climate_distance_model_results_mean_leave_1_out.csv: Cross-validation results for mean ecotype frequencies.
- climate_distance_model_results_leave_1_out_ecotype_freq_top.csv: Cross-validation results for individual samples.
- naive_climate_distance_top_ecotypes_mean.csv: Naive correlation results for mean ecotype frequencies.
- naive_climate_distance_top_ecotypes.csv: Naive correlation results for individual samples.
"""

# --- Import Libraries ---
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr, kendalltau
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import seaborn as sns
import allel
import pickle

# --- Load and Prepare Climate Data ---
# Load ecotype worldclim data and scale bio1
e_wc = pd.read_csv('../gwas/worldclim_ecotypesdata_sorted_20240517.csv', sep = '\t')[['ecotypeid','bio1']]
train_mean = np.mean(e_wc['bio1'])
train_std = np.std(e_wc['bio1'])
e_wc['bio1'] = (e_wc['bio1'] - train_mean) / train_std

# Load site climate data and scale bio1 using ecotype scaling parameters
sites_climate = pd.read_csv('../key_files/bioclimvars_experimental_sites_era5.csv')[['site', 'bio1']]
sites_climate['bio1'] = (sites_climate['bio1'] - train_mean) / train_std

# Calculate climatic distance matrix (squared difference in scaled bio1)
climate_ecotypes = e_wc['bio1'].to_numpy()
climate_sites = sites_climate['bio1'].to_numpy()
climatic_distance_matrix = (climate_ecotypes[:, np.newaxis] - climate_sites)**2

# Convert to DataFrame
distance_df = pd.DataFrame(climatic_distance_matrix, index=e_wc['ecotypeid'], columns=sites_climate['site'])

ecotypes_names = distance_df.index

# --- Load and Prepare Ecotype Frequency Data ---
# Load final generation samples and filter out '33_' samples
final_gen = pd.read_csv('../key_files/final_gen.csv')['sample_name']
final_gen = final_gen[~final_gen.str.startswith('33_')]
final_gen = final_gen.to_list()

unique_sites = pd.Series(final_gen).str.split('_').str[0].unique()

# Load delta and merged ecotype frequency data, keeping only final generation samples
delta_ef = pd.read_csv('../key_files/delta_ecotype_freq.txt', sep = '\t', usecols = final_gen)
ef = pd.read_csv('../key_files/merged_ecotype_frequency.txt', sep = '\t', usecols = final_gen)
delta_ef.index = ecotypes_names
ef.index = ecotypes_names

# --- Load Cross-validation Splits ---
# Path to pickle file containing sample splits
file_path = '../jacknife_lfmm/splits_samples_last_gen_no33.pkl'

# Open and load the .pkl file
with open(file_path, 'rb') as file:
    splits_samples_first_gen = pickle.load(file)

# --- Leave-one-out Cross-validation Analysis ---
print("Performing Leave-one-out Cross-validation...")
all_splits_predictions = {}
for split in range(len(splits_samples_first_gen)):
    # Get train and test samples/sites
    train_samples = splits_samples_first_gen[split][0]
    test_samples = splits_samples_first_gen[split][1]
    
    # Get unique sites for train and test splits
    sites_u_train = pd.Series(train_samples).str.split('_').str[0].unique()
    sites_u_test = pd.Series(test_samples).str.split('_').str[0].unique()
    
    # Convert sites to int for indexing
    sites_u_train_int =  [int(i) for i in sites_u_train ]
    sites_u_test_int =  [int(i) for i in sites_u_test ]
    
    # Get ecotype frequencies and climatic distances for train and test splits
    ef_train_samples = ef[[col for col in ef.columns if any(col.startswith(f'{site}_') for site in sites_u_train_int)]]
    cd_train_sites = distance_df[sites_u_train_int]
    cd_test_sites = distance_df[sites_u_test_int]
    
    # Melt dataframes for easy manipulation
    ef_train_melted = ef_train_samples.reset_index().melt(id_vars = ['ecotypeid'])
    cd_train_melted = cd_train_sites.reset_index().melt(id_vars = ['ecotypeid'])
    cd_test_melted = cd_test_sites.reset_index().melt(id_vars=['ecotypeid'])
    
    # Clean up melted dataframes
    ef_train_melted['site'] = ef_train_melted['variable'].str.split('_').str[0]
    ef_train_melted = ef_train_melted.drop('variable',axis=1)
    ef_train_melted.columns = ['ecotypeid', 'ef', 'site']
    
    cd_train_melted.columns = ['ecotypeid', 'site', 'cd']
    cd_train_melted['site'] = cd_train_melted['site'].astype(str)
    
    # Merge train data
    train_data = ef_train_melted.merge(cd_train_melted, on =['site', 'ecotypeid'])

    # Predict each ecotype at a time for the test site
    prediction_ecotypes = {}
    for ecotype_1 in train_data['ecotypeid'].unique():
        model_data = train_data[train_data['ecotypeid'] == ecotype_1]
        
        X = model_data['cd'].values.reshape(-1, 1)  # Independent variable (climatic distance)
        y = model_data['ef'].values # Dependent variable (ecotype frequency)
        
        # Fit the linear model
        model = LinearRegression().fit(X, y)
        
        # Get climatic distance for the current ecotype in the test site
        cd_ecotype = cd_test_melted[cd_test_melted['ecotypeid'] == ecotype_1]['value'].values[0]
        
        # Predict ecotype frequency for the test site
        y_pred = model.predict(cd_ecotype.reshape(1, -1) )
        prediction_ecotypes[ecotype_1] = y_pred
    
    # Store predictions for the current test site
    prediction_ecotypes_df = pd.DataFrame(prediction_ecotypes).T
    all_splits_predictions[sites_u_test[0]] = prediction_ecotypes_df

# Concatenate predictions from all splits
all_splits_predictions = pd.concat(all_splits_predictions,axis=1)
all_splits_predictions.columns = all_splits_predictions.columns.get_level_values(0)

# Calculate results for mean ecotype frequency
results_mean_cv = {}
for split in range(len(splits_samples_first_gen)):
    test_samples = splits_samples_first_gen[split][1]
    sites_u_test = pd.Series(test_samples).str.split('_').str[0].unique()
    sites_u_test_int =  [int(i) for i in sites_u_test ]
    
    # Get observed mean ecotype frequency for the test site
    ef_test_samples = ef[[col for col in ef.columns if any(col.startswith(f'{site}_') for site in sites_u_test_int)]]
    ef_mean_observed = ef_test_samples.mean(axis=1)
    ef_mean_observed.name = 'mean'
    
    # Select corresponding predictions
    prediction_to_test = all_splits_predictions[sites_u_test[0]]

    # Combine observed and predicted data
    conc = pd.concat([ef_mean_observed, prediction_to_test],axis=1)
    conc = conc.dropna() # Drop rows with NaN values
    
    # Calculate correlation metrics
    X_ranked = conc[sites_u_test[0]]
    y_ranked = conc['mean']
    
    sp_correlation, _ = spearmanr(X_ranked, y_ranked)
    pearsonr_value = pearsonr(X_ranked, y_ranked)[0]
    kendall_tau, _ = kendalltau(X_ranked, y_ranked)

    # Calculate R-squared
    X = conc[sites_u_test[0]].values.reshape(-1, 1)
    y = conc['mean'].values
    model = LinearRegression().fit(X, y)
    y_pred = model.predict(X)
    r_squared = r2_score(y, y_pred)
    
    # Calculate top ecotype overlap percentages
    res_list = [sp_correlation, pearsonr_value, r_squared, kendall_tau]
    for i in [10,20,30,40,50]:
        top_observed_indices = conc['mean'].sort_values().tail(i).index
        top_predicted_indices = conc[sites_u_test[0]].sort_values().tail(i).index
        common_top_indices = len(set(top_observed_indices).intersection(set(top_predicted_indices)))
        top_match_percentage = (common_top_indices / i) * 100
        res_list.append(top_match_percentage)
        
    results_mean_cv[sites_u_test[0]] = res_list

# Convert results to DataFrame and clean up columns
results_mean_cv = pd.DataFrame(results_mean_cv).T.reset_index()
results_mean_cv.columns = ['site', 'sp_correlation', 'pearsonr', 'r_squared','kendall_tau', 'top_10','top_20','top_30','top_40','top_50']

# Save results for mean ecotype frequency (cross-validation)
results_mean_cv.to_csv('climate_distance_model_results_mean_leave_1_out.csv',index=None)
print(f"Saved cross-validation results for mean ecotype frequency to climate_distance_model_results_mean_leave_1_out.csv")

# Calculate results for individual samples
results_cv = {}
for split in range(len(splits_samples_first_gen)):
    test_samples = splits_samples_first_gen[split][1]
    sites_u_test = pd.Series(test_samples).str.split('_').str[0].unique()
    sites_u_test_int =  [int(i) for i in sites_u_test ]
    
    # Get observed ecotype frequencies for the test site samples
    ef_test_samples = ef[[col for col in ef.columns if any(col.startswith(f'{site}_') for site in sites_u_test_int)]]

    # Select corresponding predictions
    prediction_to_test = all_splits_predictions[sites_u_test[0]]
    
    for sample in ef_test_samples.columns:
        ef_sample_observed = ef_test_samples[sample]
        
        # Combine observed and predicted data
        conc = pd.concat([ef_sample_observed, prediction_to_test],axis=1)
        conc = conc.dropna() # Drop rows with NaN values

        # Calculate correlation metrics
        X_ranked = conc[sites_u_test[0]]
        y_ranked = conc[sample]

        sp_correlation, _ = spearmanr(X_ranked, y_ranked)
        pearsonr_value = pearsonr(X_ranked, y_ranked)[0]
        kendall_tau, _ = kendalltau(X_ranked, y_ranked)

        # Calculate R-squared
        X = conc[sites_u_test[0]].values.reshape(-1, 1)
        y = conc[sample].values
        model = LinearRegression().fit(X, y)
        y_pred = model.predict(X)
        r_squared = r2_score(y, y_pred)
        
        # Calculate top ecotype overlap percentages
        res_list = [sp_correlation, pearsonr_value, r_squared, kendall_tau]
        for i in [10,20,30,40,50]:
            top_observed_indices = conc[sample].sort_values().tail(i).index
            top_predicted_indices = conc[sites_u_test[0]].sort_values().tail(i).index
            common_top_indices = len(set(top_observed_indices).intersection(set(top_predicted_indices)))
            top_match_percentage = (common_top_indices / i) * 100
            res_list.append(top_match_percentage)
            
        results_cv[sample] = res_list

# Convert results to DataFrame and clean up columns
results_cv = pd.DataFrame(results_cv).T.reset_index()
results_cv.columns = ['index', 'sp_correlation', 'pearsonr', 'r_squared','kendall_tau', 'top_10','top_20','top_30','top_40','top_50']
results_cv['site'] = results_cv['index'].str.split('_').str[0]
results_cv['plot'] = results_cv['index'].str.split('_').str[2]

# Save results for individual samples (cross-validation)
results_cv.to_csv('climate_distance_model_results_leave_1_out_ecotype_freq_top.csv',index=None)
print(f"Saved cross-validation results for individual samples to climate_distance_model_results_leave_1_out_ecotype_freq_top.csv")

# --- Naive Correlation Analysis ---
print("Performing Naive Correlation Analysis...")
# Calculate naive correlation results for mean ecotype frequency
results_mean_naive = {}
for site in unique_sites:
    # Get mean ecotype frequency and climatic distance for the site
    ef1 = ef[[col for col in ef.columns if col.startswith(f'{site}_')]]
    cd = distance_df[int(site)]
    ef_mean = ef1.mean(axis=1)
    ef_mean.name = 'mean'
    
    # Combine data and drop NaN
    conc = pd.concat([ef_mean, cd],axis=1)
    conc = conc.sort_values('mean')
    conc = conc.dropna()
    
    # Calculate correlation metrics
    X_ranked = conc['mean']
    y_ranked = conc[int(site)]

    sp_correlation, _ = spearmanr(X_ranked, y_ranked)
    pearsonr_value = pearsonr(X_ranked, y_ranked)[0]
    kendall_tau, _ = kendalltau(X_ranked, y_ranked)

    # Calculate R-squared
    X = conc[int(site)].values.reshape(-1, 1)
    y = conc['mean'].values
    model = LinearRegression().fit(X, y)
    y_pred = model.predict(X)
    r_squared = r2_score(y, y_pred)
    
    # Calculate top ecotype overlap percentages
    res_list = [sp_correlation, pearsonr_value, r_squared, kendall_tau]
    for i in [10,20,30,40,50]:
        # Note: for naive correlation, top ecotypes are expected to be associated with smaller climate distance
        top_observed_indices = conc['mean'].sort_values().tail(i).index
        top_predicted_indices = conc[int(site)].sort_values().head(i).index
        common_top_indices = len(set(top_observed_indices).intersection(set(top_predicted_indices)))
        top_match_percentage = (common_top_indices / i) * 100
        res_list.append(top_match_percentage)
        
    results_mean_naive[site] = res_list
    
# Convert results to DataFrame and clean up columns
results_mean_naive = pd.DataFrame(results_mean_naive).T.reset_index()
results_mean_naive.columns = ['site', 'sp_correlation', 'pearsonr', 'r_squared','kendall_tau', 'top_10','top_20','top_30','top_40','top_50']

# Invert correlation signs for interpretation (higher value = better match)
results_mean_naive['sp_correlation'] = results_mean_naive['sp_correlation'] *-1
results_mean_naive['pearsonr'] = results_mean_naive['pearsonr']  *-1
results_mean_naive['kendall_tau'] = results_mean_naive['kendall_tau']  *-1

# Save results for mean ecotype frequency (naive correlation)
results_mean_naive.to_csv('naive_climate_distance_top_ecotypes_mean.csv',index=None)
print(f"Saved naive correlation results for mean ecotype frequency to naive_climate_distance_top_ecotypes_mean.csv")

# Calculate naive correlation results for individual samples
results_naive = {}
for site in unique_sites:
    # Get ecotype frequencies and climatic distance for the site samples
    ef1 = ef[[col for col in ef.columns if col.startswith(f'{site}_')]]
    cd = distance_df[int(site)]
    
    for sample in ef1.columns:
        ef_sample = ef1[sample]
        
        # Combine data and drop NaN
        conc = pd.concat([ef_sample, cd],axis=1)
        conc = conc.sort_values(sample)
        conc = conc.dropna()

        # Calculate correlation metrics
        X_ranked = conc[sample]
        y_ranked = conc[int(site)]

        sp_correlation, _ = spearmanr(X_ranked, y_ranked)
        pearsonr_value = pearsonr(X_ranked, y_ranked)[0]
        kendall_tau, _ = kendalltau(X_ranked, y_ranked)

        # Calculate R-squared
        X = conc[int(site)].values.reshape(-1, 1)
        y = conc[sample].values
        model = LinearRegression().fit(X, y)
        y_pred = model.predict(X)
        r_squared = r2_score(y, y_pred)
        
        # Calculate top ecotype overlap percentages
        res_list = [sp_correlation, pearsonr_value, r_squared, kendall_tau]
        for i in [10,20,30,40,50]:
            # Note: for naive correlation, top ecotypes are expected to be associated with smaller climate distance
            top_observed_indices = conc[sample].sort_values().tail(i).index
            top_predicted_indices = conc[int(site)].sort_values().head(i).index
            common_top_indices = len(set(top_observed_indices).intersection(set(top_predicted_indices)))
            top_match_percentage = (common_top_indices / i) * 100
            res_list.append(top_match_percentage)
            
        results_naive[sample] = res_list
        
# Convert results to DataFrame and clean up columns
results_naive = pd.DataFrame(results_naive).T.reset_index()
results_naive.columns = ['index', 'sp_correlation', 'pearsonr', 'r_squared','kendall_tau', 'top_10','top_20','top_30','top_40','top_50']

# Invert correlation signs for interpretation (higher value = better match)
results_naive['sp_correlation'] = results_naive['sp_correlation'] *-1
results_naive['pearsonr'] = results_naive['pearsonr']  *-1
results_naive['kendall_tau'] = results_naive['kendall_tau']  *-1
results_naive['site'] = results_naive['index'].str.split('_').str[0]
results_naive['plot'] = results_naive['index'].str.split('_').str[2]

# Save results for individual samples (naive correlation)
results_naive.to_csv('naive_climate_distance_top_ecotypes.csv',index=None)
print(f"Saved naive correlation results for individual samples to naive_climate_distance_top_ecotypes.csv")
