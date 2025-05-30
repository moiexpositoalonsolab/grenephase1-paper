import pandas as pd
import os 
import random
import numpy as np
import scipy.stats
import statsmodels.api as sm
import matplotlib.pyplot as plt
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Process Kendall tau correlations for a specific run.')
parser.add_argument('run', type=str, help='The run identifier to process')

# Parse arguments
args = parser.parse_args()
run = int(args.run)

key_files = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/key_files/'
runs = pd.read_csv(key_files + 'runs_split_linages.csv')
clim_sites_during_exp = pd.read_csv('../key_files/bioclimvars_sites_era5_year_2018.csv')
env = runs.reset_index().merge(clim_sites_during_exp).drop(['index'],axis=1)
runs = env['run'].unique()

print(run)

samples_run = env[env['run']==run]['sample'].to_list()
env_run = env[env['run']==run]['bio1'].values

af = pd.read_csv('../key_files/delta_p_maf05_mincount05_firstgensamples.csv', usecols = samples_run)

af_ss = af[samples_run].copy()

kendall = {}
# Iterate through each row in af
for index, row in af_ss.iterrows():
    # Apply kendalltau to the current row and bio1
    geno_k_tau, geno_k_tau_p_value = scipy.stats.kendalltau(env_run, row)
    kendall[index] = [geno_k_tau,geno_k_tau_p_value]
    
kendall = pd.DataFrame(kendall).T
kendall.columns = ["K_tau","K_tau_p"]
kendall.to_csv(f'kendall_DELTAP_{run}.csv',index=False)