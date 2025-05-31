#!/usr/bin/env python
# coding: utf-8

"""
Required input files:
- '../key_files/final_gen.csv': Contains sample names for the last generation
- '/carnegie/nobackup/scratch/tbellagio/grene/data/bioclimvars_experimental_sites_era5.csv': Contains bioclimatic variables for experimental sites
- '../baypass_lastgen/individual_gfiles_last_gen/partition_*.txt': Contains partition data for each partition
- '../baypass_lastgen/individual_gfiles_last_gen/column_names_partition_*': Contains column names for each partition
- '../baypass_lastgen/individual_gfiles_last_gen/loci_partition_*': Contains SNP loci information for each partition

High-level outline:
This script performs binomial regression analysis for multiple bioclimatic variables (bio1-bio19).
It processes each bioclimatic variable separately, standardizes the environmental data,
and creates SLURM job scripts to run the analysis in parallel across different partitions.
The script handles all 19 bioclimatic variables (bio1 through bio19) and can be modified
to process any subset of these variables.

Bioclimatic Variables:
The script can process any of any of the 19 bioclimatic variables
"""

import pandas as pd
import pickle
import statsmodels.api as sm
from statsmodels.formula.api import glm
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import numpy as np
import matplotlib.pyplot as pl
import os
import random
import subprocess
import time

# List of all available bioclimatic variables
biovars = ['bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8',
           'bio9', 'bio10', 'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16',
           'bio17', 'bio18', 'bio19']

# Read sample names and remove site 33
last_gen = pd.read_csv('../key_files/final_gen.csv')['sample_name']
last_gen = last_gen[~last_gen.str.startswith('33_')]
samples = last_gen.to_list()

# Read and process environmental data
clim_sites_during_exp = pd.read_csv('/carnegie/nobackup/scratch/tbellagio/grene/data/bioclimvars_experimental_sites_era5.csv')
sites_af = pd.Series(samples).str.split('_').str[0].astype(int)
sites_af.name = 'site'
env = sites_af.reset_index().merge(clim_sites_during_exp).drop(['index'], axis=1)

# Process each bioclimatic variable
for biovar in biovars:
    # Get and standardize the environmental variable
    env_variable = env[biovar]
    scaler = StandardScaler()
    env_variable_scaled = scaler.fit_transform(env_variable.values.reshape(-1, 1))
    
    # Save standardized environmental variable
    pd.DataFrame(env_variable_scaled).to_csv(f'env_{biovar}.csv', index=None)

# Get list of partitions
files = os.listdir('../baypass_lastgen/individual_gfiles_last_gen/')
partitions = [int(file.split('_')[1].replace('.txt', '')) for file in files if '.txt' in file]
partitions.sort()

# Create and submit SLURM jobs for each bioclimatic variable and partition
shfiles = []
for biovar in biovars:
    for partition in partitions:
        seed = random.randint(1, 100000000)
        file = f'shfiles/partition_{partition}_{biovar}.sh'
        cmd = f'python run_partition_binomial_reg_last_gen.py {partition} {biovar}'
        
        # Create SLURM script
        text = f'''#!/bin/bash
#SBATCH --job-name=run_partition_binomial_reg{partition}
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --output=run_partition_binomial_reg{partition}_{biovar}_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

module load python/3.11_conda
conda activate /home/tbellagio/miniforge3/envs/pipeline_snakemake
export LD_LIBRARY_PATH="/home/tbellagio/miniforge3/envs/run_baypass/lib:$LD_LIBRARY_PATH"
cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/binomial_regression_lastgen
{cmd}
'''
        with open(file, 'w') as o:
            o.write("%s" % text)
        shfiles.append(file)

# Submit jobs in batches of 100
for i in range(0, len(shfiles), 100):
    batch = shfiles[i:i+100]
    for shfile in batch:
        subprocess.run(["sbatch", shfile], check=True)
    if i + 100 < len(shfiles):
        time.sleep(120)

# Process results for each bioclimatic variable
for biovar in biovars:
    print(f"Processing results for {biovar}")
    partitions_r = {}
    for i in range(len(partitions)):
        # Load SNP loci information
        pickle_file_path = f'../baypass_lastgen/individual_gfiles_last_gen/loci_partition_{i}'
        with open(pickle_file_path, 'rb') as file:
            loci_f = pickle.load(file)
        
        # Load and process results
        results = pd.read_csv(f'results_{biovar}/partition{i}.csv')
        results['snp_id'] = loci_f
        partitions_r[i] = results

    # Combine results from all partitions
    results = pd.concat(partitions_r).reset_index(drop=True)
    results.to_csv(f'{biovar}_binomial_reg_results_last_gen.csv', index=None)

    # Process top hits
    binomf = results.sort_values('pvalue').head(100)
    
    # Load block information
    with open('../key_files/blocks_snpsid_dict.pkl', 'rb') as file:
        dict_blocks = pickle.load(file)
    reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}
    
    # Add block information and remove duplicates
    binomf['block'] = binomf['snp_id'].map(reverse_mapping)
    binomf = binomf.drop_duplicates('block')
    
    # Save top hits
    binomf.to_csv(f'top_hits_binom_last_gen_{biovar}.csv')
