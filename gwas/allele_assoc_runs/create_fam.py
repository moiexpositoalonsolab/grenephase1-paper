#!/usr/bin/env python
# coding: utf-8

"""
Prepare Input Files and Submit Jobs for Climate GWAS

This script prepares input files for GEMMA and BSLMM analyses using bioclimatic variables
as phenotypes. It creates necessary directories, prepares .fam files, and submits jobs
to the cluster.

Required Files:
- bioclimvars_ecotypes_era5.csv: Climate data for ecotypes
- 1001g_grenet_climate.fam: Original .fam file with sample information
- 1001g_grenet_climate.bed: PLINK .bed file
- 1001g_grenet_climate.bim: PLINK .bim file
- 1001g_grenet_climateLDpruned_05maf.cXX.txt: Kinship matrix for GEMMA
- pc_1000.txt: Principal components for covariate adjustment

Outputs:
- For each bioclimatic variable (bio1-bio19):
  - A directory in lmm_gemma/ containing:
    - Modified .fam file with bioclimatic variable as phenotype
    - Symbolic links to .bed and .bim files
    - GEMMA job script
  - A directory in bslmm/ containing:
    - Modified .fam file with bioclimatic variable as phenotype
    - Symbolic links to .bed and .bim files
    - BSLMM job script
"""

# --- Import Libraries ---
import pandas as pd
import os
import subprocess 
import random

# --- Configuration ---
# Path to climate data file
climate_file = '../../data/bioclimvars_ecotypes_era5.csv'

# --- Load Climate Data ---
print("Loading climate data...")
climate = pd.read_csv(climate_file)

# Rename columns to include all bioclimatic variables
climate.columns = ['ecotype', 'bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 
                  'bio7', 'bio8', 'bio9', 'bio10', 'bio11', 'bio12', 'bio13', 
                  'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19']

# --- Load Original .fam File ---
print("Loading original .fam file...")
og_fam = pd.read_csv('1001g_grenet_climate.fam', sep = ' ', header=None)

# --- Set Up Paths ---
path = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/allele_assoc_runs'
gemma_path = path + '/lmm_gemma/'
bslmm_path = path + '/bslmm/'

# Get list of bioclimatic variables
biovars = [i for i in climate.columns if 'bio' in i]

# Define paths to genotype files
bed_file = path + '/1001g_grenet_climate.bed'
bim_file = path + '/1001g_grenet_climate.bim'

# --- Create Directories and Prepare Files for Each Bioclimatic Variable ---
print("Creating directories and preparing files for each bioclimatic variable...")
for biovar in biovars:
    print(f"Processing {biovar}...")
    
    # Create paths for this bioclimatic variable
    biovar_gemma_path = gemma_path + biovar
    biovar_bslmm_path = bslmm_path + biovar

    # Create directories
    os.makedirs(biovar_gemma_path, exist_ok=True)
    os.makedirs(biovar_bslmm_path, exist_ok=True)

    # Prepare .fam file with bioclimatic variable as phenotype
    biovar1 = climate[['ecotype', biovar]]
    fam = og_fam.merge(biovar1, left_on=0, right_on='ecotype', how='left')
    fam = fam.drop([5, 'ecotype'], axis=1)

    # Save .fam files
    fam.to_csv(biovar_gemma_path + '/1001g_grenet_climate.fam', index=None, header=None, sep=' ')
    fam.to_csv(biovar_bslmm_path + '/1001g_grenet_climate.fam', index=None, header=None, sep=' ')
    
    # Create symbolic links for genotype files
    os.symlink(bed_file, os.path.join(biovar_gemma_path, "1001g_grenet_climate.bed"))
    os.symlink(bim_file, os.path.join(biovar_gemma_path, "1001g_grenet_climate.bim"))
    os.symlink(bed_file, os.path.join(biovar_bslmm_path, "1001g_grenet_climate.bed"))
    os.symlink(bim_file, os.path.join(biovar_bslmm_path, "1001g_grenet_climate.bim"))

# --- Create and Submit GEMMA Jobs ---
print("Creating and submitting GEMMA jobs...")
shfiles = []

kinship_path = path + '/1001g_grenet_climateLDpruned_05maf.cXX.txt'
pcs_path = path + '/pc_1000.txt'

for biovar in biovars:
    biovar_gemma_path = gemma_path + biovar + '/'
    seed = random.randint(1, 100000000)
    file = biovar_gemma_path + f'gemma_{biovar}_1001g.sh'
    
    # Create GEMMA job script
    text = f'''#!/bin/bash
#SBATCH --job-name=gemma_{biovar}_1001g
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=120gb
#SBATCH --cpus-per-task=2
#SBATCH --output=gemma_{biovar}_1001g_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

source /home/tbellagio/miniforge3/etc/profile.d/conda.sh
export PATH="${{PATH}}:/home/username/bin"
cd {biovar_gemma_path}
conda activate /home/tbellagio/miniforge3/envs/gwas

gemma \\
-bfile 1001g_grenet_climate \\
-maf 0.05 \\
-lmm -k {kinship_path} \\
-c {pcs_path} \\
-o {biovar}
'''

    with open(file, 'w') as o:
        o.write(text)
    shfiles.append(file)

# Submit GEMMA jobs
print("Submitting GEMMA jobs...")
for file in shfiles:
    subprocess.run(['sbatch', file])

# --- Create and Submit BSLMM Jobs ---
print("Creating and submitting BSLMM jobs...")
shfiles = []

for biovar in biovars:
    biovar_bslmm_path = bslmm_path + biovar + '/'
    seed = random.randint(1, 100000000)
    file = biovar_bslmm_path + f'bslmm_{biovar}.sh'
    
    # Create BSLMM job script
    text = f'''#!/bin/bash
#SBATCH --job-name=bslmm_{biovar}
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=120gb
#SBATCH --cpus-per-task=2
#SBATCH --output=bslmm_{biovar}_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

source /home/tbellagio/miniforge3/etc/profile.d/conda.sh
cd {biovar_bslmm_path}
conda activate /home/tbellagio/miniforge3/envs/gwas

# Add BSLMM command here
'''

    with open(file, 'w') as o:
        o.write(text)
    shfiles.append(file)

# Submit BSLMM jobs
print("Submitting BSLMM jobs...")
for file in shfiles:
    subprocess.run(['sbatch', file])

print("Script execution completed.")
