#!/usr/bin/env python
# coding: utf-8

"""
Prepare .fam files and Submit GEMMA Jobs for Fitness GWAS

This script prepares input files (`.fam`) for the GEMMA software using average ecotype
frequencies per site as phenotypes. It then creates and submits SBATCH job scripts
to run GEMMA for each experimental site.

Inputs:
- ../../final_gen.csv: List of sample names for the final generation.
- ../../delta_ecotype_freq.txt: Delta ecotype frequency data.
- /carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/greneNet_final_v1.1.recode.fam: Original .fam file with sample information.
- /carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/greneNet_final_v1.1.recode.bed: PLINK .bed file (symlinked).
- /carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/greneNet_final_v1.1.recode.bim: PLINK .bim file (symlinked).
- /carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/kin/greneNet_final_v1.1.recode.cXX.txt: Kinship matrix for GEMMA.
- /carnegie/nobackackup/scratch/tbellagio/gea_grene-net/gwas/20cov_gwas: File containing principal components or covariates.

Outputs:
- For each unique site:
  - A directory named `site_<site_id>` containing:
    - greneNet_final_v1.1.recode.fam: Modified .fam file with average ecotype frequency as phenotype.
    - greneNet_final_v1.1.recode.bed: Symbolic link to the main .bed file.
    - greneNet_final_v1.1.recode.bim: Symbolic link to the main .bim file.
    - gemma_site_<site_id>.sh: SBATCH script for running GEMMA for the site.
- Submitted GEMMA jobs on the cluster via sbatch.
"""

# --- Import Libraries ---
import pandas as pd
import os
import subprocess
import random

# --- Load Data ---
# Load final generation sample names
final_samples = pd.read_csv('../../final_gen.csv')['sample_name']

# Load delta ecotype frequency data for final samples
# Assuming the file has ecotypes as index and samples as columns
ef = pd.read_csv('../../delta_ecotype_freq.txt', sep = '\t', usecols = final_samples)

# Load original .fam file
og_fam = pd.read_csv('/carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/greneNet_final_v1.1.recode.fam', sep = ' ', header=None)

# --- Prepare Ecotype Frequency Data for Phenotype ---
# Transpose, reset index, extract site, and calculate mean frequency per site
# This results in a DataFrame with ecotype IDs as index and site IDs as columns,
# where values are the average ecotype frequencies for each site.
ef = ef.T.reset_index()
ef['site'] = ef['index'].str.split('_').str[0]
# Group by site and calculate the mean ecotype frequency across samples for each site
ef = ef.groupby('site')[ef.columns[:-1]].mean()
# Transpose back to have ecotype IDs as columns and site IDs as index (as needed for fam file alignment?)
ef = ef.T

# Get unique site IDs (columns of the processed ef DataFrame)
unique_sites = ef.columns

# --- Configuration for File Paths and GEMMA ---
# Define base path for output directories
output_base_path = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/fitness_gwas/'

# Define paths to input PLINK files and GEMMA related files
bed_file = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/greneNet_final_v1.1.recode.bed'
bim_file = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/greneNet_final_v1.1.recode.bim'
kinship_path = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/kin/greneNet_final_v1.1.recode.cXX.txt'
pcs_path = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/20cov_gwas'

# --- Create .fam files and Symlinks per Site ---
print("Creating site directories, .fam files, and symlinks...")
for site in unique_sites:
    site_path = output_base_path + 'site_' + site

    # Create site-specific directory if it doesn't exist
    os.makedirs(site_path, exist_ok=True)

    # Create the site-specific .fam file
    # The original fam file columns are typically: Family ID, Individual ID, Paternal ID, Maternal ID, Sex, Phenotype
    # We replace the Phenotype column (index 5, 0-based) with the average ecotype frequencies for the current site.
    # Assuming original fam has at least 6 columns.
    # Concatenate original fam (excluding original phenotype) with new phenotype column
    fam_df = pd.concat([og_fam.iloc[:, :5], ef[site]], axis=1)

    # Save the modified .fam file
    fam_df.to_csv(os.path.join(site_path, 'greneNet_final_v1.1.recode.fam'), index=None, header=None, sep = ' ')
    print(f"Created .fam file for site {site} at {site_path}/")
    
    # Create symbolic links for .bed and .bim files in the site directory
    # Use absolute paths for symlinks to avoid issues with changing directories
    bed_symlink_path = os.path.join(site_path, "greneNet_final_v1.1.recode.bed")
    bim_symlink_path = os.path.join(site_path, "greneNet_final_v1.1.recode.bim")
    
    # Remove existing symlinks if they exist before creating new ones
    if os.path.exists(bed_symlink_path):
        os.remove(bed_symlink_path)
    if os.path.exists(bim_symlink_path):
        os.remove(bim_symlink_path)
        
    os.symlink(bed_file, bed_symlink_path)
    os.symlink(bim_file, bim_symlink_path)
    print(f"Created symlinks for site {site} at {site_path}/")

# --- Create and Submit GEMMA Job Scripts per Site ---
print("Creating and submitting GEMMA job scripts...")
shfiles = []

## submit gemma jobs
for site in unique_sites:
    site_path = output_base_path + 'site_' + site + '/'
    seed = random.randint(1, 100000000)
    file_name = f'gemma_site_{site}.sh'
    file_path = os.path.join(site_path, file_name)
    
    # Define SBATCH script content
    # Ensure paths in the script are correct for the execution environment
    text = f'''#!/bin/bash
#SBATCH --job-name=gemma_site_{site}
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=120gb
#SBATCH --cpus-per-task=2
#SBATCH --output=gemma_site_{site}_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

# Load necessary modules or activate conda environment
source /home/tbellagio/miniforge3/etc/profile.d/conda.sh
conda activate /home/tbellagio/miniforge3/envs/gwas

# Change to the site-specific directory
cd {site_path}

# Run GEMMA command
# -bfile: prefix for the genotype files (.bed, .bim, .fam)
# -maf: minor allele frequency threshold
# -lmm: perform linear mixed model analysis
# -k: path to the kinship matrix
# -c: path to the covariate file (PCs)
# -o: output file name prefix
gemma \
-bfile greneNet_final_v1.1.recode \
-maf 0.05 \
-lmm -k {kinship_path} \
-c {pcs_path} \
-o {site}

'''

    # Write the SBATCH script to file
    with open(file_path, 'w') as o:
        o.write(text)
    shfiles.append(file_path)
    print(f"Created SBATCH script: {file_path}")

# Submit SBATCH jobs
print("Submitting SBATCH jobs...")
for i, file in enumerate(shfiles):
    try:
        # Use check=True to raise an exception if the command fails
        subprocess.run(['sbatch', file], check=True)
        print(f"Submitted job {i+1}/{len(shfiles)}: {file}")
    except subprocess.CalledProcessError as e:
        print(f"Error submitting job {i+1}/{len(shfiles)} {file}: {e}")
    except FileNotFoundError:
        print(f"Error: sbatch command not found. Make sure Slurm is in your PATH.")

print("SBATCH job submission process completed.")
