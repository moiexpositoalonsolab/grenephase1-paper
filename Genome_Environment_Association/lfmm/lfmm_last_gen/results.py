#!/usr/bin/env python3

"""
LFMM Analysis Results Processor (Last Generation)

This script processes and visualizes the results from the LFMM analysis for the last generation data.
It creates Manhattan plots and QQ plots for each bioclimatic variable, identifies significant SNPs,
and saves the results to appropriate directories.

Required Input Files:
1. lfmm_fullresults_all_k/w_calibration_pvalue_full_genome_*.csv
   - Calibrated p-values for each bioclimatic variable
2. lfmm_fullresults_all_k/effect_sizes_simple_full_genome_*.csv
   - Effect sizes for each bioclimatic variable

Output Files:
1. lfmm_fullresults_all_k/manhattan_plots/*.pdf
   - Manhattan plots for each bioclimatic variable
2. lfmm_fullresults_all_k/qq_plots/*.pdf
   - QQ plots for each bioclimatic variable
3. lfmm_fullresults_all_k/significant_snps/*.csv
   - Lists of significant SNPs for each bioclimatic variable
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

# --- Configuration ---
# Base directories
BASE_DIR = '.'
RESULTS_DIR = os.path.join(BASE_DIR, 'lfmm_fullresults_all_k')

# Create output directories if they don't exist
for subdir in ['manhattan_plots', 'qq_plots', 'significant_snps']:
    os.makedirs(os.path.join(RESULTS_DIR, subdir), exist_ok=True)

# List of bioclimatic variables
BIOVARS = ['bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10',
           'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19']

# Significance threshold
SIG_THRESHOLD = 0.05

def create_manhattan_plot(p_values, biovar, k):
    """Create a Manhattan plot for a given bioclimatic variable."""
    # Create figure
    plt.figure(figsize=(15, 5))
    
    # Plot p-values
    plt.scatter(range(len(p_values)), -np.log10(p_values), alpha=0.6)
    
    # Add significance threshold line
    plt.axhline(y=-np.log10(SIG_THRESHOLD), color='r', linestyle='--')
    
    # Customize plot
    plt.title(f'Manhattan Plot - {biovar} (K={k})')
    plt.xlabel('SNP Position')
    plt.ylabel('-log10(p-value)')
    
    # Save plot
    plt.savefig(os.path.join(RESULTS_DIR, 'manhattan_plots', f'manhattan_{biovar}_k{k}.pdf'))
    plt.close()

def create_qq_plot(p_values, biovar, k):
    """Create a QQ plot for a given bioclimatic variable."""
    # Create figure
    plt.figure(figsize=(10, 10))
    
    # Calculate expected p-values
    expected = -np.log10(np.linspace(0, 1, len(p_values) + 1)[1:])
    observed = -np.log10(np.sort(p_values))
    
    # Plot
    plt.scatter(expected, observed, alpha=0.6)
    plt.plot([0, max(expected)], [0, max(expected)], 'r--')
    
    # Customize plot
    plt.title(f'QQ Plot - {biovar} (K={k})')
    plt.xlabel('Expected -log10(p-value)')
    plt.ylabel('Observed -log10(p-value)')
    
    # Save plot
    plt.savefig(os.path.join(RESULTS_DIR, 'qq_plots', f'qq_{biovar}_k{k}.pdf'))
    plt.close()

def find_significant_snps(p_values, effect_sizes, biovar, k):
    """Find significant SNPs and save to CSV."""
    # Find significant SNPs
    sig_mask = p_values < SIG_THRESHOLD
    sig_snps = pd.DataFrame({
        'SNP': range(len(p_values)),
        'p_value': p_values,
        'effect_size': effect_sizes
    }).loc[sig_mask]
    
    # Save to CSV
    sig_snps.to_csv(
        os.path.join(RESULTS_DIR, 'significant_snps', f'significant_snps_{biovar}_k{k}.csv'),
        index=False
    )
    
    return len(sig_snps)

def process_biovar(biovar, k):
    """Process results for a single bioclimatic variable."""
    print(f"Processing {biovar}...")
    
    # Load p-values and effect sizes
    pvalue_file = os.path.join(RESULTS_DIR, f'w_calibration_pvalue_full_genome_{biovar}_k{k}.csv')
    effect_file = os.path.join(RESULTS_DIR, f'effect_sizes_simple_full_genome_{biovar}_k{k}.csv')
    
    p_values = pd.read_csv(pvalue_file).values.flatten()
    effect_sizes = pd.read_csv(effect_file).values.flatten()
    
    # Create plots
    create_manhattan_plot(p_values, biovar, k)
    create_qq_plot(p_values, biovar, k)
    
    # Find significant SNPs
    n_sig = find_significant_snps(p_values, effect_sizes, biovar, k)
    print(f"Found {n_sig} significant SNPs for {biovar}")

def main():
    """Main function to process results for all bioclimatic variables."""
    # Get optimal K value from variance explained results
    var_explained_file = os.path.join(RESULTS_DIR, 'variance_explained_k.csv')
    df = pd.read_csv(var_explained_file)
    optimal_k = df.loc[df['Variance_Explained'].diff().abs() < 0.01, 'K'].iloc[0]
    
    print(f"Using optimal K value: {optimal_k}")
    
    # Process each bioclimatic variable
    for biovar in BIOVARS:
        process_biovar(biovar, optimal_k)

if __name__ == "__main__":
    main()

