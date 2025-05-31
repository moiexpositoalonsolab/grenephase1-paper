#!/usr/bin/env python3

"""
LFMM Analysis Results Processing

This script processes and visualizes the results from the LFMM analysis for multiple bioclimatic variables.
It reads the output files from run_lfmm_full_genome_ridge_multiplebiovars.r and creates various plots
and summary statistics.

Required Input Files:
1. /carnegie/nobackup/scratch/tbellagio/gea_grene-net/lfmm_full/lfmm_fullresults_all_k/w_calibration_pvalue_full_genome_*.csv
   - Calibrated p-values for each bioclimatic variable
2. /carnegie/nobackup/scratch/tbellagio/gea_grene-net/lfmm_full/lfmm_fullresults_all_k/effect_sizes_simple_full_genome_*.csv
   - Effect sizes for each bioclimatic variable

Output Files:
1. Various plots and summary statistics in the results directory
2. Combined results in a single CSV file
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# --- Configuration ---
BASE_DIR = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/lfmm_full'
RESULTS_DIR = os.path.join(BASE_DIR, 'lfmm_fullresults_all_k')

# List of bioclimatic variables
BIOVARS = ['bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10',
           'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19']

def load_results(biovar):
    """Load results for a specific bioclimatic variable."""
    # Load p-values
    pvalue_file = os.path.join(RESULTS_DIR, f'w_calibration_pvalue_full_genome_{biovar}.csv')
    pvalues = pd.read_csv(pvalue_file)
    
    # Load effect sizes
    effect_file = os.path.join(RESULTS_DIR, f'effect_sizes_simple_full_genome_{biovar}.csv')
    effects = pd.read_csv(effect_file)
    
    return pvalues, effects

def create_qq_plot(pvalues, biovar, output_dir):
    """Create QQ plot for p-values."""
    plt.figure(figsize=(10, 10))
    
    # Calculate expected quantiles
    expected = -np.log10(np.linspace(0, 1, len(pvalues) + 1)[1:])
    observed = -np.log10(np.sort(pvalues))
    
    # Create QQ plot
    plt.scatter(expected, observed, alpha=0.5)
    plt.plot([0, max(expected)], [0, max(expected)], 'r--')
    
    plt.xlabel('Expected -log10(p-value)')
    plt.ylabel('Observed -log10(p-value)')
    plt.title(f'QQ Plot for {biovar}')
    
    # Save plot
    plt.savefig(os.path.join(output_dir, f'qq_plot_{biovar}.pdf'))
    plt.close()

def create_manhattan_plot(pvalues, effects, biovar, output_dir):
    """Create Manhattan plot showing p-values and effect sizes."""
    plt.figure(figsize=(15, 5))
    
    # Plot p-values
    plt.scatter(range(len(pvalues)), -np.log10(pvalues), 
               c=effects, cmap='RdBu', alpha=0.6)
    
    # Add significance threshold
    threshold = -np.log10(0.05 / len(pvalues))  # Bonferroni correction
    plt.axhline(y=threshold, color='r', linestyle='--')
    
    plt.xlabel('SNP Index')
    plt.ylabel('-log10(p-value)')
    plt.title(f'Manhattan Plot for {biovar}')
    plt.colorbar(label='Effect Size')
    
    # Save plot
    plt.savefig(os.path.join(output_dir, f'manhattan_plot_{biovar}.pdf'))
    plt.close()

def process_all_results():
    """Process results for all bioclimatic variables."""
    # Create output directory for plots
    plots_dir = os.path.join(RESULTS_DIR, 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    
    # Process each bioclimatic variable
    for biovar in BIOVARS:
        print(f"Processing {biovar}...")
        
        # Load results
        pvalues, effects = load_results(biovar)
        
        # Create plots
        create_qq_plot(pvalues, biovar, plots_dir)
        create_manhattan_plot(pvalues, effects, biovar, plots_dir)
        
        # Calculate summary statistics
        significant = pvalues < (0.05 / len(pvalues))  # Bonferroni correction
        n_significant = significant.sum()
        
        print(f"{biovar}: {n_significant} significant SNPs")

if __name__ == "__main__":
    process_all_results()
