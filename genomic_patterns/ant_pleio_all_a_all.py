#!/usr/bin/env python
# coding: utf-8

# This script analyzes the relationship between allele frequency change (represented by the slope from binomial regression) and WZA Z-values across different experimental sites.
# It generates hexbin plots comparing slopes between pairs of sites and saves them as PNG files.

# Import necessary libraries
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns # Imported but not explicitly used for plotting in the final loop, though might be used indirectly by matplotlib styles
from itertools import combinations

# Set matplotlib font type for PDF output (although PNGs are saved, this setting might be a remnant or for other potential outputs)
matplotlib.rcParams['pdf.fonttype'] = 42

# --- Required Input Files ---
# Ensure these files are present relative to the script's location:
# - ../../grene/data/bioclimvars_experimental_sites_era5.csv
# - Files in the 'binom_reg/' directory with the pattern 'wza_site_*_pr.csv'

# --- Data Loading and Preparation ---

# Read bioclimatic variables for experimental sites
# This file contains climate information, including bio1 (Annual Mean Temperature), for different sites.
bios = pd.read_csv('../../grene/data/bioclimvars_experimental_sites_era5.csv')[['site', 'bio1']]

# List all files in the 'binom_reg/' directory
# This directory is expected to contain WZA and binomial regression results per site.
list_files = os.listdir(f'binom_reg/')

# Filter the list of files to get those with 'pr' in their name (likely indicating processed files)
final_all_sites = [i for i in list_files if 'pr' in i]

# Create a DataFrame mapping filenames to site numbers
# This extracts the site number from the filename assuming a specific format (e.g., 'wza_site_6_pr.csv').
sites_names = pd.DataFrame({
    'file': pd.Series(final_all_sites),  # Original filenames
    'site': pd.Series(final_all_sites).str.split('_').str[2].astype(int)  # Extract and convert site number
})

# Merge site names with bioclimatic data to associate sites with climate variables
sites_names = sites_names.merge(bios)


# --- Functions for Data Processing and Plotting ---

# Function to read and preprocess data for a given site file
def prepare_site_data(filename):
    """Read and preprocess data for a given file."""
    df = pd.read_csv(f'binom_reg/{filename}')
    site_number = filename.split('_')[2]  # Extract the site number from the filename
    # Select relevant columns and rename them to include the site number
    df = df[['block', 'Z_pVal', 'slope', 'snp_origin_bio1']]
    df.columns = ['block', f'Z_pVal_site{site_number}', f'slope_site{site_number}', 'snp_origin_bio1']
    return df, site_number

# Function to create and save comparison figures for pairs of sites
def create_comparison_figures(site_files):
    # This loop iterates through each site file as the primary site for comparison
    for primary_file in site_files:
        primary_data, primary_site = prepare_site_data(primary_file)

        # Create a figure with subplots for comparing the primary site to all other sites
        # The grid size (5x6) should match the number of comparison sites + potentially empty plots.
        # Assuming there are 26 site files in total, comparing one site to the other 25 requires 25 subplots. A 5x6 grid has 30 plots.
        fig, axs = plt.subplots(5, 6, figsize=(30, 25))
        axs = axs.flatten()  # Flatten the 2D array of axes for easier indexing

        idx = 0 # Initialize subplot index
        # This inner loop iterates through all site files for comparison with the primary site
        for compare_file in site_files:
            # Avoid comparing a site to itself
            if primary_file != compare_file:
                compare_data, compare_site = prepare_site_data(compare_file)
                # Merge data from the primary site and the comparison site based on 'block'
                merged = primary_data.merge(compare_data, on='block', how='inner')

                # Select the current subplot axes
                ax = axs[idx]
                # Create a hexbin plot comparing the slopes of the two sites
                hb = ax.hexbin(
                    data=merged, x=f'slope_site{primary_site}', y=f'slope_site{compare_site}',
                    gridsize=50, cmap='Greys_r', bins='log', mincnt=1
                )
                # Add horizontal and vertical lines at zero for reference
                ax.axhline(0, color='grey', linewidth=2, zorder=5)
                ax.axvline(0, color='grey', linewidth=2, zorder=5)
                # Set labels for the subplot axes
                ax.set_xlabel(f'Slope Site {primary_site}')
                ax.set_ylabel(f'Slope Site {compare_site}')
                idx += 1 # Move to the next subplot index

        # Remove any unused subplots (if the number of comparisons is less than the grid size)
        for ax in axs[idx:]:
            ax.remove()

        # Adjust layout to prevent overlap of subplots
        fig.tight_layout()
        # Save the figure to a PNG file named after the primary site
        # These are the primary output files generated by this script.
        plt.savefig(f'Comparison_Site_{primary_site}.png')
        plt.close(fig)  # Close the figure to free up memory

# --- Script Execution ---

# List of site files to process for comparison figures
site_files = [
    'wza_site_6_pr.csv', 'wza_site_32_pr.csv', 'wza_site_24_pr.csv', 'wza_site_48_pr.csv',
    'wza_site_13_pr.csv', 'wza_site_53_pr.csv', 'wza_site_45_pr.csv', 'wza_site_4_pr.csv',
    'wza_site_9_pr.csv', 'wza_site_11_pr.csv', 'wza_site_2_pr.csv', 'wza_site_57_pr.csv',
    'wza_site_43_pr.csv', 'wza_site_55_pr.csv', 'wza_site_54_pr.csv', 'wza_site_42_pr.csv',
    'wza_site_1_pr.csv', 'wza_site_49_pr.csv', 'wza_site_52_pr.csv', 'wza_site_12_pr.csv',
    'wza_site_25_pr.csv', 'wza_site_28_pr.csv', 'wza_site_46_pr.csv', 'wza_site_10_pr.csv',
    'wza_site_27_pr.csv', 'wza_site_5_pr.csv'
]

# Create and save the comparison figures
create_comparison_figures(site_files)

# Note: There were commented-out sections in the original notebook for generating
# individual site pair hexbin plots and calculating Pearson correlation.
# These have been removed in this streamlined script to focus on the main output (comparison PNGs).
# If you need those analyses or plots, please refer to the original notebook or ask to add them back.

