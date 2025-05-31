#!/usr/bin/env python
# coding: utf-8

# This script analyzes and visualizes allele frequency changes for specific genomic blocks, particularly the haplotype block containing CAM5
# over time and across different experimental sites, categorized by temperature (hot vs. cold).

# The main goals of the script are:
# 1. Prepare allele frequency data for a manageable subset of SNPs/blocks.
# 2. Plot allele frequency change across experimental sites (spatial plot - output file: adaptive_snps_af_2_1265_final_gen.svg).
# 3. Plot allele frequency change over time, comparing sites categorized as cold vs. warm
#    (temporal plot - output file: colder_Vs_warmer_block_2_1264_nolegend.pdf).

# --- Required Input Files ---
# - ../binomial_regression_lastgen/binomial_reg_results_last_gen.csv: Binomial regression results.
# - ../signficant_intersection/genes_info_bonferronicor_tair10.csv: Information about significant blocks.
# - ../key_files/blocks_snpsid_dict.pkl: Pickle file containing a dictionary mapping blocks to SNP IDs.
# - ../key_files/var_pos_grenenet.csv: SNP position and information.
# - ../baypass_first_gen/merged_hapFIRE_allele_frequency_indexed.csv: Indexed merged allele frequency data.
# - ../key_files/p0_average_seed_mix.csv: Average allele frequencies in the initial seed mix (generation 0).
# - ../key_files/bioclimvars_sites_era5_year_2018.csv: Bioclimatic variables for experimental sites.
# - cyp7_block.csv: Specific allele frequency data file (likely temporary or specific to a block).

# --- Output Files ---
# - af_blocks_to_plot_bonf_corr_all_sign.csv: Temporary CSV with filtered allele frequency data.
# - plot_data_<blocki>.csv: CSV file containing prepared plot data for a specific block.
# - adaptive_snps_af_2_1265_final_gen.svg: SVG plot of allele frequency change across space.
# - warm_vs_cold.pdf: PDF plot of allele frequency change over time in warm vs. cold sites (intermediate plot).
# - only_barcolor.pdf: PDF of the color bar legend (for bio1 gradient).
# - colder_Vs_warmer_block_2_1264_nolegend.pdf: Final PDF plot of allele frequency change over time in warm vs. cold sites.

# --- Import Libraries ---

import pandas as pd
import pickle
import dask.dataframe as dd
import matplotlib.pyplot as plt
import seaborn as sns
import math
import itertools
import subprocess
from scipy.stats import linregress # For linear regression analysis
import matplotlib # Import matplotlib for font settings
import matplotlib.colors as mcolors # For creating custom colormaps
import numpy as np

# In[1]: # Initial data loading and setup


# Load binomial regression results and significant block information.
binom_reg = pd.read_csv('../binomial_regression_lastgen/binomial_reg_results_last_gen.csv')
blocks = pd.read_csv('../signficant_intersection/genes_info_bonferronicor_tair10.csv')

# Load the dictionary mapping blocks to SNP IDs and create a reverse mapping.
dict_blocks_path = '../key_files/blocks_snpsid_dict.pkl'
with open(dict_blocks_path, 'rb') as file:
    dict_blocks = pickle.load(file)
reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}

# Add block information to the binomial regression results.
binom_reg['block'] = binom_reg['snp_id'].map(reverse_mapping)

# Filter binomial regression results to include only the significant blocks.
binom_reg_blocks = binom_reg[binom_reg['block'].isin(blocks['block_id'].unique())]

# Display mean and min slope for filtered blocks (likely for inspection).
# print(binom_reg_blocks.groupby('block')['slope'].mean())
# print(binom_reg_blocks.groupby('block')['slope'].min())

# Example of filtering for a specific block (commented out).
# print(binom_reg_blocks[binom_reg_blocks['block'] == '5_2227'])

# Define the list of blocks to be plotted.
blocks_to_plot = ['1_168', '1_4899', '4_801','4_2115', '5_571', '5_2227']


# In[12]: # Loading SNP dictionary and preparing for allele frequency extraction


# Load SNP dictionary information.
snps_dict = pd.read_csv('../key_files/var_pos_grenenet.csv')

# Add block information to the SNP dictionary.
snps_dict['block'] = snps_dict['id'].map(reverse_mapping)


# In[ ]: # awk command for extracting allele frequencies for selected blocks


# Construct an awk command to filter the large allele frequency file.
# It selects lines corresponding to the SNP IDs within the defined blocks_to_plot.
# The index values are adjusted by +2, likely to account for header lines in the source file.
awk_conditions = []
for block_id in blocks_to_plot:
    # Find the starting and ending line numbers for the SNPs in the current block.
    # Assumes the SNP IDs in snps_dict are in the same order as the rows in the allele frequency file.
    block_indices = snps_dict[snps_dict['block'] == block_id].index.values
    if len(block_indices) > 0:
        start_idx = block_indices[0] + 2
        end_idx = block_indices[-1] + 2
        awk_conditions.append(f"NR >= {start_idx} && NR <= {end_idx}")

# Join the conditions with ' || ' to match the awk syntax for selecting multiple ranges.
awk_filter = ' || '.join(awk_conditions)

# Construct the final awk command.
awk_command = f"awk '{awk_filter}' ../baypass_first_gen/merged_hapFIRE_allele_frequency_indexed.csv > af_blocks_to_plot_bonf_corr_all_sign.csv"

# Run the awk command using subprocess to extract the data.
# This writes the filtered allele frequency data to a temporary CSV file.
try:
    subprocess.run(awk_command, shell=True, check=True, text=True, capture_output=True)
    print("Allele frequency data extracted successfully.")
except subprocess.CalledProcessError as e:
    print(f"Error running awk command: {e}")
    print(f"Stdout: {e.stdout}")
    print(f"Stderr: {e.stderr}")
    # Handle the error appropriately, maybe exit or raise it.


# In[ ]sections below are likely remnants of notebook conversion and are removed.

# In[23]: # Loading and processing filtered allele frequency data

## allele counts

# Load the allele frequency data extracted by the awk command.
 af = pd.read_csv('af_blocks_to_plot_bonf_corr_all_sign.csv',header=None)


# Load the header from the original allele frequency file to assign column names.
usecols = pd.read_csv('../baypass_first_gen/merged_hapFIRE_allele_frequency_indexed.csv',nrows =1).columns 
af.columns = usecols

# Add block information to the allele frequency DataFrame.
af['block'] = af['0'].map(reverse_mapping)


# In[25]: # Further data processing and preparation for plotting


# Load average allele frequencies from the seed mix (generation 0).
p0_average_seed_mix = pd.read_csv('../key_files/p0_average_seed_mix.csv')
# Concatenate with snps_dict, likely to get SNP IDs aligned with initial frequencies.
p0_average_seed_mix = pd.concat([p0_average_seed_mix,snps_dict ], axis=1)

# The following loop appears to repeat the data loading and processing steps for block '4_2115'.
# This might be redundant if the goal is to plot a single block (like 2_1265, as seen later).
# Keeping it commented out as it might have been for testing or a different analysis flow.

# for blocki in ['4_2115']:
#     print(blocki)
#     ## filter the allles
#     aff = af[af['block'].isin([blocki])].drop('block', axis=1).copy()
#     ## filter the snpsid
#     snps_dict_filt = snps_dict[snps_dict['block'].isin([blocki])]
#     snps_dict_filt = snps_dict_filt.reset_index(drop=True)

#     causal = pd.concat([aff, snps_dict_filt],axis=1)

#     causal = causal[causal['maf05filter'].notna()]

#     causal_det = causal[['pos', 'chrom', 'maf05filter', 'total_alleles05filter_firstgen', 'block']].copy()

#     causal =causal.drop(['pos', 'chrom', 'maf05filter', 'total_alleles05filter_firstgen', 'total_alleles05filter_lastgen', 'block','Unnamed: 745', '0'],axis=1).set_index('id')

#     causal = causal.T.reset_index()

#     causal['site'] = causal['index'].str.split('_').str[0].astype(int)
#     causal['gen'] = causal['index'].str.split('_').str[1].astype(int)
#     causal['plot'] = causal['index'].str.split('_').str[2].astype(int)

#     clim_sites_during_exp = pd.read_csv('../key_files/bioclimvars_sites_era5_year_2018.csv')
#     causal = causal.merge(clim_sites_during_exp[['site', 'bio1']]).drop('index',axis=1)#.drop(['index'],axis=1)

#     # List to hold the data for plotting
#     plot_data = {}

#     # Assuming 'cold' is your DataFrame and 'p0_average_seed_mix' is defined
#     for site, group in causal.groupby(['site']):
#         for snp in group.columns[:-4]:
#             to_plot = group[[snp, 'site', 'gen', 'plot']].copy()
#             initial_freq = p0_average_seed_mix[p0_average_seed_mix['id'] == snp]['0'].values[0]
#             row = to_plot.iloc[0, :].copy()
#             row[snp] = initial_freq
#             row['gen'] = 0
#             to_plot = pd.concat([to_plot, row.to_frame().T], ignore_index=True)
#             to_plot.columns = ['freq', 'site', 'gen', 'plot']
#             # Append to the list
#             plot_data[str(site[0]) + '-' + snp] = to_plot

#     plot_data = pd.concat(plot_data,axis=0).reset_index().drop('level_1',axis=1)

#     plot_data['snp'] = plot_data['level_0'].str.split('-').str[1]

#     plot_data = plot_data.drop('level_0',axis=1)

#     plot_data = plot_data.merge(clim_sites_during_exp[['site', 'bio1']])

#     plot_data.to_csv(f'plot_data_{blocki}.csv',index=None)

# Manually setting the block to plot, overriding the blocks_to_plot list for this part.
# This indicates the script is currently focused on plotting this specific block.
blocki = '2_1265'

# Filter allele frequency data and SNP dictionary for the selected block.
aff = af[af['block'].isin([blocki])].drop('block', axis=1).copy()
snps_dict_filt = snps_dict[snps_dict['block'].isin([blocki])]
snps_dict_filt = snps_dict_filt.reset_index(drop=True)

# Combine filtered allele frequency data and SNP details.
causal = pd.concat([aff, snps_dict_filt],axis=1)

# Filter for rows with valid maf05filter (quality control step).
causal = causal[causal['maf05filter'].notna()]

# Separate detailed SNP information and allele frequency data.
causal_det = causal[['pos', 'chrom', 'maf05filter', 'total_alleles05filter_firstgen', 'block']].copy()
# Keep only allele frequency columns and reshape the DataFrame.
causal = causal.drop(['pos', 'chrom', 'maf05filter', 'total_alleles05filter_firstgen', 'total_alleles05filter_lastgen', 'block','Unnamed: 745', '0'],axis=1).set_index('id')
causal = causal.T.reset_index()

# Extract site, generation, and plot information from the index.
causal['site'] = causal['index'].str.split('_').str[0].astype(int)
causal['gen'] = causal['index'].str.split('_').str[1].astype(int)
causal['plot'] = causal['index'].str.split('_').str[2].astype(int)

# Load climate data and merge with the allele frequency data.
clim_sites_during_exp = pd.read_csv('../key_files/bioclimvars_sites_era5_year_2018.csv')
causal = causal.merge(clim_sites_during_exp[['site', 'bio1']]).drop('index',axis=1)

# Prepare data for plotting allele frequency over time for each SNP and site.
plot_data = {}
for site, group in causal.groupby(['site']):
    for snp in group.columns[:-4]: # Iterate over SNP columns (excluding the last 4 added columns)
        to_plot = group[[snp, 'site', 'gen', 'plot']].copy()
        # Get initial frequency from the seed mix for the current SNP.
        initial_freq = p0_average_seed_mix[p0_average_seed_mix['id'] == snp]['0'].values[0]
        # Create a row for generation 0 with the initial frequency.
        row = to_plot.iloc[0, :].copy()
        row[snp] = initial_freq
        row['gen'] = 0
        # Add the generation 0 data to the plot data for this SNP and site.
        to_plot = pd.concat([to_plot, row.to_frame().T], ignore_index=True)
        # Rename columns for clarity.
        to_plot.columns = ['freq', 'site', 'gen', 'plot']
        # Store the data for this SNP and site.
        plot_data[str(site[0]) + '-' + snp] = to_plot

# Combine data for all SNPs and sites into a single DataFrame.
plot_data = pd.concat(plot_data,axis=0).reset_index().drop('level_1',axis=1)

# Extract SNP ID from the concatenated index.
plot_data['snp'] = plot_data['level_0'].str.split('-').str[1]

# Drop the old index level.
plot_data = plot_data.drop('level_0',axis=1)

# Merge with climate data again (potentially redundant if already merged with 'causal').
plot_data = plot_data.merge(clim_sites_during_exp[['site', 'bio1']])

# Save the prepared plot data to a CSV file.
plot_data.to_csv(f'plot_data_{blocki}.csv',index=None)


# In[14]: # Loading saved plot data (if needed)


# Load the prepared plot data from the CSV file.
# This step allows the script to be run starting from the saved CSV if the data preparation part was already done.
plot_data = pd.read_csv(f'plot_data_{blocki}.csv')


# In[45]: # Plotting Allele Frequency Change Over Time (Cold vs. Warm Sites)


# Define custom colors for plotting, likely representing a temperature gradient.
custom_colors = [
    '#b2182b', '#b2182b', '#bc2b34', '#bc2b34', '#c53e3d', '#cf5246', '#cf5246', '#d86551', 
    '#e0775f', '#e0775f', '#e8896d', '#f09c7b', '#f09c7b', '#f5ac8b', '#f8bb9e', '#f8bb9e', 
    '#fac9b0', '#fac9b0', '#fcd7c2', '#fce1d1', '#fce1d1', '#fae8dd', '#f9f0ea', '#f9f0ea', 
    '#f7f7f7', '#f7f7f7', '#d9e9f1', '#abd2e5', '#72b1d3', '#3c8abe', '#2166ac'
]

# Set global parameters for matplotlib plots.
plt.rcParams['axes.axisbelow'] = True
fontsize = 12 # Define font size
plt.rc('font', family='sans-serif', size=fontsize, weight='normal')
plt.rc('axes', titlesize=fontsize, labelsize=fontsize)
plt.rc('xtick', labelsize=fontsize)
plt.rc('ytick', labelsize=fontsize)
dark_grey = '#4D4D4D'  # Define dark grey color

# Filter datasets based on 'bio1' using a specific cutoff.
warm = plot_data[plot_data['bio1'] > 9.541]
cold = plot_data[plot_data['bio1'] <= 9.541]

# Calculate mean allele frequencies per generation per site for both warm and cold sites.
mean_warm = warm.groupby(['gen', 'site'], as_index=False).mean()
mean_cold = cold.groupby(['gen', 'site'], as_index=False).mean()

# Create color palettes for warm (Reds) and cold (Blues_r) sites, sorted by bio1.
warm_color_map = {}
if len(warm['site'].unique()) > 0:
    warm_sites_with_bio1 = warm.groupby('site')['bio1'].mean()
    sorted_warm_sites = warm_sites_with_bio1.sort_values().index
    warm_colors = sns.color_palette("Reds", len(sorted_warm_sites))
    warm_color_map = dict(zip(sorted_warm_sites, warm_colors))

cold_color_map = {}
if len(cold['site'].unique()) > 0:
    cold_sites_with_bio1 = cold.groupby('site')['bio1'].mean()
    sorted_cold_sites = cold_sites_with_bio1.sort_values().index
    cold_colors = sns.color_palette("Blues_r", len(sorted_cold_sites))
    cold_color_map = dict(zip(sorted_cold_sites, cold_colors))

# Create subplots: one for cold sites and one for warm sites.
fig, axes = plt.subplots(ncols=2, figsize=(12, 8), sharex=True)

# Function to apply styling to plot axes (defined multiple times in the original script, using the last version).
def apply_styling(ax):
    # Apply to spines (hide them).
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Set tick parameters and enable grid on x-axis.
    plt.tick_params(axis='both', colors=dark_grey)
    ax.grid(True, color='lightgrey', alpha=0.7, which='both', axis='x')
    ax.set_xticks([0, 1, 2, 3]) # Ensure x-ticks are set explicitly.

# Plot individual data points and lines for cold sites.
sns.scatterplot(data=cold, x='gen', y='freq', hue='site', palette=cold_color_map, ax=axes[0], alpha=0.2, s=40, legend=False, edgecolor='none')
sns.lineplot(data=cold, x='gen', y='freq', hue='site', ax=axes[0], palette=cold_color_map, alpha=0.3,legend=False)
# Plot mean values for cold sites.
for site in mean_cold['site'].unique():
    site_data = mean_cold[mean_cold['site'] == site]
    sns.scatterplot(data=site_data, x='gen', y='freq', ax=axes[0], s=150, color=cold_color_map[site], alpha=0.9,legend=False, edgecolor='none')

# Add regression p-value to the cold sites plot.
# add_regression_p_value(cold, axes[0]) # Commented out as the function was defined but not used in the final plot section in the original script.

# Customize the cold sites subplot.
axes[0].set_title('Cold Gardens')
axes[0].set_ylabel('Allele Frequency')
axes[0].set_xlabel('Generations')
apply_styling(axes[0])

# Plot individual data points and lines for warm sites.
sns.scatterplot(data=warm, x='gen', y='freq', hue='site', palette=warm_color_map, ax=axes[1], alpha=0.2, s=40, legend=False, edgecolor='none')
sns.lineplot(data=warm, x='gen', y='freq', hue='site', ax=axes[1], palette=warm_color_map, alpha=0.3,legend=False)
# Plot mean values for warm sites.
for site in mean_warm['site'].unique():
    site_data = mean_warm[mean_warm['site'] == site]
    sns.scatterplot(data=site_data, x='gen', y='freq', ax=axes[1], s=150, color=warm_color_map[site], alpha=0.9,legend=False, edgecolor='none')

# Add regression p-value to the warm sites plot.
# add_regression_p_value(warm, axes[1]) # Commented out for consistency with the original script's final plot section.

# Customize the warm sites subplot.
axes[1].set_title('Warm Gardens')
axes[1].set_xlabel('Generations')
axes[1].set_ylabel('Allele Frequency')
apply_styling(axes[1])

# Create data for the legend, including site number and average bio1, sorted by bio1.
bio1_data = []
for site in mean_cold['site'].unique():
    avg_bio1 = mean_cold[mean_cold["site"] == site]["bio1"].mean()
    bio1_data.append((site, avg_bio1, cold_color_map[site]))
for site in mean_warm['site'].unique():
    avg_bio1 = mean_warm[mean_warm["site"] == site]["bio1"].mean()
    bio1_data.append((site, avg_bio1, warm_color_map[site]))
bio1_df = pd.DataFrame(bio1_data, columns=['Site', 'Avg_Bio1', 'Color'])
bio1_df.sort_values(by='Avg_Bio1', inplace=True)

# Prepare legend handles and labels based on sorted bio1 values.
handles = []
labels = []
for _, row in bio1_df.iterrows():
    handles.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=row['Color'], markersize=10))
    labels.append(f'Site {int(row["Site"])} bio1: {row["Avg_Bio1"]:.2f}')

# Adjust layout and save the plot to PDF.
plt.tight_layout(rect=[0, 0, 1, 0.9]) # Adjust layout to make space for the legend (if added separately)
plt.savefig('colder_Vs_warmer_block_2_1264_nolegend.pdf', bbox_inches='tight') # Save the plot.
# plt.savefig('colder_Vs_warmer_block_2_1264_nolegend.png', bbox_inches='tight', dpi=900) # Example saving as PNG.

plt.show() # Display the plot.

# --- Linear Regression and Selection Coefficient Calculation (for cold sites) ---

# The following sections perform linear regression of allele frequency on generation
# for each SNP within the cold sites and calculate selection coefficients.

# Example of iterating through cold SNPs and printing regression results (commented out).
# for i in cold['snp'].unique():
#     initial_freq = p0_average_seed_mix[p0_average_seed_mix['id'] == i]['0']
#     print(initial_freq.values[0])
#     data = cold[cold['snp'] == i]
#     slope, intercept, r_value, p_value, std_err = linregress(data['gen'], data['freq'])
#     print(slope, intercept, r_value, p_value, std_err)

# Load initial frequencies and SNP dictionary again (redundant if already loaded).
# p0_average_seed_mix = pd.read_csv('../key_files/p0_average_seed_mix.csv')
# snps_dict = pd.read_csv('../key_files/var_pos_grenenet.csv')
# p0_average_seed_mix = pd.concat([p0_average_seed_mix,snps_dict ], axis=1)

# Store linear regression and selection coefficient results in a list.
results = []

for snp_id in cold['snp'].unique():
    # Get initial allele frequency for the current SNP.
    initial_freq_row = p0_average_seed_mix[p0_average_seed_mix['id'] == snp_id]
    if not initial_freq_row.empty:
        initial_freq = initial_freq_row['0'].values[0]
    else:
        initial_freq = np.nan # Handle cases where initial frequency is not found.

    # Extract data for the specific SNP in cold sites.
    data = cold[cold['snp'] == snp_id].copy() # Use .copy() to avoid SettingWithCopyWarning

    # Add generation 0 with initial frequency for regression.
    if not pd.isna(initial_freq):
        gen0_data = {'gen': 0, 'freq': initial_freq, 'site': np.nan, 'plot': np.nan, 'bio1': np.nan, 'snp': snp_id}
        data = pd.concat([pd.DataFrame([gen0_data]), data], ignore_index=True)
        data = data.sort_values(by='gen').reset_index(drop=True)

    # Perform linear regression on generation vs allele frequency.
    # Ensure there is enough data for regression (at least 2 points).
    if len(data['gen'].unique()) >= 2:
        # Group by generation and calculate mean frequency for regression.
        regression_data = data.groupby('gen')['freq'].mean().reset_index()
        slope, intercept, r_value, p_value, std_err = linregress(regression_data['gen'], regression_data['freq'])

        # Calculate selection coefficient s (s = slope / (p0 * (1 - p0))).
        selection_coefficient = None
        if not pd.isna(initial_freq) and initial_freq > 0 and initial_freq < 1:
             selection_coefficient = slope / (initial_freq * (1 - initial_freq))

        # Append results to the list.
        results.append({
            'snp': snp_id,
            'initial_freq': initial_freq,
            'slope': slope,
            'selection_coefficient': selection_coefficient,
            'r_value': r_value,
            'p_value': p_value,
            'std_err': std_err
        })
    else:
         print(f"Skipping regression for SNP {snp_id} in cold sites: insufficient data points.")

# Convert results to DataFrame for better readability.
s_cold = pd.DataFrame(results)

# The following sections appear to be incomplete or commented out parts of the script
# related to further analysis or plotting, possibly specific to particular blocks
# (like 2_1264 or 2_1265) or comparing different blocks.

# In[44]: # Appears to be a repeated start of the previous section.

# from scipy.stats import linregress

# Store results in a list to display as a DataFrame later
# results = [] # This variable is redefined here.
# (Rest of the code in this section is the same as the previous one).


# In[ ]: # Empty cells or commented out code for other plots or analyses.


# Code related to spatial plotting (across sites at the final generation) might be here.
# This part would generate the 'adaptive_snps_af_2_1265_final_gen.svg' plot.
# It would likely involve:
# 1. Filtering allele frequency data for the final generation.
# 2. Merging with site geographical or climatic data.
# 3. Creating a scatter plot of allele frequency vs. a spatial variable (e.g., longitude, latitude, or bio1).
# 4. Saving the plot as adaptive_snps_af_2_1265_final_gen.svg.

# Example structure for spatial plot (conceptual):
# final_gen_data = plot_data[plot_data['gen'] == 3] # Assuming generation 3 is the final generation.
# # Merge with site location data (need to load site coordinate data).
# # final_gen_data = final_gen_data.merge(site_coordinates, on='site')
# plt.figure(figsize=(10, 6))
# sns.scatterplot(data=final_gen_data, x='bio1', y='freq', hue='site', palette=site_colors)
# plt.xlabel('Mean Annual Temperature (bio1)')
# plt.ylabel('Allele Frequency (Generation 3)')
# plt.title(f'Allele Frequency vs. Bio1 for Block {blocki} (Generation 3)')
# plt.savefig(f'adaptive_snps_af_{blocki}_final_gen.svg', format='svg')
# plt.show()


# Example of generating the color bar legend separately (as seen in the original script).
# This might be used if the legend is manually added to combined plots.
# sites = range(len(custom_colors))
# bio1_values = np.linspace(plot_data['bio1'].min(), plot_data['bio1'].max(), len(sites))
# plot_data_for_colorbar = pd.DataFrame({'site': sites, 'bio1': bio1_values})
# custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", custom_colors, N=256)
# fig_cbar, ax_cbar = plt.subplots(figsize=(6, 1))
# fig_cbar.subplots_adjust(bottom=0.5)
# gradient = np.linspace(0, 1, 256).reshape(1, -1)
# ax_cbar.imshow(gradient, aspect="auto", cmap=custom_cmap)
# ax_cbar.set_axis_off()
# cbar_ax = fig_cbar.add_axes([0.05, 0.2, 0.9, 0.15])
# norm = mcolors.Normalize(vmin=plot_data_for_colorbar['bio1'].min(), vmax=plot_data_for_colorbar['bio1'].max())
# sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
# sm.set_array([])
# cbar = plt.colorbar(sm, cax=cbar_ax, orientation='horizontal')
# cbar.set_label('bio1 Values')
# plt.savefig('only_barcolor.pdf')
# plt.show()






