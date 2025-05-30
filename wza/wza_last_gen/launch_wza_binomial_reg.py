#!/usr/bin/env python
# coding: utf-8

"""
Launch WZA Analysis for Binomial Regression Results

This script performs WZA (Weighted-Z Analysis) on binomial regression results.
It processes the regression results, applies genomic inflation factor correction,
and submits the analysis to the cluster.

Required Files:
1. Input Files:
   - ../key_files/var_pos_grenenet.csv: SNP position information
   - ../key_files/blocks_snpsid_dict.pkl: Dictionary mapping SNP IDs to block IDs
   - ../key_files/maf_all_samples_last_gen.csv: Minor allele frequencies
   - ../binomial_regression_lastgen/{biovar}_binomial_reg_results_last_gen.csv: Binomial regression results

2. Script Files:
   - general_WZA_script_mod_polynomial_order7.py: WZA analysis script

Outputs:
1. Analysis Results:
   - binomial_reg_lastgen_wmaf_{biovar}.csv: Processed regression results with MAF
   - wza_binomial_regression_{biovar}_poly7.csv: WZA analysis results
   - wza_{biovar}.sh: SLURM job scripts

Script Outline:
1. Load SNP and block mapping data
2. Process binomial regression results for each bioclimatic variable
3. Apply genomic inflation factor correction
4. Create and submit SLURM jobs for WZA analysis
"""

# --- Import Libraries ---
import pandas as pd
import os 
import random
import subprocess
import dask.dataframe as dd
import seaborn as sns
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from scipy.stats import chi2
import pickle

# --- Define Helper Functions ---
def calculate_genomic_inflation_factor(p_values):
    """
    Calculate the genomic inflation factor (lambda) for p-value correction.
    
    Parameters:
    p_values (array-like): Array of p-values to correct
    
    Returns:
    float: Genomic inflation factor
    """
    chi_squared_stats = -2 * np.log(p_values)
    lambda_gif = np.median(chi_squared_stats) / chi2.ppf(0.5, 1)
    return lambda_gif

def adjust_p_values(p_values, lambda_gif):
    """
    Adjust p-values using the genomic inflation factor.
    
    Parameters:
    p_values (array-like): Array of p-values to correct
    lambda_gif (float): Genomic inflation factor
    
    Returns:
    array: Adjusted p-values
    """
    adjusted_chi_squared = -2 * np.log(p_values) / lambda_gif
    adjusted_p_values = np.exp(-adjusted_chi_squared / 2)
    return adjusted_p_values

# --- Load SNP and Block Mapping ---
print("Loading SNP and block mapping...")
snps_names = pd.read_csv('../key_files/var_pos_grenenet.csv')

with open('../key_files/blocks_snpsid_dict.pkl', 'rb') as f:
    dict_blocks = pickle.load(f)

reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}
snps_names = snps_names[snps_names['total_alleles05filter_firstgen'].notna()].reset_index(drop=True)

# --- Load MAF Data ---
print("Loading MAF data...")
maf = pd.read_csv('../key_files/maf_all_samples_last_gen.csv')

# --- Define Bioclimatic Variables ---
biovars = [
    'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9',
    'bio10', 'bio11', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio19'
]

# --- Process Binomial Regression Results ---
print("Processing binomial regression results...")
for biovar in biovars:
    print(f"Processing {biovar}...")
    binomial_reg = pd.read_csv(f'../binomial_regression_lastgen/{biovar}_binomial_reg_results_last_gen.csv')
    
    # Map SNP IDs to blocks
    binomial_reg['block'] = binomial_reg['snp_id'].map(reverse_mapping)
    
    # Apply genomic inflation factor correction
    p_values = binomial_reg['pvalue']
    lambda_gif = calculate_genomic_inflation_factor(p_values)
    adjusted_p_values = adjust_p_values(p_values, 5)
    binomial_reg['adj_pvalue'] = adjusted_p_values
    
    # Add MAF information
    binomial_reg = pd.concat([binomial_reg, maf], axis=1)
    binomial_reg.columns = ['slope', 'pvalue', 'snp_id', 'block', 'adj_pvalue', 'MAF']
    
    # Save processed results
    binomial_reg.to_csv(
        f'../binomial_regression_lastgen/binomial_reg_lastgen_wmaf_{biovar}.csv',
        index=None
    )

# --- Create and Submit SLURM Jobs ---
print("Creating SLURM job scripts...")
path = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/wza_last_gen/'
shfiles = []

for biovar in biovars:
    seed = random.randint(1, 100000000)
    file = f'wza_{biovar}.sh'
    cmd = f'''python general_WZA_script_mod_polynomial_order7.py \\
        --correlations ../binomial_regression_lastgen/binomial_reg_lastgen_wmaf_{biovar}.csv \\
        --summary_stat pvalue \\
        --window "block" \\
        --output wza_binomial_regression_{biovar}_poly7.csv \\
        --sep ","'''

    text = f'''#!/bin/bash
#SBATCH --job-name=wza_{biovar}
#SBATCH --time=1:00:00  # Time limit set to 4 days
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --output=wza_%j_{biovar}.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

module load python/3.11_conda
conda activate /home/tbellagio/miniforge3/envs/pipeline_snakemake
export LD_LIBRARY_PATH="/home/tbellagio/miniforge3/envs/run_baypass/lib:$LD_LIBRARY_PATH"
cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/wza_last_gen
{cmd}
'''

    with open(file, 'w') as o:
        o.write(text)
    shfiles.append(file)

# Submit jobs
print("Submitting SLURM jobs...")
for shfile in shfiles:
    subprocess.run(["sbatch", shfile], check=True)

print("Script execution completed.")


# In[38]:


biovars[5]


# In[ ]:





# In[ ]:


biovar=biovars[5]


# In[ ]:


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# In[ ]:





# In[ ]:





# In[51]:


sns.set_context("talk")


# In[63]:


biovars = ['bio1', 'bio2',
 'bio3',
 'bio4',
 'bio5',
 'bio6',
 'bio7',
 'bio8',
 'bio9',
 'bio10',
 'bio11',
 'bio12',
 'bio13',
 'bio14',
 'bio15',
 'bio16',
 'bio17',
 'bio18',
 'bio19']


# In[81]:


for biovar in biovars:
    wza = pd.read_csv(f'wza_binomial_regression_{biovar}_poly7.csv')

    wza['chrom'] = wza['gene'].str.split('_').str[0].astype(int)
    wza['pos'] = wza['gene'].str.split('_').str[1].astype(int)
    # Assuming 'wza' and 'biovar' are defined and properly set up as needed for your plots
    # Data preparation and calculations
    observed_quantiles = -np.log10(np.sort(wza['Z_pVal'].values))
    expected_quantiles = -np.log10(np.linspace(1 / len(wza), 1, len(wza)))
    #threshold_value = 0.05 / len(wza)
    # Apply Benjamini-Hochberg correction and find critical p-value
    _, adjusted_pvals, _, _ = multipletests(wza['Z_pVal'], alpha=0.05, method='fdr_bh')
    wza['adjusted_pval'] = adjusted_pvals

    # Find the largest raw p-value that corresponds to an adjusted p-value less than or equal to the FDR threshold
    critical_pvalue = wza.loc[wza['adjusted_pval'] <= 0.05, 'Z_pVal'].max()
    #significance_line = all[all['Bonferroni_corrected_pval'] < 0.05]['-log10(pvalue)'].min()
    # DataFrame setup for the Manhattan plot
    df = wza[['Z_pVal', 'pos', 'chrom']].copy()
    df['chromosome'] = df['chrom']
    df['position'] = df['pos']
    df['-log10(pvalue)'] = -np.log10(df['Z_pVal'])
    
    # Color setup
    colors = sns.color_palette("crest", n_colors=5)
    
    # Calculate chromosome offsets as before
    chromosome_offsets = {}
    offset = 0
    chrom_ends = {}
    for chrom in sorted(wza['chrom'].unique()):
        chromosome_offsets[chrom] = offset
        max_position = wza[wza['chrom'] == chrom]['pos'].max()
        offset += max_position + 200
        chrom_ends[offset] = (chrom, max_position)
    
    # Plotting as before
    # ...
    df['adjusted_position'] = df.apply(lambda row: row['position'] + chromosome_offsets[row['chromosome']], axis=1)
    
    # Create a figure with custom subplot widths
    fig = plt.figure(figsize=(20, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1])  # 80% to 20% width ratio
    
    # Manhattan plot on the first subplot
    ax1 = plt.subplot(gs[0])
    for chrom in sorted(df['chromosome'].unique()):
        subset = df[df['chromosome'] == chrom]
        ax1.scatter(
            subset['adjusted_position'],
            subset['-log10(pvalue)'],
            alpha=0.7,
            c=colors[chrom % len(colors)],
            s=50
        )
    
    ax1.set_xlabel('Chromosomes')
    ax1.set_ylabel('-log10(pvalue)')
    ax1.set_title(f'{biovar} Binomial GLM + WZA')

    ax1.axhline(y=-np.log10(critical_pvalue), color='grey', linestyle='dashed')

    #ax1.axhline(y=-np.log10(threshold_value), color='grey', linestyle='dashed')
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.grid(axis='y')
    
    # QQ plot on the second subplot
    ax2 = plt.subplot(gs[1])
    sns.scatterplot(x=expected_quantiles, y=observed_quantiles, edgecolor='b', facecolor='none', alpha=0.5, ax=ax2)
    ax2.plot([min(expected_quantiles), max(expected_quantiles)], [min(expected_quantiles), max(expected_quantiles)], 'r--')
    ax2.set_xlabel("Expected -log10(p-values)")
    ax2.set_ylabel("Observed -log10(p-values)")
    ax2.set_title(f'QQ Plot for {biovar} Binomial GLM + WZA')
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.grid(axis='y')
        
    # Custom function to determine actual genomic positions from adjusted positions
    def get_original_position(adjusted_pos):
        for end in sorted(chrom_ends.keys()):
            if adjusted_pos <= end:
                chrom, max_pos = chrom_ends[end]
                return f"{chrom}"
        return ""
    
    # Setting the ticks on the Manhattan plot to show actual genomic positions
    ax1.set_xticks([chromosome_offsets[chrom] + wza[wza['chrom'] == chrom]['pos'].max()/2 for chrom in sorted(wza['chrom'].unique())])  # Set ticks at the middle of each chromosome segment
    ax1.set_xticklabels([get_original_position(chromosome_offsets[chrom] + wza[wza['chrom'] == chrom]['pos'].max()/2) for chrom in sorted(wza['chrom'].unique())],)  # Use the function to get original positions as labels
    
    plt.tight_layout()
    plt.savefig(f'last_gen_{biovar}_GLMbinomWZA.png')
    plt.show()


# In[ ]:





# In[ ]:





# In[ ]:


#genes = ['2_199', '3_2730', '5_2244']
genes = ['2_1265']
# Add circles around the specific genes
for gene in genes:
    chrom, pos = map(int, gene.split('_'))
    subset_gene = df[(df['chromosome'] == chrom) & (df['position'] == pos)]
    if not subset_gene.empty:
        plt.scatter(
            subset_gene['adjusted_position'],
            subset_gene['-log10(pvalue)'],
            edgecolor='red',  # Color of the edge of the circle
            linewidth=2,      # Width of the edge line
            facecolor='none', # Facecolor of the circle (None means transparent)
            s=100,            # Size of the circle
            label=f'Gene {gene}'
        )


# In[ ]:





# In[ ]:


#genes = ['2_199', '3_2730', '5_2244']
genes = ['2_1265', '4_801']
# Add circles around the specific genes
for gene in genes:
    chrom, pos = map(int, gene.split('_'))
    subset_gene = df[(df['chromosome'] == chrom) & (df['position'] == pos)]
    if not subset_gene.empty:
        plt.scatter(
            subset_gene['adjusted_position'],
            subset_gene['-log10(pvalue)'],
            edgecolor='red',  # Color of the edge of the circle
            linewidth=2,      # Width of the edge line
            facecolor='none', # Facecolor of the circle (None means transparent)
            s=100,            # Size of the circle
            label=f'Gene {gene}'
        )


# In[54]:


threshold_value = 0.05 / len(wza)
biovar='bio1'
#sm.qqplot(pvalues['pvalue'], line ='45') 
#py.show() 

df = wza[['Z_pVal', 'index']].copy()

colors = sns.color_palette("crest", n_colors = 5)

# Parsing chromosome number and position
df['chromosome'] = 1
df['position'] = df['index']
df['-log10(pvalue)'] = -np.log10(df['Z_pVal'])

# Calculate the offset for each chromosome to prevent overlap
chromosome_offsets = {}
offset = 0
for chrom in sorted(df['chromosome'].unique()):
    chromosome_offsets[chrom] = offset
    max_position = df[df['chromosome'] == chrom]['position'].max()
    offset += max_position + 1000000  # Adding 1 million as a buffer between chromosomes

# Apply offsets to positions
df['adjusted_position'] = df.apply(lambda row: row['position'] + chromosome_offsets[row['chromosome']], axis=1)

# Creating the Manhattan plot
plt.figure(figsize=(20, 6))

for chrom in sorted(df['chromosome'].unique()):
    subset = df[df['chromosome'] == chrom]
    plt.scatter(subset['adjusted_position'], subset['-log10(pvalue)'], c=colors[chrom % len(colors)], label=f'Chr {chrom}', s=10)

# Aesthetics
plt.xlabel('Adjusted Position')
plt.ylabel('-log10(pvalue)')
#plt.title('Manhattan Plot')
#plt.grid(axis='y')
#plt.legend(title="Chromosome", bbox_to_anchor=(1.05, 1), loc='upper left')
ax = plt.gca()  # Get current axes
ax.spines['top'].set_visible(False)  # Remove the top spine
ax.spines['right'].set_visible(False)
# Threshold line (optional)
threshold = -np.log10(threshold_value)
plt.axhline(y=threshold, color='grey', linestyle='dashed')
plt.title(f'{biovar} WZA')  # Set the title

# Show the plot
plt.tight_layout()
#plt.savefig(f'manhattan_{biovar}.png')
plt.show()


# In[ ]:





# In[50]:


threshold_value = 0.05 / len(wza)
biovar='bio1'
#sm.qqplot(pvalues['pvalue'], line ='45') 
#py.show() 

df = wza[['top_candidate_p', 'index']].copy()

colors = sns.color_palette("crest", n_colors = 5)

# Parsing chromosome number and position
df['chromosome'] = 1
df['position'] = df['index']
df['-log10(pvalue)'] = -np.log10(df['top_candidate_p'])

# Calculate the offset for each chromosome to prevent overlap
chromosome_offsets = {}
offset = 0
for chrom in sorted(df['chromosome'].unique()):
    chromosome_offsets[chrom] = offset
    max_position = df[df['chromosome'] == chrom]['position'].max()
    offset += max_position + 1000000  # Adding 1 million as a buffer between chromosomes

# Apply offsets to positions
df['adjusted_position'] = df.apply(lambda row: row['position'] + chromosome_offsets[row['chromosome']], axis=1)

# Creating the Manhattan plot
plt.figure(figsize=(20, 6))

for chrom in sorted(df['chromosome'].unique()):
    subset = df[df['chromosome'] == chrom]
    plt.scatter(subset['adjusted_position'], subset['-log10(pvalue)'], c=colors[chrom % len(colors)], label=f'Chr {chrom}', s=10)

# Aesthetics
plt.xlabel('Adjusted Position')
plt.ylabel('-log10(pvalue)')
#plt.title('Manhattan Plot')
#plt.grid(axis='y')
#plt.legend(title="Chromosome", bbox_to_anchor=(1.05, 1), loc='upper left')
ax = plt.gca()  # Get current axes
ax.spines['top'].set_visible(False)  # Remove the top spine
ax.spines['right'].set_visible(False)
# Threshold line (optional)
threshold = -np.log10(threshold_value)
plt.axhline(y=threshold, color='grey', linestyle='dashed')
plt.title(f'{biovar} WZA')  # Set the title

# Show the plot
plt.tight_layout()
plt.savefig(f'manhattan_{biovar}.png')
plt.show()


# In[ ]:




