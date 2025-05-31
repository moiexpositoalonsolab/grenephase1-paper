#!/usr/bin/env python
# coding: utf-8

# Required input files:
# - '../key_files/generation_1_sample_names.txt': Contains sample names for the first generation.
# - '../key_files/bioclimvars_sites_era5_year_2018.csv': Contains bioclimatic variables for experimental sites.
# - '../key_files/allele_freq_maf05_mincount05_firstgensamples.csv': Contains allele frequencies for the first generation samples.
# - '../key_files/delta_p_maf05_mincount05_firstgensamples.csv': Contains delta p values for the first generation samples.
# - '../key_files/var_pos_grenenet.csv': Contains SNP names and positions.
# - '../key_files/blocks_snpsid_dict.pkl': Contains a dictionary mapping blocks to SNP IDs.

# High-level outline:
# This script calculates Kendall's tau correlation between bioclimatic variables and allele frequencies for the first generation samples.
# It reads sample names, bioclimatic variables, allele frequencies, and delta p values from the required input files.
# It then calculates the mean allele frequency and minor allele frequency (MAF).
# It iterates through each row of allele frequencies and calculates Kendall's tau correlation with the bioclimatic variable.
# The results are stored in a DataFrame and saved to a CSV file.

# In[21]:

import numpy as np
import pandas as pd
import gzip
import argparse
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns

# In[22]:

get_ipython().system('pwd')

# In[23]:

biovar = 'bio18'

# In[24]:

# Read sample names for the first generation
first_gen_samples = pd.read_csv('../key_files/generation_1_sample_names.txt', header=None)[0]
samples = first_gen_samples.to_list()

# Read bioclimatic variables for experimental sites
clim_sites_during_exp = pd.read_csv('../key_files/bioclimvars_sites_era5_year_2018.csv')

# Extract site IDs from sample names
sites_af = pd.Series(samples).str.split('_').str[0].astype(int)
sites_af.name = 'site'

# Merge site IDs with bioclimatic variables
env = sites_af.reset_index().merge(clim_sites_during_exp).drop(['index'], axis=1)

# In[25]:

bio = env[biovar]

# In[26]:

# Read allele frequencies
af = pd.read_csv('../key_files/allele_freq_maf05_mincount05_firstgensamples.csv')

# In[27]:

# Calculate mean allele frequency
p_bar = af.mean(axis=1)

# Calculate minor allele frequency (MAF)
maf = p_bar.apply(lambda x: 1 - x if x > 0.5 else x)

# In[28]:

# maf.to_csv('../key_files/maf_all_sample_first_gen.csv', index=None)

# In[29]:

# Read delta p values
deltap = pd.read_csv('../key_files/delta_p_maf05_mincount05_firstgensamples.csv')

# In[30]:

deltap.head()

# In[31]:

kendall = {}

# Iterate through each row in af and calculate Kendall's tau correlation with bioclimatic variable
for index, row in af.iterrows():
    geno_k_tau, geno_k_tau_p_value = scipy.stats.kendalltau(bio, row)
    kendall[index] = [geno_k_tau, geno_k_tau_p_value]

# In[32]:

kendall = pd.DataFrame(kendall).T

# In[33]:

kendall['MAF'] = maf

# In[34]:

kendall.columns = ["K_tau", "K_tau_p", "MAF"]

# In[35]:

kendall

# In[36]:

kendall

# In[37]:

# Read SNP names and positions
snps_names = pd.read_csv('../key_files/var_pos_grenenet.csv')

# In[38]:

# Read SNP names and positions
snps_names = pd.read_csv('../key_files/var_pos_grenenet.csv')

# Load dictionary mapping blocks to SNP IDs
import pickle
with open('../key_files/blocks_snpsid_dict.pkl', 'rb') as f:
    dict_blocks = pickle.load(f)

# Create reverse mapping from SNP IDs to blocks
reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}

# Filter SNP names to include only those with valid total alleles
snps_names = snps_names[snps_names['total_alleles05filter_firstgen'].notna()].reset_index(drop=True)

# In[39]:

# Concatenate SNP names and positions with Kendall's tau results
kendall = pd.concat([snps_names[['id', 'pos', 'chrom']], kendall], axis=1)

# Map SNP IDs to blocks
kendall['block'] = kendall['id'].map(reverse_mapping)

# In[40]:

# kendall.to_csv(f'kendall_corr_{biovar}.csv', index=None)

# In[53]:

kendall = pd.read_csv(f'kendall_corr_{biovar}.csv')

# In[54]:

kendall

# In[57]:

kendall = kendall.sort_values('K_tau_p').head(100)

# In[58]:

import pickle
dict_blocks = '../key_files/blocks_snpsid_dict.pkl'

with open(dict_blocks, 'rb') as file:
    dict_blocks = pickle.load(file)

reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}

# In[59]:

kendall['block'] = kendall['id'].map(reverse_mapping)

# In[60]:

kendall = kendall.drop_duplicates('block')

# In[61]:

biovar = 'bio18'

# In[62]:

kendall

# In[63]:

kendall.to_csv(f'top_hits_binom_first_gen_{biovar}.csv')

# In[64]:

kendall

# In[53]:

sns.histplot(kendall['K_tau_p'])

# In[55]:

biovar='bio1'

# In[25]:

# Sort the p-values in ascending order
observed_quantiles = -np.log10(np.sort(kendall['K_tau_p'].values))

# Expected quantiles from the uniform distribution
expected_quantiles = -np.log10(np.linspace(1 / len(kendall), 1, len(kendall)))

# QQ plot
sns.scatterplot(x = expected_quantiles, y = observed_quantiles, edgecolor='b', facecolor='none', alpha=0.5)
plt.plot([min(expected_quantiles), max(expected_quantiles)], [min(expected_quantiles), max(expected_quantiles)], 'r--')

plt.xlabel("Expected -log10(p-values)")
plt.ylabel("Observed -log10(p-values)")
plt.title(f'QQ Plot for {biovar} Kendall tau corr')

plt.show()

# In[55]:

threshold_value = 0.05 / len(kendall)

#sm.qqplot(pvalues['pvalue'], line ='45') 
#py.show() 

df = kendall[['K_tau_p', 'pos', 'chrom']].copy()

colors = sns.color_palette("crest", n_colors = 5)

# Parsing chromosome number and position
df['chromosome'] = df['chrom']
df['position'] = df['pos']
df['-log10(pvalue)'] = -np.log10(df['K_tau_p'])

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
plt.title(f'{biovar} Kendall tau correlation')  # Set the title

# Show the plot
plt.tight_layout()
plt.savefig(f'manhattan_kendall_{biovar}.png')
plt.show()
