#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import os 
import random
import subprocess

import dask.dataframe as dd
import seaborn as sns
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt

import pandas as pd
import scipy.stats
import numpy as np
import sys, argparse
from scipy.stats import norm
from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype


# In[2]:


maf = pd.read_csv('../key_files/maf_all_samples.csv')


# In[3]:


maf.columns = ['MAF']


# In[4]:


wd = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/lfmm_full/lfmm_fullresults_all_k/'


# In[5]:


pwd - P


# In[6]:


biovar = 'bio1'


# In[17]:


pvalues_file = f'/carnegie/nobackup/scratch/tbellagio/gea_grene-net/lfmm_full/lfmm_fullresults_all_k/lfmm_{biovar}_k25_results.csv'


# In[ ]:





# In[ ]:





# In[ ]:





# In[279]:


lfmm = pd.read_csv(pvalues_file)


# In[20]:


lfmm[lfmm['p_value'].isna()]


# In[22]:


lfmm.dropna()


# In[663]:


lfmm[lfmm['block'] == '5_2922']


# In[287]:


lfmm.groupby('block').size().sort_values().tail(21)


# In[291]:


lfmm[lfmm['block'].isin(nan_blocks_0_min)].groupby('block').size().sort_values()


# In[77]:


lfmm[lfmm['block'].isin(blocks_na)].groupby('block').size()


# In[ ]:





# In[ ]:





# In[23]:


lfmm = pd.concat([lfmm, maf['MAF']],axis=1)


# In[25]:


lfmm.to_csv('lfmm_results_all_samples_k25_wmaf.csv',index=None)


# In[25]:


lfmm.dropna()


# In[49]:


blocks_na


# In[46]:


filter = lfmm[lfmm['MAF'] >= 0.05]


# In[48]:


filter[filter['block'].isin(blocks_na)]


# In[ ]:





# In[59]:


path = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/wza'


# In[209]:


## create a dir 

#        --sample_snps 1729 \
#        --resample 50 \
# create sbatch files to submit on cedar server
shfiles = []

seed = random.randint(1,100000000)
file = 'wza.sh'
cmd = f'python general_WZA_script_mod.py \
        --correlations lfmm_results_all_samples_k25_wmaf.csv \
        --summary_stat p_value \
        --sample_snps 0 \
        --window block \
        --output wza_results_lfmm_bio1.csv \
        --sep ","'
text = f'''#!/bin/bash
#SBATCH --job-name=wza
#SBATCH --time=1:00:00  # Time limit set to 4 days
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --output=wza_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

module load python/3.11_conda
conda activate /home/tbellagio/miniforge3/envs/pipeline_snakemake
export LD_LIBRARY_PATH="/home/tbellagio/miniforge3/envs/run_baypass/lib:$LD_LIBRARY_PATH"
cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/wza
{cmd}


'''
with open(file, 'w') as o:
    o.write("%s" % text)
shfiles.append(file)


# In[210]:


## now run the shfiles
for shfile in shfiles:
    # Submit each sbatch script to the SLURM scheduler
    subprocess.run(["sbatch", shfile], check=True)


# In[57]:


biovar = 'bio1'


# In[ ]:





# In[280]:


wza_df = pd.read_csv('before_filtering_wza_df.csv').drop('Unnamed: 0', axis=1)


# In[ ]:





# In[281]:


wza_t.shape


# In[282]:


# remove null Z values - they won't help us
wza_t =  wza_df[ ~wza_df.Z.isnull() ].reset_index(drop=True)
wza_s = wza_t.sort_values('SNPs').reset_index(drop=True)


# In[283]:


wza_s


# In[284]:


roller=50
minEntries=5
poly_deg = 12


# In[285]:


rolled_Z_vars = wza_s.Z.rolling(window=roller, min_periods=minEntries).var()


# In[286]:


masking_array = ~rolled_Z_vars.isnull()


# In[287]:


rolled_Z_sd = np.sqrt(rolled_Z_vars)[masking_array]


# In[288]:


rolled_Z_means = wza_s.Z.rolling(window=roller, min_periods=minEntries).mean()[masking_array]
rolled_mean_SNP_number = wza_s.SNPs.rolling(window=roller, min_periods=minEntries).mean()[masking_array]


# In[289]:


# Generating weights for polynomial function with degree=2 - standard deviation
sd_weights = np.polyfit(rolled_mean_SNP_number, rolled_Z_sd, deg=poly_deg)


# In[290]:


sd_polynomial_model = np.poly1d(sd_weights)


# In[291]:


# Generating weights for polynomial function with degree=2 - mean
mean_weights = np.polyfit(rolled_mean_SNP_number, rolled_Z_means, deg=poly_deg)
mean_polynomial_model = np.poly1d(mean_weights)

sd_predictions = sd_polynomial_model(wza_df["SNPs"])
#min_nonnegative_value = sd_predictions[sd_predictions >= 0].min()
#sd_predictions = np.clip(sd_predictions, a_min=min_nonnegative_value, a_max=None)
mean_predictions = mean_polynomial_model(wza_df["SNPs"])


# In[292]:


sd_predictions


# In[293]:


mean_predictions


# In[294]:


wza_p_values = [1 - norm.cdf(wza_df["Z"][i], loc=mean_predictions[i], scale=sd_predictions[i]) for i in range(wza_df.shape[0])]


# In[295]:


wza_df['Z_pVal'] = wza_p_values


# In[296]:


wza_df


# In[ ]:





# In[297]:


sns.scatterplot(data=wza_df, x='SNPs', y='Z_pVal', alpha = 0.2).set_yscale("log")
sns.rugplot(data=wza_df, x="SNPs", y="Z_pVal", height=.1, alpha=0.2).set_yscale("log")


# In[ ]:





# In[298]:


wza_p_values = [1 - norm.cdf(wza_df["Z"][i]) for i in range(wza_df.shape[0])]


# In[299]:


wza_df['Z_pVal_noweight'] = wza_p_values


# In[300]:


# Create a new column 'highlight' to identify rows where the gene is one of the specified ones
highlight_genes = [ '2_970']
wza_df['highlight'] = np.where(wza_df['gene'].isin(highlight_genes), 'red', 'grey')


sns.scatterplot(data=wza_df[wza_df['highlight'] == 'grey'], x='SNPs', y='Z_pVal_noweight', alpha = 0.2, color='grey',).set_yscale("log")
sns.scatterplot(data=wza_df[wza_df['highlight'] == 'red'], x='SNPs', y='Z_pVal_noweight', alpha = 0.2, color='red',).set_yscale("log")

sns.rugplot(data=wza_df, x="SNPs", y="Z_pVal_noweight", height=.1, alpha=0.2).set_yscale("log")


# In[ ]:





# In[203]:


# Create a new column 'highlight' to identify rows where the gene is one of the specified ones
highlight_genes = [ '2_970']
wza_df['highlight'] = np.where(wza_df['gene'].isin(highlight_genes), 'red', 'grey')


sns.scatterplot(data=wza_df[wza_df['highlight'] == 'grey'], x='SNPs', y='Z_pVal', alpha = 0.2, color='grey',).set_yscale("log")
sns.scatterplot(data=wza_df[wza_df['highlight'] == 'red'], x='SNPs', y='Z_pVal', alpha = 0.2, color='red',).set_yscale("log")

sns.rugplot(data=wza_df, x="SNPs", y="Z_pVal", height=.1, alpha=0.2).set_yscale("log")


# In[ ]:





# In[ ]:





# In[232]:


wza_df


# In[238]:


plt.figure(figsize=(10, 8))
sns.scatterplot(data = wza_df, x = 'Z_pVal', y = 'Z_pVal_noweight', hue = 'SNPs', alpha = 0.3)


# In[263]:


import seaborn as sns
import matplotlib.pyplot as plt

g = sns.JointGrid(data=wza_df, x = 'Z_pVal', y = 'Z_pVal_noweight')
g.plot(sns.scatterplot, sns.histplot)


# In[264]:


large_w = wza_df[wza_df['SNPs'] > 50]


# In[269]:


plt.figure(figsize=(10, 8))
sns.scatterplot(data = large_w, x = 'Z_pVal', y = 'Z_pVal_noweight', hue = 'SNPs', alpha = 0.5, palette='Spectral_r', legend=False)#.set_yscale("log").set_xscale("log")
sns.rugplot(data=large_w, x="Z_pVal", y="Z_pVal_noweight", height=.1, alpha=0.2,hue = 'SNPs',palette='Spectral_r', legend=False)#.set_yscale("log").set_xscale("log")


# In[273]:


large_w = wza_df[wza_df['SNPs'] > 500]


# In[301]:


plt.figure(figsize=(10, 8))
sns.scatterplot(data = large_w, x = 'Z_pVal', y = 'Z_pVal_noweight', hue = 'SNPs', alpha = 0.8, palette='Spectral_r', legend=False)#.set_yscale("log").set_xscale("log")
sns.rugplot(data=large_w, x="Z_pVal", y="Z_pVal_noweight", height=.1, alpha=0.8,hue = 'SNPs',palette='Spectral_r', legend=False)#.set_yscale("log").set_xscale("log")


# In[ ]:





# In[ ]:





# In[304]:


wza_df['scale'] = sd_predictions


# In[305]:


wza_df['mean'] = mean_predictions


# In[196]:


sns.scatterplot(data = wza_df, x = 'SNPs', y = 'mean')


# In[ ]:





# In[197]:


sns.scatterplot(data = wza_df, x = 'mean', y = 'scale')


# In[ ]:





# In[303]:


wza_df


# In[306]:


sns.scatterplot(data = wza_df, x = 'SNPs', y = 'scale', alpha=0.2).set_yscale("log")
sns.rugplot(data=wza_df, x="SNPs", y="scale", height=.1, alpha=0.2).set_yscale("log")


# In[ ]:





# In[307]:


# Create a new column 'highlight' to identify rows where the gene is one of the specified ones
highlight_genes = [ '2_970']
wza_df['highlight'] = np.where(wza_df['gene'].isin(highlight_genes), 'red', 'grey')

# Plot the non-highlighted (grey) points first
sns.scatterplot(
    data=wza_df[wza_df['highlight'] == 'grey'], 
    x='SNPs', 
    y='scale', 
    color='grey',  # Color for non-highlighted points
    alpha=0.5  # Transparency
).set_yscale("log")

# Now plot the highlighted (red) points on top
sns.scatterplot(
    data=wza_df[wza_df['highlight'] == 'red'], 
    x='SNPs', 
    y='scale', 
    color='red',  # Color for highlighted points
    alpha=1,  # Full opacity for better visibility
    zorder=3  # Ensure these points are drawn on top
).set_yscale("log")


# In[ ]:





# In[25]:


wza_df[wza_df['scale']>=0]['SNPs'].max()


# In[26]:


.set_yscale("log")


# In[29]:


sns.scatterplot(data = wza_df, x = 'top_candidate_p', y = 'Z_pVal')


# In[30]:


sns.histplot(wza_df['Z_pVal'])


# In[31]:


sns.histplot(wza_df['top_candidate_p'])


# In[ ]:





# In[ ]:





# In[27]:


sns.scatterplot(data = wza_df, x = 'scale', y = 'Z_pVal')


# In[160]:


wza_df


# In[ ]:





# In[73]:


len(wza_df[wza_df['scale']<0])


# In[ ]:





# In[631]:


#wza_df[wza_df['sd']<0].sort_values('sd')


# In[632]:


#pd.Series(sdev)[pd.Series(sdev) < 0]


# In[633]:


rolled_mean_SNP_number


# In[108]:


loc = mean_predictions[i]
loc


# In[122]:


i


# In[123]:


wza_df["SNPs"][i]


# In[161]:


scale = sd_predictions[i]
scale


# In[174]:


scale = 1


# In[175]:


z = wza_df["Z"][i]
z


# In[176]:


p = 1 - norm.cdf(z, loc = loc, scale =scale)


# In[177]:


p


# In[ ]:





# In[89]:


gotneg = []
for i in range(wza_df.shape[0]):
    loc = mean_predictions[i]
    scale = sd_predictions[i]
    p = 1 - norm.cdf(wza_df["Z"][i], loc = loc, scale =scale)
    if np.isnan(p):
        print(loc, scale, wza_df["Z"][i], rolled_mean_SNP_number[i])


# In[541]:


gotneg


# In[ ]:





# In[592]:


wza_results_lfmm_bio1_minentry20 = pd.read_csv('wza_results_lfmm_bio1_minentry0roll2_fix_log.csv')


# In[594]:


wza_results_lfmm_bio1_minentry20[wza_results_lfmm_bio1_minentry20['Z_pVal'].isna()]


# In[ ]:





# In[278]:


nan_blocks_0_min = wza_results_lfmm_bio1_minentry20[wza_results_lfmm_bio1_minentry20['Z_pVal'].isna()]['gene'].values


# In[ ]:





# In[81]:


wza_lfmm_test = pd.read_csv('wza_results_lfmm_bio1_test.csv')


# In[240]:


wza_lfmm = pd.read_csv('wza_results_lfmm_bio1.csv')


# In[241]:


wza_lfmm


# In[268]:


wza_lfmm


# In[244]:


wza_lfmm[wza_lfmm['Z_pVal'].isna()]['gene'].values


# In[89]:


problematic_windows = pd.read_csv('problematic_windows.csv')


# In[ ]:





# In[112]:


problematic_windows.sort_values('Z')


# In[101]:


problematic_windows['gene'].unique()


# In[103]:


snps_prob = lfmm[lfmm['block'].isin(problematic_windows['gene'].unique())]


# In[106]:


snps_prob.groupby('block').size().sort_values()


# In[ ]:





# In[99]:


lfmm.groupby('block').size().sort_values()[lfmm.groupby('block').size().sort_values() <= 40]


# In[ ]:





# In[83]:


wza_lfmm_test


# In[72]:


(wza_lfmm == wza_lfmm_test).all()


# In[73]:


wza_lfmm


# In[75]:


wza_lfmm_test[wza_lfmm_test['Z_pVal'].isna()]


# In[ ]:





# In[ ]:





# In[ ]:





# In[270]:


wza_lfmm = pd.read_csv('wza_results_lfmm_bio1_nocorrsnps.csv')


# In[397]:


wza_lfmm[wza_lfmm['Z_pVal'].isna()]


# In[ ]:





# In[398]:


## results
biovar = 'bio1'


# In[667]:


wza_lfmm = pd.read_csv('wza_results_lfmm_bio1_minentry0roll2_fix_log.csv').reset_index()


# In[668]:


wza_lfmm[wza_lfmm['Z_pVal'].isna()]


# In[637]:


'2_970'


# In[638]:


'1_702'


# In[ ]:





# In[ ]:





# In[639]:


threshold_value = 0.05 / len(wza_lfmm)


# In[640]:


#wza_lfmm[wza_lfmm['Z_pVal']< threshold_value].to_csv('wza_lfmm_significant_blocks.csv',index=None)


# In[641]:


observed_quantiles = -np.log10(np.sort(wza_lfmm['Z_pVal'].values))

# Expected quantiles from the uniform distribution
expected_quantiles = -np.log10(np.linspace(1 / len(wza_lfmm), 1, len(wza_lfmm)))

# QQ plot
sns.scatterplot(x = expected_quantiles, y = observed_quantiles, edgecolor='b', facecolor='none', alpha=0.5)
plt.plot([min(expected_quantiles), max(expected_quantiles)], [min(expected_quantiles), max(expected_quantiles)], 'r--')

plt.xlabel("Expected -log10(p-values)")
plt.ylabel("Observed -log10(p-values)")
plt.title(f'QQ Plot for {biovar}) WZA')

plt.show()


# In[642]:


nan_blocks_0_min


# In[ ]:





# In[669]:


old = pd.read_csv('wza_results_lfmm_bio1.csv')


# In[670]:


old = old[old['Z_pVal'] < threshold_value]


# In[671]:


new = wza_lfmm[wza_lfmm['Z_pVal'] < threshold_value]


# In[672]:


old


# In[673]:


old.merge(new, on ='gene', how='outer')


# In[674]:


new


# In[675]:


5_2922


# In[ ]:





# In[676]:


import pickle
dict_blocks = '../key_files/blocks_snpsid_dict.pkl'

with open(dict_blocks, 'rb') as file:
    dict_blocks = pickle.load(file)

reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}


# In[677]:


dict_blocks['5_2922'][-1] #5:24100585-24286474


# In[ ]:





# In[678]:


blocks_di


# In[644]:


set(wza_lfmm.loc[wza_lfmm['Z_pVal'] == 0]['gene']).intersection(set(nan_blocks_0_min))


# In[645]:


wza_lfmm.loc[wza_lfmm['Z_pVal'] == 0]


# In[679]:


wza_lfmm['chrom'] = wza_lfmm['gene'].str.split('_').str[0].astype(int)
wza_lfmm['pos'] = wza_lfmm['gene'].str.split('_').str[1].astype(int)


# In[680]:


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Significance threshold
threshold_value = 0.05 / len(wza_lfmm)
threshold = -np.log10(threshold_value)
biovar = 'bio1'

# Create chrom_pos in wza_lfmm by combining 'chrom' and 'pos'
wza_lfmm['chrom_pos'] = wza_lfmm['chrom'].astype(str) + '_' + wza_lfmm['pos'].astype(str)

# Copy the relevant columns for plotting
df = wza_lfmm[['Z_pVal', 'pos', 'chrom', 'chrom_pos']].copy()

# Parse chromosome number and position
df['chromosome'] = df['chrom']
df['position'] = df['pos']
df['-log10(pvalue)'] = -np.log10(df['Z_pVal'])

# Define colors for each chromosome
colors = sns.color_palette("crest", n_colors=5)

# Calculate chromosome offsets to prevent overlap
chromosome_offsets = {}
offset = 0
for chrom in sorted(df['chromosome'].unique()):
    chromosome_offsets[chrom] = offset
    max_position = df[df['chromosome'] == chrom]['position'].max()
    offset += max_position + 200  # Add buffer to prevent overlap

# Apply offsets to the position
df['adjusted_position'] = df.apply(lambda row: row['position'] + chromosome_offsets[row['chromosome']], axis=1)

# Create the Manhattan plot
plt.figure(figsize=(20, 6))

# Plot each chromosome separately
for chrom in sorted(df['chromosome'].unique()):
    subset = df[df['chromosome'] == chrom]
    plt.scatter(
        subset['adjusted_position'],
        subset['-log10(pvalue)'],
        alpha=0.7, 
        c=colors[chrom % len(colors)], 
        label=f'Chr {chrom}',
        s=20
    )

# Aesthetics
plt.xlabel('Adjusted Position')
plt.ylabel('-log10(pvalue)')
plt.title(f'{biovar} Manhattan Plot')
plt.grid(axis='y')

# Add a threshold line for significance
plt.axhline(y=threshold, color='grey', linestyle='dashed')

# Identify significant blocks in wza_lfmm
significant_blocks = df[df['-log10(pvalue)'] >= threshold]


# Highlight specific genes if needed
genes = ['2_1265', '4_801']  # List of specific genes to highlight
for gene in genes:
    chrom, pos = map(int, gene.split('_'))
    subset_gene = df[(df['chromosome'] == chrom) & (df['position'] == pos)]
    if not subset_gene.empty:
        plt.scatter(
            subset_gene['adjusted_position'],
            subset_gene['-log10(pvalue)'],
            edgecolor='red',
            linewidth=2,
            facecolor='none',
            s=100,
            label=f'Gene {gene}'
        )

plt.tight_layout()
plt.show()


# In[ ]:





# In[44]:


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Significance threshold
threshold_value = 0.05 / len(wza_lfmm)
threshold = -np.log10(threshold_value)
biovar = 'bio1'

# Create chrom_pos in wza_lfmm by combining 'chrom' and 'pos'
wza_lfmm['chrom_pos'] = wza_lfmm['chrom'].astype(str) + '_' + wza_lfmm['pos'].astype(str)

# Copy the relevant columns for plotting
df = wza_lfmm[['Z_pVal', 'pos', 'chrom', 'chrom_pos']].copy()

# Parse chromosome number and position
df['chromosome'] = df['chrom']
df['position'] = df['pos']
df['-log10(pvalue)'] = -np.log10(df['Z_pVal'])

# Define colors for each chromosome
colors = sns.color_palette("crest", n_colors=5)

# Calculate chromosome offsets to prevent overlap
chromosome_offsets = {}
offset = 0
for chrom in sorted(df['chromosome'].unique()):
    chromosome_offsets[chrom] = offset
    max_position = df[df['chromosome'] == chrom]['position'].max()
    offset += max_position + 200  # Add buffer to prevent overlap

# Apply offsets to the position
df['adjusted_position'] = df.apply(lambda row: row['position'] + chromosome_offsets[row['chromosome']], axis=1)

# Create the Manhattan plot
plt.figure(figsize=(20, 6))

# Plot each chromosome separately
for chrom in sorted(df['chromosome'].unique()):
    subset = df[df['chromosome'] == chrom]
    plt.scatter(
        subset['adjusted_position'],
        subset['-log10(pvalue)'],
        alpha=0.7, 
        c=colors[chrom % len(colors)], 
        label=f'Chr {chrom}',
        s=20
    )

# Aesthetics
plt.xlabel('Adjusted Position')
plt.ylabel('-log10(pvalue)')
plt.title(f'{biovar} Manhattan Plot')
plt.grid(axis='y')

# Add a threshold line for significance
plt.axhline(y=threshold, color='grey', linestyle='dashed')

# Identify significant blocks in wza_lfmm
significant_blocks = df[df['-log10(pvalue)'] >= threshold]

# Merge significant blocks with the annotation dataframe based on chrom_pos and block
annotated_blocks = significant_blocks.merge(small, left_on='chrom_pos', right_on='block', how='inner')

# Annotate the significant points with 'description1' and space annotations vertically
for chrom_pos, group in annotated_blocks.groupby('chrom_pos'):
    # Sort the group to ensure consistent spacing
    group = group.sort_values(by='description1')
    vertical_offset = 0  # Start the vertical offset at 0
    for i, row in group.iterrows():
        plt.annotate(
            row['description1'], 
            (row['adjusted_position'], row['-log10(pvalue)']),
            textcoords="offset points",  # Specify the offset point for the text
            xytext=(0, 5 + vertical_offset),  # Increment vertical offset for each annotation
            ha='center',    # Align horizontally to center
            fontsize=8,
            color='black'
        )
        vertical_offset += 15  # Increase vertical offset for the next annotation

# Highlight specific genes if needed
genes = ['2_1265', '4_801']  # List of specific genes to highlight
for gene in genes:
    chrom, pos = map(int, gene.split('_'))
    subset_gene = df[(df['chromosome'] == chrom) & (df['position'] == pos)]
    if not subset_gene.empty:
        plt.scatter(
            subset_gene['adjusted_position'],
            subset_gene['-log10(pvalue)'],
            edgecolor='red',
            linewidth=2,
            facecolor='none',
            s=100,
            label=f'Gene {gene}'
        )

plt.tight_layout()
plt.show()


# In[ ]:





# In[45]:


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Significance threshold
threshold_value = 0.05 / len(wza_lfmm)
threshold = -np.log10(threshold_value)
biovar = 'bio1'

# Create chrom_pos in wza_lfmm by combining 'chrom' and 'pos'
wza_lfmm['chrom_pos'] = wza_lfmm['chrom'].astype(str) + '_' + wza_lfmm['pos'].astype(str)

# Copy the relevant columns for plotting
df = wza_lfmm[['Z_pVal', 'pos', 'chrom', 'chrom_pos']].copy()

# Parse chromosome number and position
df['chromosome'] = df['chrom']
df['position'] = df['pos']
df['-log10(pvalue)'] = -np.log10(df['Z_pVal'])

# Define colors for each chromosome or block
colors = sns.color_palette("crest", n_colors=len(df['chromosome'].unique()))  # Assign a unique color per chromosome

# Calculate chromosome offsets to prevent overlap
chromosome_offsets = {}
offset = 0
for chrom in sorted(df['chromosome'].unique()):
    chromosome_offsets[chrom] = offset
    max_position = df[df['chromosome'] == chrom]['position'].max()
    offset += max_position + 200  # Add buffer to prevent overlap

# Apply offsets to the position
df['adjusted_position'] = df.apply(lambda row: row['position'] + chromosome_offsets[row['chromosome']], axis=1)

# Create the Manhattan plot
plt.figure(figsize=(20, 6))

# Plot each chromosome separately
for chrom in sorted(df['chromosome'].unique()):
    subset = df[df['chromosome'] == chrom]
    plt.scatter(
        subset['adjusted_position'],
        subset['-log10(pvalue)'],
        alpha=0.7, 
        c=colors[chrom % len(colors)], 
        label=f'Chr {chrom}',
        s=20
    )

# Aesthetics
plt.xlabel('Adjusted Position')
plt.ylabel('-log10(pvalue)')
plt.title(f'{biovar} Manhattan Plot')
plt.grid(axis='y')

# Add a threshold line for significance
plt.axhline(y=threshold, color='grey', linestyle='dashed')

# Identify significant blocks in wza_lfmm
significant_blocks = df[df['-log10(pvalue)'] >= threshold]

# Merge significant blocks with the annotation dataframe based on chrom_pos and block
annotated_blocks = significant_blocks.merge(small, left_on='chrom_pos', right_on='block', how='inner')

# Annotate the significant points with 'description1', space annotations, and color-code them
for chrom_pos, group in annotated_blocks.groupby('chrom_pos'):
    group = group.sort_values(by='description1')
    vertical_offset = 0  # Start the vertical offset at 0
    
    # Color for this block
    color = colors[chromosome_offsets[group['chromosome'].iloc[0]] % len(colors)]

    for i, row in group.iterrows():
        # Add the annotation
        plt.annotate(
            row['description1'], 
            (row['adjusted_position'], row['-log10(pvalue)']),
            textcoords="offset points",  
            xytext=(0, 5 + vertical_offset),  # Offset vertically for each annotation
            ha='center',
            fontsize=8,
            color=color  # Use the same color for annotation
        )
        
        # Draw a line connecting the point to the annotation
        plt.plot(
            [row['adjusted_position'], row['adjusted_position']],  # X coordinates (same X)
            [row['-log10(pvalue)'], row['-log10(pvalue)'] + (vertical_offset / 20)],  # Y coordinates
            color=color, linewidth=1, linestyle='dotted'
        )
        
        vertical_offset += 15  # Increase vertical offset for the next annotation

# Highlight specific genes if needed
genes = ['2_1265', '4_801']  # List of specific genes to highlight
for gene in genes:
    chrom, pos = map(int, gene.split('_'))
    subset_gene = df[(df['chromosome'] == chrom) & (df['position'] == pos)]
    if not subset_gene.empty:
        plt.scatter(
            subset_gene['adjusted_position'],
            subset_gene['-log10(pvalue)'],
            edgecolor='red',
            linewidth=2,
            facecolor='none',
            s=100,
            label=f'Gene {gene}'
        )

plt.tight_layout()
plt.show()


# In[ ]:





# In[27]:


threshold_value = 0.05 / len(wza_lfmm)
biovar='bio1'
#sm.qqplot(pvalues['pvalue'], line ='45') 
#py.show() 

df = wza_lfmm[['Z_pVal','pos','chrom']].copy()


# Parsing chromosome number and position
df['chromosome'] = df['chrom']
df['position'] = df['pos']
df['-log10(pvalue)'] = -np.log10(df['Z_pVal'])

colors = sns.color_palette("crest", n_colors = 5)

# Calculate the offset for each chromosome to prevent overlap
chromosome_offsets = {}
offset = 0
for chrom in sorted(df['chromosome'].unique()):
    chromosome_offsets[chrom] = offset
    max_position = df[df['chromosome'] == chrom]['position'].max()
    offset += max_position + 200  # Buffer to prevent overlap

# Apply offsets to positions
df['adjusted_position'] = df.apply(lambda row: row['position'] + chromosome_offsets[row['chromosome']], axis=1)

# Normalize sizes for better visualization
size_transform = 2  # Adjust this factor as needed


# Create a color map based on `n_est`
#df['color'] = df['n_est'].map(lambda x: cmap(norm(x)))

# Creating the Manhattan plot
plt.figure(figsize=(20, 6))

for chrom in sorted(df['chromosome'].unique()):
    subset = df[df['chromosome'] == chrom]
    plt.scatter(
        subset['adjusted_position'],
        subset['-log10(pvalue)'],
        alpha=0.7,  # Transparency for better visibility
        c=colors[chrom % len(colors)], 
        label=f'Chr {chrom}',
            s= 20)

# Aesthetics
plt.xlabel('Adjusted Position')
plt.ylabel('-log10(pvalue)')
plt.title(f'{biovar} Manhattan Plot')  # Set the title
plt.grid(axis='y')


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
# Create a legend for the number of estimated lineages
#handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=cmap(norm(n)), markersize=10, label=f'Lineages {n}') for n in sorted(n_est_unique)]
#plt.legend(handles=handles, title="Estimated Lineages", bbox_to_anchor=(1.05, 1), loc='upper left')

# Threshold line
threshold = -np.log10(threshold_value)
plt.axhline(y=threshold, color='grey', linestyle='dashed')

# Show the plot
plt.tight_layout()
plt.show()


# In[ ]:





# In[ ]:





# In[51]:


wza_kendall = pd.read_csv('wza_kendalltau_results_bio1.csv').reset_index()


# In[52]:


wza_kendall['chrom'] = wza_kendall['gene'].str.split('_').str[0].astype(int)
wza_kendall['pos'] = wza_kendall['gene'].str.split('_').str[1].astype(int)


# In[53]:


threshold_value = 0.05 / len(wza_kendall)


# In[54]:


wza_kendall[wza_kendall['Z_pVal']< threshold_value].to_csv('wza_kendall_significant_blocks.csv',index=None)


# In[ ]:


wza_kendall


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


gene - the name of the window
SNPs - the number of SNPs in this window
hits - the number of SNPs in the 99th percentile (not used for anything, just good to know)
Z - the Z score calculated for the gene
top_candidate_p - the result of the top-candidate method of Yeaman et al (2016 - Science)
LA - an indicator of whether the gene is causal for local adaptation
position - the average position of all SNPs in the window
Z_pVal - the p-value of the Z score (This is the WZA score)


# In[48]:


observed_quantiles = -np.log10(np.sort(wza_kendall['Z_pVal'].values))

# Expected quantiles from the uniform distribution
expected_quantiles = -np.log10(np.linspace(1 / len(wza_kendall), 1, len(wza_kendall)))

# QQ plot
sns.scatterplot(x = expected_quantiles, y = observed_quantiles, edgecolor='b', facecolor='none', alpha=0.5)
plt.plot([min(expected_quantiles), max(expected_quantiles)], [min(expected_quantiles), max(expected_quantiles)], 'r--')

plt.xlabel("Expected -log10(p-values)")
plt.ylabel("Observed -log10(p-values)")
plt.title(f'QQ Plot for {biovar}) WZA')

plt.show()


# In[ ]:





# In[ ]:





# In[57]:


threshold_value = 0.05 / len(wza_kendall)
biovar='bio1'
#sm.qqplot(pvalues['pvalue'], line ='45') 
#py.show() 

df = wza_kendall[['Z_pVal','pos','chrom']].copy()


# Parsing chromosome number and position
df['chromosome'] = df['chrom']
df['position'] = df['pos']
df['-log10(pvalue)'] = -np.log10(df['Z_pVal'])

colors = sns.color_palette("crest", n_colors = 5)

# Calculate the offset for each chromosome to prevent overlap
chromosome_offsets = {}
offset = 0
for chrom in sorted(df['chromosome'].unique()):
    chromosome_offsets[chrom] = offset
    max_position = df[df['chromosome'] == chrom]['position'].max()
    offset += max_position + 200  # Buffer to prevent overlap

# Apply offsets to positions
df['adjusted_position'] = df.apply(lambda row: row['position'] + chromosome_offsets[row['chromosome']], axis=1)

# Normalize sizes for better visualization
size_transform = 2  # Adjust this factor as needed


# Create a color map based on `n_est`
#df['color'] = df['n_est'].map(lambda x: cmap(norm(x)))

# Creating the Manhattan plot
plt.figure(figsize=(20, 6))

for chrom in sorted(df['chromosome'].unique()):
    subset = df[df['chromosome'] == chrom]
    plt.scatter(
        subset['adjusted_position'],
        subset['-log10(pvalue)'],
        alpha=0.7,  # Transparency for better visibility
        c=colors[chrom % len(colors)], 
        label=f'Chr {chrom}',
            s= 20)

# Aesthetics
plt.xlabel('Adjusted Position')
plt.ylabel('-log10(pvalue)')
plt.title(f'{biovar} Manhattan Plot')  # Set the title
plt.grid(axis='y')

# Create a legend for the number of estimated lineages
#handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=cmap(norm(n)), markersize=10, label=f'Lineages {n}') for n in sorted(n_est_unique)]
#plt.legend(handles=handles, title="Estimated Lineages", bbox_to_anchor=(1.05, 1), loc='upper left')
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


# Threshold line
threshold = -np.log10(threshold_value)
plt.axhline(y=threshold, color='grey', linestyle='dashed')

# Show the plot
plt.tight_layout()
plt.show()


# In[ ]:





# In[ ]:





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




