#!/usr/bin/env python
# coding: utf-8

"""
Required input files:
- '../wza_last_gen/wza_results_lfmm_{biovar}_poly7.csv': LFMM results
- '../wza_last_gen/wza_binomial_regression_{biovar}_poly7.csv': Binomial regression results
- '../wza_last_gen/wza_kendalltau_results_{biovar}_poly7.csv': Kendall tau results
- 'genes_info_BH_tair10_{biovar}.csv': Gene annotations

High-level outline:
This script generates Manhattan plots for GEA analysis results, specifically:
1. A plot showing all significant hits (gea_bh)
2. A plot showing one gene per block (gea_final_gen_bh_onegeneperblock)

The script processes results from three different methods (LFMM, Binomial Regression, and Kendall Tau)
and creates visualizations with proper chromosome spacing and significance thresholds.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
from adjustText import adjust_text

# Set the bioclimatic variable
biovar = 'bio18'  # This script is specifically for bio18

# Load and process GWAS data
def load_gwas_data(file_path):
    df = pd.read_csv(file_path)
    df['chrom'] = df['gene'].str.split('_').str[0].astype(int)
    df['pos'] = df['gene'].str.split('_').str[1].astype(int)
    df['chrom_pos'] = df['chrom'].astype(str) + '_' + df['pos'].astype(str)
    df['Bonferroni_corrected_pval'] = multipletests(df['Z_pVal'], method='fdr_bh')[1]
    df['-log10(pvalue)'] = -np.log10(df['Z_pVal'])
    return df[['Z_pVal', 'pos', 'chrom', 'chrom_pos', 'Bonferroni_corrected_pval', '-log10(pvalue)']]

# Load all three datasets
df1 = load_gwas_data(f'../wza_last_gen/wza_results_lfmm_{biovar}_poly7.csv')
df2 = load_gwas_data(f'../wza_last_gen/wza_binomial_regression_{biovar}_poly7.csv')
df3 = load_gwas_data(f'../wza_last_gen/wza_kendalltau_results_{biovar}_poly7.csv')

# Calculate chromosome offsets
chromosome_offsets = {}
offset = 0
for chrom in sorted(pd.concat([df1['chrom'], df2['chrom'], df3['chrom']]).unique()):
    chromosome_offsets[chrom] = offset
    max_position = max(df1[df1['chrom'] == chrom]['pos'].max(),
                      df2[df2['chrom'] == chrom]['pos'].max(),
                      df3[df3['chrom'] == chrom]['pos'].max())
    offset += max_position + 300

# Apply offsets
for df in [df1, df2, df3]:
    df['adjusted_position'] = df.apply(lambda row: row['pos'] + chromosome_offsets[row['chrom']], axis=1)

# Load annotations
annot = pd.read_csv(f'genes_info_BH_tair10_{biovar}.csv')
annot = annot[annot['gen'] != 'first_gen']
annot = annot.merge(df3[['chrom_pos', 'adjusted_position']], left_on='block_id', right_on='chrom_pos')

# Get significant points
significant_df1 = df1[df1['Bonferroni_corrected_pval'] <= 0.05]
significant_df2 = df2[df2['Bonferroni_corrected_pval'] <= 0.05]
significant_df3 = df3[df3['Bonferroni_corrected_pval'] <= 0.05]

# Create the first plot (gea_bh)
plt.figure(figsize=(20, 6))
colors = ['#2aad2a', '#208420', '#006400']  # Green shades for different methods

# Plot all points
for df, color in zip([df1, df2, df3], colors):
    plt.scatter(df['adjusted_position'], df['-log10(pvalue)'], 
               c=color, alpha=0.5, s=10, label=f'{df.name if hasattr(df, "name") else "Method"}')

# Add significance threshold
threshold = -np.log10(0.05)
plt.axhline(y=threshold, color='grey', linestyle='--', alpha=0.5)

# Customize plot
plt.xlabel('Chromosome Position')
plt.ylabel('-log10(p-value)')
plt.title(f'Manhattan Plot - {biovar}')
plt.legend()
plt.tight_layout()
plt.savefig(f'gea_bh_{biovar}.png', dpi=300, bbox_inches='tight')
plt.close()

# Create the second plot (gea_final_gen_bh_onegeneperblock)
plt.figure(figsize=(20, 6))

# Plot significant points
for df, color in zip([significant_df1, significant_df2, significant_df3], colors):
    plt.scatter(df['adjusted_position'], df['-log10(pvalue)'], 
               c=color, alpha=0.5, s=10, label=f'{df.name if hasattr(df, "name") else "Method"}')

# Add annotations
texts = []
for _, row in annot.iterrows():
    texts.append(plt.text(row['adjusted_position'], row['-log10(pvalue)'], 
                         row['gene_name'], fontsize=8))

# Adjust text positions to avoid overlap
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', lw=0.5))

# Add significance threshold
plt.axhline(y=threshold, color='grey', linestyle='--', alpha=0.5)

# Customize plot
plt.xlabel('Chromosome Position')
plt.ylabel('-log10(p-value)')
plt.title(f'Manhattan Plot with Annotations - {biovar}')
plt.legend()
plt.tight_layout()
plt.savefig(f'gea_final_gen_bh_onegeneperblock_{biovar}.png', dpi=300, bbox_inches='tight')
plt.close()


# In[ ]:




