#!/usr/bin/env python
# coding: utf-8

# This script analyzes the effect of a specific SNP, likely associated with or within
# the TSF gene region, on flowering time phenotypes
# (measured at 10°C and 16°C growth chamber conditions).

# The script performs the following main steps:
# 1. Loads genotype data for a subset of SNPs (including the target SNP) and phenotype data.
# 2. Extracts genotype information for the specific SNP of interest.
# 3. Merges genotype data with phenotype data based on ecotype/sample identifiers.
# 4. Generates boxplots to visualize the relationship between the target SNP genotype
#    and flowering time phenotypes.

# --- Required Input Files ---
# - key_files/founder_ecotype_names.csv: List of founder ecotype names.
# - key_files/subset_TSF.recode.vcf: VCF file containing genotype data for a subset of SNPs (should include the target SNP).
# - key_files/atlas_phenotype_matrix_withid.csv: CSV file containing flowering time and other phenotype data with ecotype/sample IDs.

# --- Output Files ---
# - The script generates plots directly using `matplotlib.pyplot.show()`, which displays the plots.
#   If saving the plots to files is desired, `plt.savefig()` calls would need to be added.
#   Based on the analysis, it will produce two boxplots: one for FT10 and one for FT16.

# --- Import Libraries ---

import pandas as pd
import allel # For reading and working with VCF files
import seaborn as sns # For enhanced visualizations
import matplotlib.pyplot as plt # For plotting

# In[1]: # Import pandas (already done above, kept for notebook structure)


# In[2]: # Import allel (already done above, kept for notebook structure)


# In[30]: # Load founder ecotype names


founder_ecotype_names = pd.read_csv('key_files/founder_ecotype_names.csv')['0']


# In[3]: # Load VCF and phenotype data


# Load the subset VCF file containing genotype data.
tsf = allel.read_vcf('key_files/subset_TSF.recode.vcf')

# Load the phenotype matrix.
pheno = pd.read_csv('key_files/atlas_phenotype_matrix_withid.csv')

# Print the number of non-missing values for FT10 and FT16 (likely for data inspection).
# print(pheno['FT10'].notna().sum())
# print(pheno['FT16'].notna().sum())


# In[51]: # Extract genotype data for the target SNP


# Extract genotype data for all samples from the VCF.
genotype_data = tsf['calldata/GT']  # Assuming GT field contains genotype information

# Define the target SNP ID.
target_snp_id = '4_10999396'

# Find the index of the target SNP in the VCF file.
snp_index = list(tsf['variants/ID']).index(target_snp_id)

# Get genotypes for the specific SNP across all samples.
# The shape will be (number_of_samples, number_of_alleles), typically (samples, 2) for diploid data.
genotypes = genotype_data[:, snp_index, :]


# In[52]: # Calculate allele counts


# Calculate allele counts for each sample at all SNPs in the VCF.
# Summing across the last axis (alleles) gives the number of non-reference alleles for each sample at each SNP.
allele_counts = genotype_data.sum(axis=2)


# In[53]: # Transpose allele counts


# Transpose the allele counts matrix to have SNPs as columns and samples as rows.
allele_counts = allele_counts.T


# In[54]: # Convert allele counts to DataFrame


# Convert the transposed allele counts to a pandas DataFrame.
# Assign SNP IDs as column names.
df = pd.DataFrame(allele_counts, columns=tsf['variants/ID'])


# In[55]: # Add ecotype/sample IDs to the DataFrame


# Add a column for sample IDs (ecotype names) to the DataFrame.
df['id'] = founder_ecotype_names


# In[56]: # Select genotype data for the target SNP and ID


# Select only the column for the target SNP and the 'id' column.
# This DataFrame 'x4' now contains the genotype (allele count) at the target SNP
# and the corresponding sample ID.
x4 = df[[target_snp_id, 'id']]


# In[58]: # Merge phenotype and genotype data


# Merge the phenotype DataFrame ('pheno') with the target SNP genotype DataFrame ('x4').
# The merge is performed based on the 'id' column, linking phenotype data to genotypes.
merged_df = pd.merge(pheno, x4, on='id')


# In[60]: # Display merged DataFrame columns (likely for inspection)


# print(merged_df.columns) # Print column names of the merged DataFrame.


# In[61]: # Create and display boxplot for FT10


# Set the seaborn plot style.
sns.set(style="whitegrid")

# Create a figure for the plot.
plt.figure(figsize=(6, 6))

# Create a boxplot showing the distribution of FT10 for each genotype category at the target SNP.
# showcaps=True: Show the caps of the whiskers.
# fliersize=0: Do not show outlier points (fliers).
# width=0.5: Set the width of the boxes.
# boxprops=dict(alpha=0.6): Set transparency for the boxes.
sns.boxplot(x=target_snp_id, y='FT10', data=merged_df, showcaps=True, fliersize=0, width=0.5, boxprops=dict(alpha=0.6))

# Add jittered data points (strip plot) overlaid on the boxplot.
# jitter=True: Add random noise to prevent points from overlapping.
# alpha=0.4: Set transparency for the points.
# color='black': Set point color to black.
# edgecolor=None: Remove edges from points.
sns.stripplot(x=target_snp_id, y='FT10', data=merged_df, jitter=True, alpha=0.4, color='black', edgecolor=None)

# Label axes.
plt.xlabel(f"SNP {target_snp_id} Genotype (Allele Count)") # Dynamic label using the target SNP ID.
plt.ylabel("Flowering time in growth chamber (10°C)")

# Adjust layout to prevent labels/titles from overlapping.
plt.tight_layout()

# Display the plot.
plt.show()

# In[62]: # Create and display boxplot for FT16


# Set the seaborn plot style (repeated, but fine).
sns.set(style="whitegrid")

# Create a figure for the plot.
plt.figure(figsize=(6, 6))

# Create a boxplot for FT16, similar to the FT10 plot.
sns.boxplot(x=target_snp_id, y='FT16', data=merged_df, showcaps=True, fliersize=0, width=0.5, boxprops=dict(alpha=0.6))
sns.stripplot(x=target_snp_id, y='FT16', data=merged_df, jitter=True, alpha=0.4, color='black', edgecolor=None)

# Label axes.
plt.xlabel(f"SNP {target_snp_id} Genotype (Allele Count)") # Dynamic label.
plt.ylabel("Flowering time in growth chamber (16°C)") # Corrected temperature in label.

# Adjust layout and display the plot.
plt.tight_layout()
plt.show()


# In[65]: # Inspect data with non-missing FT16 (likely for debugging/verification)


# Display a subset of the merged DataFrame where FT16 values are not missing.
# This is likely a step for data inspection during development.
# print(merged_df[['4_10999396', 'FT16']].dropna())


# In[ ] sections are likely empty cells or remnants of notebook conversion and are removed.


