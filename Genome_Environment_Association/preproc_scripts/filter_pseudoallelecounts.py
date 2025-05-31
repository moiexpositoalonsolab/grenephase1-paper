import pandas as pd
import dask.dataframe as dd

# Read the final generation samples
last_gen = pd.read_csv('key_files/final_gen.csv')['sample_name']

# Read allele frequencies for last generation samples
allele_freq = dd.read_csv('key_files/merged_hapFIRE_allele_frequency.txt', 
                         usecols=last_gen, 
                         sep='\t')

allele_freq = allele_freq.compute()
allele_freq = allele_freq.reset_index(drop=True)

# Read sample table for flower counts
fc = pd.read_csv('key_files/merged_sample_table.csv')

# Read SNP dictionary and apply MAF filter
dict_snps = pd.read_csv('key_files/var_pos_grenenet.csv')
mask_maf_05 = dict_snps['maf05filter'].notna()
dict_snps_filter = dict_snps[dict_snps['maf05filter'].notna()].drop(['maf05filter', 'total_alleles05filter'], axis=1).reset_index(drop=True)

# Filter allele frequencies by MAF
allele_freq_f = allele_freq[mask_maf_05]
allele_freq_f = allele_freq_f.reset_index(drop=True)

# Calculate allele counts
allele_counts = {}
for sample in allele_freq_f.columns:
    flower_counts = fc[fc['sample_name'] == sample]['total_flower_counts'].values[0]
    allele_counts[sample] = allele_freq_f[sample] * flower_counts * 2

allele_counts = pd.concat(allele_counts, axis=1)

# Calculate minimum count threshold
total_flowers_collected_last_gen = fc[fc['sample_name'].isin(last_gen)]['total_flower_counts'].sum()
total_genomes_collected_last_gen = total_flowers_collected_last_gen * 2
min_freq = 0.05
min_count = total_genomes_collected_last_gen * min_freq

# Filter by minimum count
total_alleles_sampled = allele_counts.sum(axis=1)
mask_min_counts = total_alleles_sampled > min_count
allele_freq_f_maf_min_count = allele_freq_f[mask_min_counts]

# Save filtered allele frequencies
allele_freq_f_maf_min_count.to_csv('key_files/allele_freq_maf05_mincount05_lastgensamples.csv', index=None) 