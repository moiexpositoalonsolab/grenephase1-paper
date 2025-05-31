#!/usr/bin/env python
# coding: utf-8

"""
Analysis and Annotation of CAM5 Variants and Ecotype Frequencies

This script processes VCF files for the CAM5 genomic region, annotates variants with functional
information and climate data, and visualizes relationships between variants, annotations,
climate, and ecotype frequencies using heatmaps.

Key steps include:
- Filtering 1001 Genomes and Grenenet VCFs for the CAM5 region.
- Annotating variants with effect, impact, and gene feature type.
- Mapping variants to ecotype climate information.
- Analyzing ecotype frequencies in relation to both ecotype origin climate and experimental
  site climate.
- Generating heatmaps to visualize genotype patterns and ecotype frequency distribution.

Inputs:
- cam5_annot.csv: Annotation data for CAM5 variants.
- ../key_files/1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv: Climate and ecotype info for 1001 Genomes accessions.
- ../gwas/1001g/1001g_grenet_climate.recode.vcf: 1001 Genomes VCF file.
- ../gwas/greneNet_final_v1.1.recode.vcf.gz: Grenenet VCF file (gzipped).
- ../key_files/merged_sample_table.csv: Sample metadata for Grenenet.
- ../key_files/merged_ecotype_frequency.txt: Merged ecotype frequency data for Grenenet.
- ../key_files/bioclimvars_experimental_sites_era5.csv: ERA5 climate data for experimental sites.
- Data containing genomic feature coordinates for CAM5 (loaded as exons_utr_df).

Outputs:
- filtered_cam5_1001g.vcf: Filtered 1001 Genomes VCF for CAM5 region.
- filtered_cam5_grenenet.vcf: Filtered Grenenet VCF for CAM5 region.
- Interactive heatmaps visualizing variant data and ecotype frequency patterns.
"""

# --- Import Libraries ---
import pandas as pd
import allel
import seaborn as sns
import numpy as np
import subprocess
import matplotlib.pyplot as plt

# --- Configuration ---
# Define genomic region for CAM5 (Chromosome 2, approximate coordinates)
CAM5_CHROM = '2'
CAM5_START_POS = 11531967
CAM5_END_POS = 11534358

# --- Load Data ---
# Load CAM5 annotation data
cam5_annot = pd.read_csv('cam5_annot.csv')

# Load 1001 Genomes ecotype climate info
clim1001 = pd.read_csv('../key_files/1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv')
clim1001['ecotypeid'] = clim1001['ecotypeid'].astype(str)
ecotypeid_to_bio1 = clim1001.set_index('ecotypeid')['bio1'].to_dict()

# Load sample metadata for Grenenet
samples_meta = pd.read_csv('../key_files/merged_sample_table.csv')[['site', 'plot', 'generation', 'total_flower_counts']]
samples_meta = samples_meta.groupby(['site', 'plot', 'generation'])['total_flower_counts'].sum().reset_index()
samples_meta['code'] = samples_meta['site'].astype(str) + '_'  + samples_meta['generation'].astype(str) + '_' + samples_meta['plot'].astype(str)

# Load merged ecotype frequency data
ecotype_freq_raw = pd.read_csv('../key_files/merged_ecotype_frequency.txt', sep = '\t')

# Load experimental site climate data
clim_sites_during_exp = pd.read_csv('/carnegie/nobackup/scratch/tbellagio/grene/data/worldclim_sitesdata.csv')[['site', 'bio1']]
clim_sites_during_exp['site'] = clim_sites_during_exp['site'].astype(str)

# Assume exons_utr_df is loaded elsewhere or define a placeholder if needed
# Placeholder/Example structure based on usage in original script:
exons_utr_df = pd.DataFrame({
    'seqid': [CAM5_CHROM, CAM5_CHROM, CAM5_CHROM, CAM5_CHROM, CAM5_CHROM, CAM5_CHROM],
    'feature_type': ['exon', 'UTR', 'exon', 'UTR', 'exon', 'UTR'],
    'start': [11531967, 11531967, 11532500, 11532500, 11533000, 11533000],
    'end': [11532100, 11532100, 11532650, 11532650, 11533150, 11533150],
    'annotation': ['Parent=AT2G27030.1', 'Parent=AT2G27030.2', 'Parent=AT2G27030.3', 'Parent=AT2G27030.1', 'Parent=AT2G27030.2', 'Parent=AT2G27030.3']
})
# Convert start and end into an IntervalIndex
exons_utr_df['interval'] = pd.IntervalIndex.from_arrays(exons_utr_df['start'], exons_utr_df['end'], closed='both')
exons_utr_df.set_index('interval', inplace=True)

# --- Filter VCF files for CAM5 Region ---
print(f"Filtering VCF files for CAM5 region ({CAM5_CHROM}:{CAM5_START_POS}-{CAM5_END_POS})...")

# Filter 1001 Genomes VCF
command_1001g = f"grep \"^#\" ../gwas/1001g/1001g_grenet_climate.recode.vcf > filtered_cam5_1001g.vcf && grep \"^{CAM5_CHROM}\\s\" ../gwas/1001g/1001g_grenet_climate.recode.vcf | awk \'$2 >= {CAM5_START_POS} && $2 <= {CAM5_END_POS}\' >> filtered_cam5_1001g.vcf"
subprocess.run(command_1001g, shell=True, executable='/bin/bash', check=True)
print("Filtered 1001 Genomes VCF saved to filtered_cam5_1001g.vcf")

# Filter Grenenet VCF
command_grene = f"zgrep \"^#\" ../gwas/greneNet_final_v1.1.recode.vcf.gz> filtered_cam5_grenenet.vcf && zgrep \"^{CAM5_CHROM}\\s\" ../gwas/greneNet_final_v1.1.recode.vcf.gz | awk \'$2 >= {CAM5_START_POS} && $2 <= {CAM5_END_POS}\' >> filtered_cam5_grenenet.vcf"
subprocess.run(command_grene, shell=True, executable='/bin/bash', check=True)
print("Filtered Grenenet VCF saved to filtered_cam5_grenenet.vcf")

# --- Process 1001 Genomes VCF and Annotate ---
print("Processing 1001 Genomes VCF...")

# Read filtered 1001 Genomes VCF using allel
vcf1001 = allel.read_vcf('filtered_cam5_1001g.vcf')

# Extract genotype data and convert to DataFrame
geno1001 = allel.create_genotype_array(vcf1001['calldata/GT']).to_allele_counts()[:].sum(axis=2)
geno1001_df = pd.DataFrame(geno1001, columns=vcf1001['samples'])

# Combine variant position with genotype data
vcf1001_df = pd.DataFrame({'chrom': vcf1001['variants/CHROM'], 'pos': vcf1001['variants/POS']})
vcf1001_df = pd.concat([vcf1001_df, geno1001_df], axis=1)

# Merge with CAM5 annotation
vcf_1001_cam5annot = cam5_annot.merge(vcf1001_df, on = 'pos')
vcf_1001_cam5annot = vcf_1001_cam5annot.drop('chrom_x',axis=1)

# Map ecotype IDs to bio1 climate values for columns
bios_maped_cols = vcf_1001_cam5annot.columns[5:].map(ecotypeid_to_bio1)
new_columns = ['pos', 'EFF', 'LOF', 'NMD', 'chrom_y'] + bios_maped_cols.tolist()
vcf_1001_cam5annot.columns = new_columns

# Parse 'EFF' column for effect and impact
eff_columns = vcf_1001_cam5annot['EFF'].str.split('|', expand=True)
eff_columns.columns = ['Effect', 'Impact', 'Position', 'Change', 'Gene', 'Biotype', 'Rank', 'HGVS.c', 'HGVS.p', 'cdna_pos', 'Other']
eff_columns = eff_columns.drop(['Position', 'Change', 'Gene', 'Biotype', 'cdna_pos', 'Other'], axis=1)
vcf_1001_cam5annot = pd.concat([eff_columns, vcf_1001_cam5annot], axis=1)

# Annotate variants within the significant block region
SIGN_BLOCK_START = 11533904
SIGN_BLOCK_END = 11534263
vcf_1001_cam5annot.loc[(vcf_1001_cam5annot['pos']>= SIGN_BLOCK_START) & (vcf_1001_cam5annot['pos']<= SIGN_BLOCK_END), 'block'] = 'sign_block'

# Map SNP position to gene feature type based on different isoforms
# Filter exons_utr_df for each isoform and map positions
exons_utr_df_var1 = exons_utr_df[exons_utr_df['annotation'] == 'Parent=AT2G27030.1']
vcf_1001_cam5annot['feature_type_var1'] = vcf_1001_cam5annot['pos'].apply(
    lambda x: exons_utr_df_var1.loc[exons_utr_df_var1.index.contains(x)]['feature_type'].values[0] if any(exons_utr_df_var1.index.contains(x)) else 'intergenic'
)

exons_utr_df_var2 = exons_utr_df[exons_utr_df['annotation'] == 'Parent=AT2G27030.2']
vcf_1001_cam5annot['feature_type_var2'] = vcf_1001_cam5annot['pos'].apply(
    lambda x: exons_utr_df_var2.loc[exons_utr_df_var2.index.contains(x)]['feature_type'].values[0] if any(exons_utr_df_var2.index.contains(x)) else 'intergenic'
)

exons_utr_df_var3 = exons_utr_df[exons_utr_df['annotation'] == 'Parent=AT2G27030.3']
vcf_1001_cam5annot['feature_type_var3'] = vcf_1001_cam5annot['pos'].apply(
    lambda x: exons_utr_df_var3.loc[exons_utr_df_var3.index.contains(x)]['feature_type'].values[0] if any(exons_utr_df_var3.index.contains(x)) else 'intergenic'
)

# Filter to keep only variants in the significant block
vcf_1001_cam5annot_filtered = vcf_1001_cam5annot[vcf_1001_cam5annot['block'] == 'sign_block'].copy()

# --- Visualize 1001 Genomes Data (Heatmaps) ---
print("Generating heatmaps for 1001 Genomes data...")

# Define columns to use as index for heatmaps
heatmap_indices_1001g = ['block', 'feature_type_var1', 'feature_type_var3', 'Impact']

for index_col in heatmap_indices_1001g:
    if index_col in vcf_1001_cam5annot_filtered.columns:
        # Prepare data for heatmap
        for_hm = vcf_1001_cam5annot_filtered.set_index(index_col)
        
        # Select genotype columns (assuming they start from index 5, exclude last 4 annotation columns)
        genotype_cols_start = 5 # Adjust if columns before genotype change
        annotation_cols_end = 4 # Number of annotation columns at the end (block, feature_types, Impact)
        genotype_data = for_hm.iloc[:, genotype_cols_start:-annotation_cols_end].T
        
        # Round index (bio1 values) for clarity
        genotype_data.index = np.round(genotype_data.index.astype(float), 2)
        for_hm_final = genotype_data.T

        # Create and display heatmap
        plt.figure(figsize=(16, 8))
        sns.heatmap(for_hm_final, cmap="viridis", cbar=False)
        plt.title(f'Genotype Heatmap (1001 Genomes) indexed by {index_col}')
        plt.xlabel('Ecotype Origin Bio1')
        plt.ylabel(index_col)
        plt.tight_layout()
        plt.show()
    else:
        print(f"Warning: Index column \'{index_col}\' not found in filtered 1001 Genomes data.")

# --- Process Grenenet VCF and Annotate ---
print("Processing Grenenet VCF...")

# Read filtered Grenenet VCF using allel
vcfgrene = allel.read_vcf('filtered_cam5_grenenet.vcf')

# Extract genotype data and convert to DataFrame
geno_grene = allel.create_genotype_array(vcfgrene['calldata/GT']).to_allele_counts()[:].sum(axis=2)
geno_grene_df = pd.DataFrame(geno_grene, columns=vcfgrene['samples'])

# Combine variant position with genotype data
vcfgrene_df = pd.DataFrame({'chrom': vcfgrene['variants/CHROM'], 'pos': vcfgrene['variants/POS']})
vcfgrene_df = pd.concat([vcfgrene_df, geno_grene_df], axis=1)

# Merge with CAM5 annotation
vcf_grene_cam5annot = cam5_annot.merge(vcfgrene_df, on = 'pos')
vcf_grene_cam5annot = vcf_grene_cam5annot.drop('chrom_x',axis=1)

# Map sample IDs to bio1 climate values for columns (using experimental site climate)
# This mapping is different from the 1001 Genomes data
sample_to_bio1_exp = samples_meta.set_index('code')['site'].map(clim_sites_during_exp.set_index('site')['bio1']).to_dict()

bios_maped_cols_grene = vcf_grene_cam5annot.columns[5:].map(sample_to_bio1_exp)
new_columns_grene = ['pos', 'EFF', 'LOF', 'NMD', 'chrom_y'] + bios_maped_cols_grene.tolist()
vcf_grene_cam5annot.columns = new_columns_grene

# Parse 'EFF' column for effect and impact
eff_columns_grene = vcf_grene_cam5annot['EFF'].str.split('|', expand=True)
eff_columns_grene.columns = ['Effect', 'Impact', 'Position', 'Change', 'Gene', 'Biotype', 'Rank', 'HGVS.c', 'HGVS.p', 'cdna_pos', 'Other']
eff_columns_grene = eff_columns_grene.drop(['Position', 'Change', 'Gene', 'Biotype', 'cdna_pos', 'Other'], axis=1)
vcf_grene_cam5annot = pd.concat([eff_columns_grene, vcf_grene_cam5annot], axis=1)

# Annotate variants within the significant block region (using same coordinates as 1001G)
vcf_grene_cam5annot.loc[(vcf_grene_cam5annot['pos']>= SIGN_BLOCK_START) & (vcf_grene_cam5annot['pos']<= SIGN_BLOCK_END), 'block'] = 'sign_block'

# Map SNP position to gene feature type based on different isoforms (using same exons_utr_df)
vcf_grene_cam5annot['feature_type_var1'] = vcf_grene_cam5annot['pos'].apply(
    lambda x: exons_utr_df_var1.loc[exons_utr_df_var1.index.contains(x)]['feature_type'].values[0] if any(exons_utr_df_var1.index.contains(x)) else 'intergenic'
)
vcf_grene_cam5annot['feature_type_var2'] = vcf_grene_cam5annot['pos'].apply(
    lambda x: exons_utr_df_var2.loc[exons_utr_df_var2.index.contains(x)]['feature_type'].values[0] if any(exons_utr_df_var2.index.contains(x)) else 'intergenic'
)
vcf_grene_cam5annot['feature_type_var3'] = vcf_grene_cam5annot['pos'].apply(
    lambda x: exons_utr_df_var3.loc[exons_utr_df_var3.index.contains(x)]['feature_type'].values[0] if any(exons_utr_df_var3.index.contains(x)) else 'intergenic'
)

# Filter to keep only variants in the significant block
vcf_grene_cam5annot_filtered = vcf_grene_cam5annot[vcf_grene_cam5annot['block'] == 'sign_block'].copy()

# --- Visualize Grenenet Data (Heatmaps) ---
print("Generating heatmaps for Grenenet data...")

# Define columns to use as index for heatmaps (using a subset for Grenenet)
heatmap_indices_grene = ['Impact', 'feature_type_var1', 'feature_type_var3'] # Exclude 'block' as it's constant after filtering

for index_col in heatmap_indices_grene:
     if index_col in vcf_grene_cam5annot_filtered.columns:
        # Prepare data for heatmap
        for_hm_grene = vcf_grene_cam5annot_filtered.set_index(index_col)
        
        # Select genotype columns (assuming they start from index 5, exclude last annotation columns)
        genotype_cols_start = 5 # Adjust if columns before genotype change
        annotation_cols_end = 4 # Number of annotation columns at the end (block, feature_types, Impact)
        genotype_data_grene = for_hm_grene.iloc[:, genotype_cols_start:-annotation_cols_end].T
        
        # Round index (bio1 values) for clarity
        genotype_data_grene.index = np.round(genotype_data_grene.index.astype(float), 2)
        for_hm_grene_final = genotype_data_grene.T

        # Create and display heatmap
        plt.figure(figsize=(16, 8))
        sns.heatmap(for_hm_grene_final, cmap="viridis", cbar_label='Genotype Count') # Add cbar label
        plt.title(f'Genotype Heatmap (Grenenet) indexed by {index_col}')
        plt.xlabel('Experimental Site Bio1')
        plt.ylabel(index_col)
        plt.tight_layout()
        plt.show()
     else:
        print(f"Warning: Index column \'{index_col}\' not found in filtered Grenenet data.")

# --- Analyze Ecotype Frequencies vs. Climate ---
print("Analyzing ecotype frequencies vs. climate...")

# Process raw ecotype frequency data
# The original script used grenenet_ecotypes as index, but merged_ecotype_frequency.txt
# seems to have samples as columns. Assuming samples as columns for melting.

ecotype_freq = ecotype_freq_raw.T.reset_index()
ecotype_freq.columns = ecotype_freq.columns.astype(str) # Ensure columns are strings for melting

# Extract generation, plot, and site from sample names
ecotype_freq['generation'] = ecotype_freq['index'].str.split('_').str[1]
ecotype_freq['plot'] = ecotype_freq['index'].str.split('_').str[2]
ecotype_freq['site'] = ecotype_freq['index'].str.split('_').str[0]

# Melt the DataFrame to long format
ecotype_freq_long = ecotype_freq.melt(id_vars = ['index', 'plot','site', 'generation'])
ecotype_freq_long.columns = ['sample', 'plot','site', 'generation', 'ecotypeid', 'freq']

# Merge with ecotype origin climate data
ecotype_freq_climate = ecotype_freq_long.merge(clim1001[['ecotypeid', 'bio1']], on = 'ecotypeid', how='left')
ecotype_freq_climate.rename(columns={'bio1': 'bio1_ecotype_origin'}, inplace=True)

# Merge with experimental site climate data
ecotype_freq_climate = ecotype_freq_climate.merge(clim_sites_during_exp, on ='site', how='left')
ecotype_freq_climate.rename(columns={'bio1': 'bio1_experimental_site'}, inplace=True)

# Round bio1 values for plotting
ecotype_freq_climate['bio1_experimental_site'] = ecotype_freq_climate['bio1_experimental_site'].round(2)
ecotype_freq_climate['bio1_ecotype_origin'] = ecotype_freq_climate['bio1_ecotype_origin'].round(2)

# Create pivot table for heatmap (Median Frequency)
bioexp_vs_bioecotype = ecotype_freq_climate.pivot_table(index = 'bio1_experimental_site', columns = 'bio1_ecotype_origin', values= 'freq' , aggfunc = 'median')

# Sort index and columns
bioexp_vs_bioecotype = bioexp_vs_bioecotype.sort_index() 
sorted_columns = sorted(bioexp_vs_bioecotype.columns)
bioexp_vs_bioecotype = bioexp_vs_bioecotype[sorted_columns]

# Map bio1 values to combined strings for clearer labels (optional, depends on desired visual output)
# Assuming clim1001 has 'ecotypeid' and 'bio1' and clim_sites_during_exp has 'site' and 'bio1'
# This part of the original script seemed intended to relabel the heatmap axes with more info,
# but the mapping logic was complex. Keeping the rounded bio1 values for simplicity for now.
# If specific string labels are needed, this section would require refinement based on data structure.

# Create and display heatmap for ecotype frequency vs. climate
print("Generating heatmap for ecotype frequency vs. climate...")
plt.figure(figsize=(12, 8))
sns.heatmap(bioexp_vs_bioecotype, cmap="viridis", cbar_label='Median Ecotype Frequency') # Add cbar label
plt.title('Median Ecotype Frequency by Experimental Site and Ecotype Origin Temperature')
plt.xlabel('Ecotype Origin Bio1')
plt.ylabel('Experimental Site Bio1')
plt.tight_layout()
plt.show()

