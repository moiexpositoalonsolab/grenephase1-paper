#!/usr/bin/env python
# coding: utf-8

"""
This script prepares data for the last generation analysis.
It processes allele frequency data, generates allele counts, and creates partitions for downstream analysis.

Required input files:
- key_files/merged_sample_table.csv: Sample information including flower counts
- key_files/final_gen.csv: List of last generation sample names
- key_files/merged_hapFIRE_allele_frequency_LDpruned.txt: LD-pruned allele frequencies
- key_files/greneNet_final_v1.1_LDpruned.recode.vcf: LD-pruned VCF file for SNP information

Output:
- Allele counts for last generation samples
- Partitions of SNPs for analysis
- SNP names and positions
"""

import pandas as pd
import os
import allel
import numpy as np
from dask import delayed
import dask.dataframe as dd
from typing import Dict, List, Union
import random

# Set up paths
KEY_FILES = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/key_files/'
OUTPUT_DIR = 'last_gen_analysis/'

def makedir(directory: str) -> str:
    """
    Create directory if it doesn't exist.
    
    Args:
        directory: Directory path to create
        
    Returns:
        Created directory path
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory

def read_node_list(filename: str) -> List[str]:
    """
    Read a list of nodes from a file.
    
    Args:
        filename: Path to file containing node list
        
    Returns:
        List of node names
    """
    with open(filename, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def generate_allele_counts(df: pd.DataFrame, num_flowers_map: Dict[str, int]) -> pd.DataFrame:
    """
    Generate allele counts from allele frequencies.
    
    Args:
        df: DataFrame containing allele frequencies
        num_flowers_map: Dictionary mapping sample names to flower counts
        
    Returns:
        DataFrame containing allele counts
    """
    allele_counts = {}
    for sample in df.columns:
        num_flowers = num_flowers_map[sample]
        # Calculate minor allele counts
        allele_counts[f'{sample}_minor'] = df[sample] * num_flowers * 2
        # Calculate major allele counts
        maj = 1 - df[sample]
        allele_counts[f'{sample}_major'] = maj * num_flowers * 2
    
    return pd.concat(allele_counts, axis=1).round().astype(int)

def save_partition(partition: pd.DataFrame, partition_index: int, output_path: str):
    """
    Save a partition to files with and without header.
    
    Args:
        partition: DataFrame containing partition data
        partition_index: Index of the partition
        output_path: Directory to save the files
    """
    # Create output directory if it doesn't exist
    makedir(output_path)
    
    # Save with and without header
    partition.to_csv(f'{output_path}/partition_{partition_index}_wheader.txt', 
                    sep='\t', index=False)
    partition.to_csv(f'{output_path}/partition_{partition_index}_nheader.txt', 
                    sep='\t', index=False, header=False)

def main():
    # Create output directory
    makedir(OUTPUT_DIR)
    
    # Load sample information
    flowers = pd.read_csv(KEY_FILES + 'merged_sample_table.csv')
    last_gen = pd.read_csv(KEY_FILES + 'final_gen.csv')['sample_name']
    
    # Filter flowers data for last generation
    flowers = flowers[flowers['sample_name'].isin(last_gen)]
    num_flowers_map = flowers.set_index('sample_name')['total_flower_counts'].to_dict()
    
    # Load allele frequency data
    allele_freq_ldp = KEY_FILES + 'merged_hapFIRE_allele_frequency_LDpruned.txt'
    df_ldp = dd.read_csv(allele_freq_ldp, sep='\t', usecols=last_gen)
    
    # Load VCF file for SNP information
    ld_prunned_vcf_file = KEY_FILES + 'greneNet_final_v1.1_LDpruned.recode.vcf'
    ld_prunned_vcf = allel.read_vcf(ld_prunned_vcf_file)
    
    # Get chromosome and position information
    ld_prunned_chrom = ld_prunned_vcf['variants/CHROM']
    ld_prunned_pos = ld_prunned_vcf['variants/POS']
    
    # Convert dask dataframe to pandas
    df_ldp = df_ldp.compute()
    
    # Generate allele counts
    allele_counts = generate_allele_counts(df_ldp, num_flowers_map)
    
    # Save allele counts
    allele_counts.to_csv(f'{OUTPUT_DIR}/allele_counts_lastgen_wheader.txt', 
                        sep='\t', index=False)
    allele_counts.to_csv(f'{OUTPUT_DIR}/allele_counts_lastgen_nheader.txt', 
                        sep='\t', index=False, header=False)
    
    # Create SNP names file
    snp_names = pd.DataFrame({
        'chrom': ld_prunned_chrom,
        'pos': ld_prunned_pos,
        'snp_id': [f"{chrom}_{pos}" for chrom, pos in zip(ld_prunned_chrom, ld_prunned_pos)]
    })
    
    # Save SNP names
    snp_names.to_csv(f'{OUTPUT_DIR}/snp_names.txt', sep='\t', index=False)
    
    # Create partitions (example: 10 partitions)
    num_partitions = 10
    partition_size = len(df_ldp) // num_partitions
    
    for i in range(num_partitions):
        start_idx = i * partition_size
        end_idx = start_idx + partition_size if i < num_partitions - 1 else len(df_ldp)
        
        partition = df_ldp.iloc[start_idx:end_idx]
        save_partition(partition, i + 1, OUTPUT_DIR)
        
        # Save corresponding SNP names for this partition
        partition_snps = snp_names.iloc[start_idx:end_idx]
        partition_snps.to_csv(f'{OUTPUT_DIR}/snp_names_partition_{i+1}.txt', 
                            sep='\t', index=False)

if __name__ == "__main__":
    main() 