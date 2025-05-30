#!/usr/bin/env python
# coding: utf-8

"""
Generate Gene Tables from Significant Blocks

This script processes significant genomic blocks to identify associated genes and generate
various output files for downstream analysis. It uses BEDTools to intersect blocks with
gene annotations and creates tables linking blocks to their associated genes.

Required Files:
1. Input Files:
   - ../key_files/TAIR10_GFF3_genes_transposons_formatted4topr_w_blocks.csv: Gene annotations
   - ../key_files/blocks_snpsid_dict.pkl: Dictionary mapping SNP IDs to block IDs
   - ../key_files/TAIR10_GFF3_genes_transposons.gff: GFF file with gene annotations
   - top_hits_first_last_gen_{biovar}.csv: Significant blocks for each bioclimatic variable

2. Output Files:
   - top_hits_sign_blocks_union_first_last_gen_BH_final_{biovar}.bed: BED file of significant blocks
   - top_hits_sign_blocks_union_first_last_gen_blocks_annotated_BH_final_{biovar}.txt: Block-gene intersections
   - top_hits_sign_blocks_union_first_last_gen_gene_ids_BH_final_{biovar}.txt: List of associated gene IDs
   - top_hits_gene_and_model_sign_blocks_union_first_last_gen_BH_{biovar}.csv: Final table with block-gene associations

Script Outline:
1. Load gene annotations and block mappings
2. Process significant blocks into BED format
3. Use BEDTools to find gene intersections
4. Filter for protein-coding genes
5. Generate final tables linking blocks to genes
"""

# --- Import Libraries ---
import pandas as pd
import os 
import subprocess
import pickle
from Bio import Entrez

# Set up Entrez email for NCBI queries
Entrez.email = "tbellg@berkeley.edu"

# --- Load Required Files ---
print("Loading gene annotations and block mappings...")
tair = pd.read_csv('../key_files/TAIR10_GFF3_genes_transposons_formatted4topr_w_blocks.csv')

with open('../key_files/blocks_snpsid_dict.pkl', 'rb') as file:
    dict_blocks = pickle.load(file)

reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}

def get_gene_info_ncbi(gene_id):
    """
    Retrieve gene information from NCBI using Entrez API.
    
    Parameters:
    gene_id (str): The gene identifier to look up
    
    Returns:
    dict: Dictionary containing gene information or None if not found
    """
    try:
        search_handle = Entrez.esearch(db="gene", term=f"{gene_id}[Gene]", retmax=1)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        if search_results["IdList"]:
            gene_ncbi_id = search_results["IdList"][0]
            fetch_handle = Entrez.efetch(db="gene", id=gene_ncbi_id, retmode="xml")
            gene_data = Entrez.read(fetch_handle)
            fetch_handle.close()
            
            gene_info = gene_data[0]
            gene_name = gene_info['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']
            gene_desc = gene_info['Entrezgene_summary']
            return {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "description": gene_desc
            }
        return None
    except Exception as e:
        print(f"Error fetching data for gene {gene_id}: {e}")
        return None

def process_bioclimatic_variable(biovar):
    """
    Process significant blocks for a given bioclimatic variable.
    
    Parameters:
    biovar (str): Name of the bioclimatic variable (e.g., 'bio18')
    """
    print(f"Processing {biovar}...")
    
    # Load significant blocks
    sign_blocks_union_first_last_gen = pd.read_csv(f'top_hits_first_last_gen_{biovar}.csv')
    df = sign_blocks_union_first_last_gen.copy()
    
    # Define output file names
    bed_file_name = f"top_hits_sign_blocks_union_first_last_gen_BH_final_{biovar}.bed"
    blocks_annotated_file = f"top_hits_sign_blocks_union_first_last_gen_blocks_annotated_BH_final_{biovar}.txt"
    gene_ids_file = f'top_hits_sign_blocks_union_first_last_gen_gene_ids_BH_final_{biovar}.txt'
    
    # Create BED file
    print("Creating BED file...")
    with open(bed_file_name, "w") as bed_file:
        for block in df['block']:
            if block in dict_blocks:
                snps = dict_blocks[block]
                chrom_positions = [snp.split('_') for snp in snps]
                chrom = f"Chr{chrom_positions[0][0]}"
                start = min(int(pos) for chrom, pos in chrom_positions)
                end = max(int(pos) for chrom, pos in chrom_positions)
                bed_file.write(f"{chrom}\t{start}\t{end}\t{block}\n")
    
    # Find gene intersections using BEDTools
    print("Finding gene intersections...")
    gff_file = "../key_files/TAIR10_GFF3_genes_transposons.gff"
    cmd = f"""
    module load BEDTools/2.29.2; \
    bedtools intersect -wb -a {bed_file_name} -b {gff_file} > {blocks_annotated_file}
    """
    
    try:
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        print(f"Intersection complete. Results saved in {blocks_annotated_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running bedtools: {e}")
        return
    
    # Process annotations
    print("Processing gene annotations...")
    annot = pd.read_csv(blocks_annotated_file, sep='\t', header=None)
    annot.columns = [
        'block_chrom', 'block_start', 'block_end', 'block_id',
        'gene_chrom', 'source', 'feature', 'gene_start', 'gene_end',
        'score', 'strand', 'frame', 'attribute'
    ]
    
    # Filter for protein-coding genes
    annot = annot[annot['attribute'].str.contains('protein_coding_gene')]
    annot['gene_id'] = annot['attribute'].str.split(';').str[0].str.replace('ID=', '')
    
    # Merge with original significant blocks
    sign_blocks_union_first_last_gen = sign_blocks_union_first_last_gen.merge(
        annot, 
        left_on='block', 
        right_on='block_id'
    )
    
    # Get gene information from NCBI
    print("Retrieving gene information from NCBI...")
    gene_info_list = []
    for index, row in sign_blocks_union_first_last_gen.iterrows():
        gene_id = row['gene_id']
        block = row['block']
        print(f"Processing gene {gene_id}")
        gene_info = get_gene_info_ncbi(gene_id.strip())
        if gene_info:
            gene_info_list.append({
                "block": block,
                "gene_id": gene_info.get('gene_id'),
                "gene_name": gene_info.get('gene_name'),
                "description": gene_info.get('description')
            })
    
    # Create gene info DataFrame and merge with main results
    df_gene_info = pd.DataFrame(gene_info_list)
    if not df_gene_info.empty:
        sign_blocks_union_first_last_gen = sign_blocks_union_first_last_gen.merge(
            df_gene_info,
            on=['block', 'gene_id']
        )
    
    # Save final results
    print("Saving results...")
    sign_blocks_union_first_last_gen.to_csv(
        f'top_hits_gene_and_model_sign_blocks_union_first_last_gen_BH_{biovar}.csv'
    )
    
    # Save gene IDs
    gene_ids = annot['gene_id'].unique()
    with open(gene_ids_file, 'w') as f:
        for gene_id in gene_ids:
            f.write(f"{gene_id}\n")
    
    print(f"Processing complete for {biovar}")

def clean_up_temporary_files(biovar):
    """
    Clean up temporary files after processing.
    
    Parameters:
    biovar (str): Name of the bioclimatic variable
    """
    print(f"Cleaning up temporary files for {biovar}...")
    temp_files = [
        f"top_hits_sign_blocks_union_first_last_gen_BH_final_{biovar}.bed",
        f"top_hits_sign_blocks_union_first_last_gen_blocks_annotated_BH_final_{biovar}.txt"
    ]
    
    for file in temp_files:
        if os.path.exists(file):
            os.remove(file)
            print(f"Removed {file}")

# --- Main Execution ---
if __name__ == "__main__":
    # Process all bioclimatic variables
    bioclimatic_variables = ['bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10',
                            'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19']
    
    for biovar in bioclimatic_variables:
        process_bioclimatic_variable(biovar)
        clean_up_temporary_files(biovar)
