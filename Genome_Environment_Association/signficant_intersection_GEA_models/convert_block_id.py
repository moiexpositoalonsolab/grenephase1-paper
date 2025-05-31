#!/usr/bin/env python
# coding: utf-8

"""
Convert Block IDs to Readable Format

This script converts block IDs from the format 'chrom_block' (e.g., '2_841') to a more readable
format that includes the genomic coordinates (e.g., '2:123456-123789'). This is useful for
visualization and interpretation of genomic regions.

Required Files:
1. Input Files:
   - ../key_files/blocks_snpsid_dict.pkl: Dictionary mapping SNP IDs to block IDs
   - genes_info_BH_tair10_{biovar}.csv: Gene information file with block IDs

2. Output Files:
   - genes_info_BH_tair10_{biovar}_w_block_id.csv: Updated gene information file with formatted block IDs

Script Outline:
1. Load block mapping dictionary
2. Load gene information file
3. Convert block IDs to genomic coordinate format
4. Save updated file
"""

# --- Import Libraries ---
import pandas as pd
import pickle

# --- Load Required Files ---
print("Loading block mapping dictionary...")
with open('../key_files/blocks_snpsid_dict.pkl', 'rb') as file:
    dict_blocks = pickle.load(file)

# Create reverse mapping for quick lookups
reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}

def format_block_id(block_id):
    """
    Convert a block ID to genomic coordinate format.
    
    Parameters:
    block_id (str): Block ID in format 'chrom_block' (e.g., '2_841')
    
    Returns:
    str: Formatted block ID with genomic coordinates (e.g., '2:123456-123789')
    """
    try:
        # Get first and last SNP positions for the block
        first_snp = dict_blocks[block_id][0].split('_')
        last_snp = dict_blocks[block_id][-1].split('_')
        
        # Format as chromosome:start-end
        return f"{first_snp[0]}:{first_snp[-1]}-{last_snp[-1]}"
    except KeyError:
        print(f"Warning: Block ID {block_id} not found in dictionary")
        return block_id
    except Exception as e:
        print(f"Error formatting block ID {block_id}: {e}")
        return block_id

def process_bioclimatic_variable(biovar):
    """
    Process gene information file for a given bioclimatic variable.
    
    Parameters:
    biovar (str): Name of the bioclimatic variable (e.g., 'bio18')
    """
    print(f"Processing {biovar}...")
    
    # Load gene information file
    input_file = f'genes_info_BH_tair10_{biovar}.csv'
    output_file = f'genes_info_BH_tair10_{biovar}_w_block_id.csv'
    
    try:
        # Read the input file
        df = pd.read_csv(input_file)
        
        # Convert block IDs
        print("Converting block IDs...")
        df['block_id'] = df['block_id'].apply(format_block_id)
        
        # Save the updated file
        print(f"Saving results to {output_file}...")
        df.to_csv(output_file, index=None)
        print(f"Processing complete for {biovar}")
        
    except FileNotFoundError:
        print(f"Error: Input file {input_file} not found")
    except Exception as e:
        print(f"Error processing {biovar}: {e}")

# --- Main Execution ---
if __name__ == "__main__":
    # Example usage with bio18
    process_bioclimatic_variable('bio18')


# In[2]:


import pickle
dict_blocks = '../key_files/blocks_snpsid_dict.pkl'

with open(dict_blocks, 'rb') as file:
    dict_blocks = pickle.load(file)

reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}


# In[20]:


bio18 = pd.read_csv('genes_info_BH_tair10_bio18.csv')


# In[21]:


dict_blocks['2_841'][0]


# In[22]:


dict_blocks['2_841'][-1]


# In[23]:


dict_blocks[i][0]


# In[24]:


import pandas as pd

# Assuming bio18 and dict_blocks are already defined
# Function to format the string based on block_id
def format_block_id(block_id):
    # Accessing the first and last elements of the list for the block_id in dict_blocks
    first_part = dict_blocks[block_id][0].split('_')
    last_part = dict_blocks[block_id][-1].split('_')
    # Formatting the string as per your provided structure
    return f"{first_part[0]}:{first_part[-1]}-{last_part[-1]}"

# Apply the function to the block_id column and create a new column
bio18['block_id'] = bio18['block_id'].apply(format_block_id)


# In[26]:


bio18.to_csv('genes_info_BH_tair10_bio18_w_block_id.csv',index=None)


# In[ ]:




