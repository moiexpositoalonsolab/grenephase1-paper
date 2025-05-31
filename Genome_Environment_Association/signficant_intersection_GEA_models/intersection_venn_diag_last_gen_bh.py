#!/usr/bin/env python
# coding: utf-8

"""
Analyze Intersections Between Different Models (Last Generation)

This script analyzes the overlap between significant genomic blocks identified by different models
in the last generation of the study. It includes additional models like GEMMA and performs
statistical tests to assess the significance of overlaps.

Required Files:
1. Input Files:
   - ../key_files/blocks_snpsid_dict.pkl: Dictionary mapping SNP IDs to block IDs
   - sign_blocks_union_first_last_gen_BH_final.csv: Combined results from first and last generation
   - ../wza_last_gen/gemma_231_wza_{biovar}.csv: GEMMA model results
   - ../wza_last_gen/wza_results_lfmm_{biovar}_poly7.csv: LFMM model results
   - ../wza_last_gen/wza_kendalltau_results_{biovar}_poly7.csv: Kendall Tau model results
   - ../wza_last_gen/wza_binomial_regression_{biovar}_poly7.csv: Binomial regression model results
   - ../linages_wza_picmin_last_gen/linage_based_kendall_wza_picmin.csv: Lineage-based results

2. Output Files:
   - Various statistical test results and visualizations
   - Venn diagrams showing overlaps between models

Script Outline:
1. Load and process results from different models
2. Apply multiple testing correction (BH)
3. Identify significant blocks
4. Analyze overlaps between models
5. Perform statistical tests on overlaps
6. Generate visualizations
"""

# --- Import Libraries ---
import pandas as pd
import os
import subprocess
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from venny4py.venny4py import venny4py
import pickle
import random

# --- Load Required Files ---
print("Loading block mapping dictionary...")
with open('../key_files/blocks_snpsid_dict.pkl', 'rb') as file:
    dict_blocks = pickle.load(file)

# Create reverse mapping for quick lookups
reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}

# Load combined results from first and last generation
print("Loading combined results from first and last generation...")
sign_blocks_union_first_last_gen_BH_final = pd.read_csv('sign_blocks_union_first_last_gen_BH_final.csv')

def process_model_results(biovar):
    """
    Process results from different models for a given bioclimatic variable.
    
    Parameters:
    biovar (str): Name of the bioclimatic variable (e.g., 'bio1')
    
    Returns:
    tuple: DataFrames containing significant blocks from each model
    """
    print(f"Processing results for {biovar}...")
    
    # Process GEMMA results
    path_gemma = f'../wza_last_gen/gemma_231_wza_{biovar}.csv'
    gemma_result = pd.read_csv(path_gemma)
    gemma_result = gemma_result.rename(columns={'gene': 'block'})
    gemma_result['BH_corrected_p'] = multipletests(gemma_result['Z_pVal'], method='fdr_bh')[1]
    sign_wza_gemma = gemma_result[gemma_result['BH_corrected_p'] < 0.05]
    
    # Process LFMM results
    wza_lfmm = pd.read_csv(f'../wza_last_gen/wza_results_lfmm_{biovar}_poly7.csv')
    wza_lfmm = wza_lfmm.rename(columns={'gene': 'block'})
    wza_lfmm['BH_corrected_p'] = multipletests(wza_lfmm['Z_pVal'], method='fdr_bh')[1]
    sign_wza_lfmm = wza_lfmm[wza_lfmm['BH_corrected_p'] < 0.05]
    
    # Process Kendall Tau results
    wza_kendall = pd.read_csv(f'../wza_last_gen/wza_kendalltau_results_{biovar}_poly7.csv')
    wza_kendall['BH_corrected_p'] = multipletests(wza_kendall['Z_pVal'], method='fdr_bh')[1]
    sign_wza_kendall = wza_kendall[wza_kendall['BH_corrected_p'] < 0.05]
    sign_wza_kendall = sign_wza_kendall.rename(columns={'gene': 'block'})
    
    # Process Binomial Regression results
    binomial_reg = pd.read_csv(f'../wza_last_gen/wza_binomial_regression_{biovar}_poly7.csv')
    binomial_reg['BH_corrected_p'] = multipletests(binomial_reg['Z_pVal'], method='fdr_bh')[1]
    sign_binomial_reg = binomial_reg[binomial_reg['BH_corrected_p'] < 0.05]
    sign_binomial_reg = sign_binomial_reg.rename(columns={'gene': 'block'})
    
    return sign_wza_gemma, sign_wza_lfmm, sign_wza_kendall, sign_binomial_reg

def generate_venn_diagram(sign_wza_gemma, sign_wza_lfmm, sign_wza_kendall, sign_binomial_reg, biovar):
    """
    Generate Venn diagram showing overlaps between models.
    
    Parameters:
    sign_wza_gemma, sign_wza_lfmm, sign_wza_kendall, sign_binomial_reg (DataFrame): 
        Significant blocks from each model
    biovar (str): Name of the bioclimatic variable
    """
    sets = {
        'Climate GWAS': set(sign_wza_gemma['block']),
        'LFMM GEA': set(sign_wza_lfmm['block']),
        'Kendall Tau GEA': set(sign_wza_kendall['block']),
        'Binomial GLM GEA': set(sign_binomial_reg['block'])
    }
    
    plt.figure(figsize=(10, 8))
    venny4py(sets=sets, colors=['blue', 'green', 'red', 'purple'], line_width=1.5, font_size=12)
    plt.title(f'Overlap of Significant Blocks Between Models ({biovar})')
    plt.savefig(f'venn_diagram_last_gen_{biovar}.png')
    plt.savefig(f'venn_diagram_last_gen_{biovar}.pdf')
    plt.close()

def permutation_test_3_sets(total1, total2, total3, sig1, sig2, sig3, observed_shared_123, n_permutations=10000):
    """
    Perform permutation test for overlap between three sets.
    
    Parameters:
    total1, total2, total3 (set): Total sets of elements
    sig1, sig2, sig3 (set): Significant elements in each set
    observed_shared_123 (int): Observed number of elements shared by all three sets
    n_permutations (int): Number of permutations to perform
    
    Returns:
    float: p-value for the observed overlap
    """
    count_more_extreme = 0
    
    for _ in range(n_permutations):
        # Randomly sample from total sets
        perm_sig1 = set(random.sample(list(total1), len(sig1)))
        perm_sig2 = set(random.sample(list(total2), len(sig2)))
        perm_sig3 = set(random.sample(list(total3), len(sig3)))
        
        # Count shared elements
        shared = len(perm_sig1.intersection(perm_sig2, perm_sig3))
        
        if shared >= observed_shared_123:
            count_more_extreme += 1
    
    return count_more_extreme / n_permutations

def permutation_test_2_sets(total1, total2, sig1, sig2, observed_shared_12, n_permutations=10000):
    """
    Perform permutation test for overlap between two sets.
    
    Parameters:
    total1, total2 (set): Total sets of elements
    sig1, sig2 (set): Significant elements in each set
    observed_shared_12 (int): Observed number of elements shared by both sets
    n_permutations (int): Number of permutations to perform
    
    Returns:
    float: p-value for the observed overlap
    """
    count_more_extreme = 0
    
    for _ in range(n_permutations):
        # Randomly sample from total sets
        perm_sig1 = set(random.sample(list(total1), len(sig1)))
        perm_sig2 = set(random.sample(list(total2), len(sig2)))
        
        # Count shared elements
        shared = len(perm_sig1.intersection(perm_sig2))
        
        if shared >= observed_shared_12:
            count_more_extreme += 1
    
    return count_more_extreme / n_permutations

# --- Main Execution ---
if __name__ == "__main__":
    # Process results for bio1
    biovar = 'bio1'
    sign_wza_gemma, sign_wza_lfmm, sign_wza_kendall, sign_binomial_reg = process_model_results(biovar)
    
    # Generate Venn diagram
    generate_venn_diagram(sign_wza_gemma, sign_wza_lfmm, sign_wza_kendall, sign_binomial_reg, biovar)
    
    print("Analysis complete. Results saved in venn_diagram_last_gen_{biovar}.png/pdf")


# In[1]:


import pandas as pd
import os

import subprocess
from statsmodels.stats.multitest import multipletests


# In[2]:


sign_blocks_union_first_last_gen_BH_final = pd.read_csv('sign_blocks_union_first_last_gen_BH_final.csv')


# In[3]:


sign_blocks_union_first_last_gen_BH_final[sign_blocks_union_first_last_gen_BH_final['gen'].isin([ 'first_gen,last_gen', 'last_gen'])].block.nunique()


# In[4]:


# wza_kendalltau_results_bio1.csv'
# wza_results_lfmm_bio1.csv


# In[5]:


import pickle
dict_blocks = '../key_files/blocks_snpsid_dict.pkl'

with open(dict_blocks, 'rb') as file:
    dict_blocks = pickle.load(file)

reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}


# In[ ]:





# In[6]:


## gemma 
biovar = 'bio1'


# In[7]:


path_gemma = f'../wza_last_gen/gemma_231_wza_{biovar}.csv'
        
gemma_result = pd.read_csv(path_gemma)
    
gemma_result = gemma_result.rename(columns = {'gene':'block'})
gemma_result['BH_corrected_p'] = multipletests(gemma_result['Z_pVal'], method='fdr_bh')[1]
sign_wza_gemma = gemma_result[gemma_result['BH_corrected_p'] < 0.05]


# In[ ]:





# In[8]:


## all samples wza on lfmm
wza_lfmm = pd.read_csv('../wza_last_gen/wza_results_lfmm_bio1_poly7.csv')


# In[9]:


wza_lfmm = wza_lfmm.rename(columns = {'gene':'block'})

wza_lfmm['BH_corrected_p'] = multipletests(wza_lfmm['Z_pVal'], method='fdr_bh')[1]
sign_wza_lfmm = wza_lfmm[wza_lfmm['BH_corrected_p'] < 0.05]


# In[ ]:





# In[10]:


## all samples wza on kendalltau
wza_kendall = pd.read_csv('../wza_last_gen/wza_kendalltau_results_bio1_poly7.csv')


# In[11]:


wza_kendall['BH_corrected_p'] = multipletests(wza_kendall['Z_pVal'], method='fdr_bh')[1]
sign_wza_kendall = wza_kendall[wza_kendall['BH_corrected_p'] < 0.05]
sign_wza_kendall = sign_wza_kendall.rename(columns = {'gene':'block'})


# In[ ]:





# In[12]:


genes = ['2_1265', '4_801']


# In[ ]:





# In[13]:


linages_wza_picmin = pd.read_csv('../linages_wza_picmin_last_gen//linage_based_kendall_wza_picmin.csv')


# In[14]:


th = 0.05/len(linages_wza_picmin)


# In[15]:


sign_linages_wza_picmin = linages_wza_picmin[linages_wza_picmin['p'] <= th]


# In[16]:


sign_linages_wza_picmin = sign_linages_wza_picmin.rename(columns = {'locus':'block'})


# In[ ]:





# In[ ]:





# In[17]:


binomial_reg = pd.read_csv('../wza_last_gen/wza_binomial_regression_bio1_poly7.csv')
binomial_reg['BH_corrected_p'] = multipletests(binomial_reg['Z_pVal'], method='fdr_bh')[1]
sign_binomial_reg = binomial_reg[binomial_reg['BH_corrected_p'] < 0.05]
sign_binomial_reg = sign_binomial_reg.rename(columns = {'gene':'block'})


# In[ ]:





# In[18]:


biovar = 'bio1'


# In[19]:


sign_wza_gemma


# In[ ]:





# In[20]:


sign_wza_lfmm.merge(sign_wza_gemma, on ='block')


# In[21]:


sign_wza_gemma.merge(sign_wza_kendall, on ='block')


# In[22]:


sign_wza_gemma.merge(sign_binomial_reg, on ='block')


# In[23]:


sign_wza_gemma


# In[ ]:





# In[24]:


sign_wza_lfmm.merge(sign_wza_kendall, on ='block')


# In[25]:


sign_wza_lfmm.merge(sign_binomial_reg, on ='block')


# In[ ]:





# In[26]:


sign_wza_kendall.merge(sign_binomial_reg, on ='block')


# In[27]:


all = sign_wza_lfmm.merge(sign_wza_kendall, on ='block', how = 'outer').merge(sign_binomial_reg, on ='block', how = 'outer')


# In[ ]:





# In[32]:


sign_wza_lfmm[sign_wza_lfmm['block'] == '4_2320']


# In[33]:


sign_wza_kendall[sign_wza_kendall['block'] == '4_2320']


# In[34]:


sign_binomial_reg[sign_binomial_reg['block'] == '4_2320']


# In[ ]:





# In[ ]:





# In[29]:


all[all['block'] == '4_2320']


# In[ ]:





# In[ ]:





# In[33]:


dict_blocks[i]


# In[55]:


blocks_union = sign_wza_lfmm.merge(sign_wza_kendall, on ='block', how = 'outer').merge(sign_binomial_reg, on ='block', how = 'outer')['block'].sort_values()

for i in blocks_union:
    if i == '4_2115':
        print(i)
        print(dict_blocks[i][0].split('_')[0]+ ':' + dict_blocks[i][0].split('_')[1] + '-' +  dict_blocks[i][-1].split('_')[1] )


# In[59]:


int = ['2_970', '3_2664', '5_602', '1_565', '2_1434', '3_209', '1_4950']


# In[60]:


for i in int:
    print(dict_blocks[i][0].split('_')[0]+ ':' + dict_blocks[i][0].split('_')[1] + '-' +  dict_blocks[i][0].split('_')[-1] )


# In[ ]:





# In[28]:


sign_wza_lfmm['model'] = 'wza_lfmm_l'
sign_wza_lfmm['gen'] = 'last_gen'

sign_wza_kendall['model'] = 'wza_kendall_l'
sign_wza_kendall['gen'] = 'last_gen'

sign_binomial_reg['model'] = 'wza_binom_reg_l'
sign_binomial_reg['gen'] = 'last_gen'

sign_blocks_union = pd.concat([sign_wza_lfmm[['block', 'model', 'gen']], 
          sign_wza_kendall[['block', 'model', 'gen']],
          sign_binomial_reg[['block', 'model', 'gen']]]).reset_index(drop=True)


# In[29]:


sign_blocks_union['model'] = sign_blocks_union.groupby('block')['model'].transform(lambda x: ','.join(x.unique().astype(str)))
sign_blocks_union['gen'] = sign_blocks_union.groupby('block')['gen'].transform(lambda x: ','.join(x.unique().astype(str)))

sign_blocks_union = sign_blocks_union.drop_duplicates()


# In[ ]:


sign_blocks_union_first_gen = pd.read_csv('sign_blocks_union_first_gen_BH.csv')

sign_blocks_union = pd.concat([sign_blocks_union_first_gen, sign_blocks_union])



sign_blocks_union['model'] = sign_blocks_union.groupby('block')['model'].transform(lambda x: ','.join(x.unique().astype(str)))
sign_blocks_union['gen'] = sign_blocks_union.groupby('block')['gen'].transform(lambda x: ','.join(x.unique().astype(str)))


# In[35]:


sign_blocks_union


# In[36]:


sign_blocks_union = sign_blocks_union.drop_duplicates()


# In[37]:


sign_blocks_union


# In[38]:


sign_blocks_union.to_csv('sign_blocks_union_first_last_gen_BH.csv', index=None)


# In[ ]:





# In[ ]:





# In[ ]:





# In[41]:


import matplotlib
from matplotlib import pyplot as plt
from upsetplot import generate_counts, plot


# In[ ]:





# In[107]:


# Extract 'block' columns and convert to sets (excluding set1/picmin)
set2 = set(sign_wza_lfmm['block'])
set3 = set(sign_wza_kendall['block'])
set4 = set(sign_binomial_reg['block'])

# Calculate the unions using set operations
counts = {
    (True, False, False): len(set2) ,  # Only in set2
    (False, True, False): len(set3),  # Only in set3
    (False, False, True): len(set4),  # Only in set4
    (True, True, False): len(set2.intersection(set3)),  # In set2 and set3 only
    (True, False, True): len(set2.intersection(set4)),  # In set2 and set4 only
    (False, True, True): len(set3.intersection(set4)),  # In set3 and set4 only
    (True, True, True): len(set2.intersection(set3).intersection(set4))  # In set2, set3, and set4
}

# Convert the dictionary to a pandas Series with a MultiIndex
index = pd.MultiIndex.from_tuples(counts.keys(), names=['LFMM+WZA', 'KendallTau+WZA', 'BinomialReg+WZA'])
data = pd.Series(counts, index=index)

# Plot the data
plot(data, show_counts=True)
plt.suptitle("GEA significant blocks Terminal Gen (Reduced to 3 Sets)")
plt.show()


# In[52]:


set2.intersection(set3)


# In[53]:


set2.intersection(set4)


# In[54]:


set3.intersection(set4)


# In[79]:


set1


# In[117]:


from venny4py.venny4py import *


# In[ ]:





# In[119]:





# In[ ]:





# In[ ]:


sets = {


# In[ ]:


set1 = set(sign_wza_gemma['block'])
set2 = set(sign_wza_lfmm['block'])
set3 = set(sign_wza_kendall['block'])
set4 = set(sign_binomial_reg['block'])


# In[124]:


from venny4py.venny4py import venny4py

# Define the sets
sets = {
    'sign_wza_gemma': set(sign_wza_gemma['block']), 
    'sign_wza_lfmm': set(sign_wza_lfmm['block']), 
    'sign_wza_kendall': set(sign_wza_kendall['block']), 
    'sign_binomial_reg': set(sign_binomial_reg['block'])
}

# Use venny4py to visualize the overlap
venny4py(sets=sets)


# In[126]:


set(sign_wza_gemma['block']).intersection(set(sign_wza_lfmm['block']))


# In[ ]:





# In[116]:


plot(data, show_counts=True)
plt.suptitle("GEA significant blocks Terminal Gen")
plt.show()


# In[ ]:





# In[75]:


index = pd.MultiIndex.from_tuples(counts.keys(), names=['GEMMA+WZA', 'LFMM+WZA', 'KendallTau+WZA', 'BinomialReg+WZA'])
data = pd.Series(counts, index=index)


# In[76]:


data


# In[78]:


import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import plot

# Define your counts as given
counts = {
    (True, False, False, False): 23,  # Only in set1
    (False, True, False, False): 15,  # Only in set2
    (False, False, True, False): 15,  # Only in set3
    (False, False, False, True): 19,  # Only in set4
    (True, True, False, False): 3,    # In set1 and set2 only
    (True, False, True, False): 2,    # In set1 and set3 only
    (True, False, False, True): 1,    # In set1 and set4 only
    (False, True, True, False): 1,    # In set2 and set3 only
    (False, True, False, True): 1,    # In set2 and set4 only
    (False, False, True, True): 4,    # In set3 and set4 only
    (True, True, True, True): 1       # In all sets
}

# Convert the dictionary to a pandas Series with a MultiIndex
index = pd.MultiIndex.from_tuples(counts.keys(), names=['GEMMA+WZA', 'LFMM+WZA', 'KendallTau+WZA', 'BinomialReg+WZA'])
data = pd.Series(counts, index=index)

# Manually sort the MultiIndex based on your specific order
sorted_order = [
    (True, False, False, False),
    (False, True, False, False),
    (False, False, True, False),
    (False, False, False, True),
    (True, True, False, False),
    (True, False, True, False),
    (True, False, False, True),
    (False, True, True, False),
    (False, True, False, True),
    (False, False, True, True),
    (True, True, True, True)
]

# Reorder the data Series according to sorted_order
data = data.reindex(sorted_order)

# Plot using upsetplot
plot(data, show_counts=True)
plt.suptitle("GEA significant blocks Terminal Gen (Correct Order)")
plt.show()


# In[71]:


import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import from_contents, plot

# Define the sets
set1 = set(sign_wza_gemma['block'])
set2 = set(sign_wza_lfmm['block'])
set3 = set(sign_wza_kendall['block'])
set4 = set(sign_binomial_reg['block'])

# Generate the data for the UpSet plot
data = {
    'GEMMA+WZA': set1,
    'LFMM+WZA': set2,
    'KendallTau+WZA': set3,
    'BinomialReg+WZA': set4
}

# Prepare the counts data using from_contents function
counts = from_contents(data)

# Define the specific combinations you want to show
desired_combinations = {
    (True, False, False, False),  # Only in set1
    (False, True, False, False),  # Only in set2
    (False, False, True, False),  # Only in set3
    (False, False, False, True),  # Only in set4
    (True, True, False, False),   # In set1 and set2 only
    (True, False, True, False),   # In set1 and set3 only
    (True, False, False, True),   # In set1 and set4 only
    (False, True, True, False),   # In set2 and set3 only
    (False, True, False, True),   # In set2 and set4 only
    (False, False, True, True),   # In set3 and set4 only
    (True, True, True, True)      # In all sets
}

# Filter the counts to only include the desired combinations
filtered_counts = counts.loc[desired_combinations]

# Plot using upsetplot
plot(filtered_counts, show_counts=True)
plt.suptitle("GEA significant blocks Terminal Gen (Selected Intersections)")
plt.show()


# In[65]:


counts


# In[ ]:





# In[63]:


len(set1.intersection(set2))


# In[48]:


set1.intersection(set4)


# In[22]:


#sign_linages_wza_picmin = sign_linages_wza_picmin[sign_linages_wza_picmin['n_est'] > 3]


# In[28]:


# Extract 'block' columns and convert to sets (excluding set1/picmin)
set2 = set(sign_wza_lfmm['block'])
set3 = set(sign_wza_kendall['block'])
set4 = set(sign_binomial_reg['block'])

# Calculate the unions using set operations
counts = {
    (True, False, False): len(set2) ,  # Only in set2
    (False, True, False): len(set3),  # Only in set3
    (False, False, True): len(set4),  # Only in set4
    (True, True, False): len(set2.intersection(set3)),  # In set2 and set3 only
    (True, False, True): len(set2.intersection(set4)),  # In set2 and set4 only
    (False, True, True): len(set3.intersection(set4)),  # In set3 and set4 only
    (True, True, True): len(set2.intersection(set3).intersection(set4))  # In set2, set3, and set4
}

# Convert the dictionary to a pandas Series with a MultiIndex
index = pd.MultiIndex.from_tuples(counts.keys(), names=['LFMM+WZA', 'KendallTau+WZA', 'BinomialReg+WZA'])
data = pd.Series(counts, index=index)

# Plot the data
plot(data, show_counts=True)
plt.suptitle("GEA significant blocks Terminal Gen (Reduced to 3 Sets)")
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[38]:


data.reset_index


# In[52]:


from matplotlib import pyplot as plt

from upsetplot import generate_counts, plot, UpSet

greyscale_colors = ['#666666', '#BBBBBB', '#666666', '#BBBBBB', '#666666', '#BBBBBB']

# Define green tones for significant points (alternating dark and light greens)
green_colors = ['#2aad2a', '#208420'] # 

# Create the UpSet plot using the UpSet class
upset = UpSet(data, facecolor="#666666",show_counts=True)

# Style subsets with specific size

upset.style_subsets(min_degree=19, max_degree=19, facecolor="green")
upset.style_subsets(min_degree=15, max_degree=15, facecolor="purple")
upset.style_subsets(min_degree=11, max_degree=11, facecolor="orange")

upset.style_subsets(min_degree=2, facecolor="#6DAEDB")
upset.style_subsets(min_degree=3, facecolor="#1D70A2")

# Plot the upset plot
plot_result = upset.plot()

plot_result["intersections"].set_ylabel("Significant blocks")
plot_result["totals"].set_xlabel("Category size")
plt.suptitle("GEA significant blocks Terminal Gen")

plt.savefig("gea_significant_blocks_terminal_gen_bh.png",bbox_inches='tight')  # Save as PNG with 300 DPI resolution
plt.savefig("gea_significant_blocks_terminal_gen_bh.svg",bbox_inches='tight')  # Save as PNG with 300 DPI resolution

# Show the plot
plt.show()


# In[ ]:





# In[ ]:


## permutation tests 


# In[ ]:


import numpy as np


# In[ ]:


# Function to perform permutation test for overlap across all three sets
def permutation_test_3_sets(total1, total2, total3, sig1, sig2, sig3, observed_shared_123, n_permutations=10000):
    overlap_count = 0
    print(total1, total2, total3)
    print(sig1, sig2, sig3)
    print(observed_shared_123)
    # Perform permutations
    for _ in range(n_permutations):
        # Randomly select significant SNPs for each set
        perm_sig1 = np.random.choice(total1, sig1, replace=False)
        perm_sig2 = np.random.choice(total2, sig2, replace=False)
        perm_sig3 = np.random.choice(total3, sig3, replace=False)

        # Find the intersection (shared SNPs) across all three sets
        shared_123 = len(set(perm_sig1) & set(perm_sig2) & set(perm_sig3))

        # Count how many times the shared SNPs meet or exceed the observed value
        if shared_123 >= observed_shared_123:
            overlap_count += 1

    # Calculate empirical p-value
    p_value = overlap_count / n_permutations
    return p_value


# In[32]:


# Extract 'block' columns and convert to sets
set1 = set(sign_linages_wza_picmin['block'])
set2 = set(sign_wza_lfmm['block'])
set3 = set(sign_wza_kendall['block'])


# In[33]:


permutation_test_3_sets(len(wza_kendall),
                        len(wza_lfmm),
                        len(linages_wza_picmin),
                        len(sign_wza_kendall),
                        len(sign_wza_lfmm),
                        len(sign_linages_wza_picmin),
                        len(set1.intersection(set2).intersection(set3)),
                        1000000
)


# In[ ]:


# Function to perform permutation test for overlap across all three sets
def permutation_test_2_sets(total1, total2, sig1, sig2, observed_shared_12, n_permutations=10000):
    overlap_count = 0
    print(total1, total2)
    print(sig1, sig2)
    print(observed_shared_12)
    # Perform permutations
    for _ in range(n_permutations):
        # Randomly select significant SNPs for each set
        perm_sig1 = np.random.choice(total1, sig1, replace=False)
        perm_sig2 = np.random.choice(total2, sig2, replace=False)

        # Find the intersection (shared SNPs) across all three sets
        shared_12 = len(set(perm_sig1) & set(perm_sig2))

        # Count how many times the shared SNPs meet or exceed the observed value
        if shared_12 >= observed_shared_12:
            overlap_count += 1

    # Calculate empirical p-value
    p_value = overlap_count / n_permutations
    return p_value


# In[ ]:


len(set1.intersection(set2))


# In[ ]:


permutation_test_2_sets(len(wza_kendall),
                        len(wza_lfmm),
                        len(sign_wza_kendall),
                        len(sign_wza_lfmm),
                        len(set2.intersection(set3)),
                        1000000
)


# In[ ]:


permutation_test_2_sets(len(wza_kendall),
                        len(linages_wza_picmin),
                        len(sign_wza_kendall),
                        len(sign_linages_wza_picmin),
                        len(set1.intersection(set3)),
                        1000000
)


# In[ ]:


permutation_test_2_sets(len(wza_lfmm),
                        len(linages_wza_picmin),
                        len(sign_wza_lfmm),
                        len(sign_linages_wza_picmin),
                        len(set1.intersection(set2)),
                        1000000
)


# In[ ]:





# In[28]:


shared_blocks = ['2_1265', '4_801']


# In[42]:


selection_atlas = pd.read_csv('../key_files/POP_EVOLUTION_selection_atlas_haplotype_sig_counts.txt',sep='\t')


# In[43]:


selection_atlas


# In[ ]:





# In[ ]:





# In[44]:


selection_atlas['block'] = selection_atlas['chr'].astype(str) + '_' + selection_atlas['start'].astype(str) 


# In[45]:


bb


# In[50]:


# Create a new column 'chr_rowindex' that combines 'chr' and a zero-based index for each 'chr' group
selection_atlas['chr_rowindex'] = selection_atlas.groupby('chr').cumcount()

# Create an ID that combines 'chr' and the 'chr_rowindex'
selection_atlas['block_id'] = selection_atlas['chr'].astype(str) + "_" + selection_atlas['chr_rowindex'].astype(str)


# In[48]:


selection_atlas


# In[ ]:





# In[51]:


selection_atlas[selection_atlas['block_id'].isin(shared_blocks)]


# In[ ]:





# In[52]:


dict_blocks['2_1265']


# In[35]:


selection_atlas


# In[ ]:




