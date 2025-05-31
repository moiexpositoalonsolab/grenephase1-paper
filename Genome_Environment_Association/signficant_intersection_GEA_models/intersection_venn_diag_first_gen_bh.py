#!/usr/bin/env python
# coding: utf-8

"""
Analyze Intersections Between Different Models (First Generation)

This script analyzes the overlap between significant genomic blocks identified by different models
in the first generation of the study. It performs statistical tests to assess the significance
of these overlaps and generates visualizations.

Required Files:
1. Input Files:
   - ../key_files/blocks_snpsid_dict.pkl: Dictionary mapping SNP IDs to block IDs
   - ../wza/wza_results_lfmm_bio1_poly7.csv: LFMM model results
   - ../wza/wza_kendalltau_results_bio1_poly7.csv: Kendall Tau model results
   - ../wza/wza_binomial_regression_bio1_poly7.csv: Binomial regression model results

2. Output Files:
   - sign_blocks_union_first_gen_BH.csv: Combined significant blocks from all models
   - Various statistical test results and visualizations

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

# In[2]:


# wza_kendalltau_results_bio1.csv'
# wza_results_lfmm_bio1.csv


# In[3]:


dict_blocks = '../key_files/blocks_snpsid_dict.pkl'

with open(dict_blocks, 'rb') as file:
    dict_blocks = pickle.load(file)

reverse_mapping = {item: key for key, values in dict_blocks.items() for item in values}


# In[4]:


#dict_blocks['2_1265']#[-1] # 2:9540996-9714921   2_11533904  2_11534263


# In[5]:


## all samples wza on lfmm
top_candidate = pd.read_csv('../wza/wza_results_lfmm_bio1_poly7.csv')


# In[6]:


top_candidate = top_candidate.rename(columns = {'gene':'block'})


# In[7]:


sig_top_candidate = top_candidate[top_candidate['top_candidate_p'] <= 0.05/len(top_candidate)]


# In[ ]:





# In[8]:


## all samples wza on lfmm
wza_lfmm = pd.read_csv('../wza/wza_results_lfmm_bio1_poly7.csv')


# In[9]:


wza_lfmm


# In[9]:


wza_lfmm = wza_lfmm.rename(columns = {'gene':'block'})


# In[10]:


wza_lfmm['BH_corrected_p'] = multipletests(wza_lfmm['Z_pVal'], method='fdr_bh')[1]
sign_wza_lfmm = wza_lfmm[wza_lfmm['BH_corrected_p'] < 0.05]


# In[11]:


## all samples wza on kendalltau
wza_kendall = pd.read_csv('../wza/wza_kendalltau_results_bio1_poly7.csv')


# In[12]:


wza_kendall['BH_corrected_p'] = multipletests(wza_kendall['Z_pVal'], method='fdr_bh')[1]
sign_wza_kendall = wza_kendall[wza_kendall['BH_corrected_p'] < 0.05]
sign_wza_kendall = sign_wza_kendall.rename(columns = {'gene':'block'})


# In[13]:


#linages_wza_picmin = pd.read_csv('../linages_wza_picmin/linage_based_kendall_wza_picmin.csv')


# In[14]:


#linages_wza_picmin[linages_wza_picmin['n_est'] ==12].sort_values('p')


# In[15]:


#th = 0.05/len(linages_wza_picmin)


# In[16]:


#sign_linages_wza_picmin = linages_wza_picmin[linages_wza_picmin['p'] <= th]


# In[17]:


#sign_linages_wza_picmin = sign_linages_wza_picmin.rename(columns = {'locus':'block'})


# In[ ]:





# In[18]:


binomial_reg = pd.read_csv('../wza/wza_binomial_regression_bio1_poly7.csv')
binomial_reg['BH_corrected_p'] = multipletests(binomial_reg['Z_pVal'], method='fdr_bh')[1]
sign_binomial_reg = binomial_reg[binomial_reg['BH_corrected_p'] < 0.05]
sign_binomial_reg = sign_binomial_reg.rename(columns = {'gene':'block'})


# In[19]:


sign_binomial_reg['block'].values


# In[20]:


sign_wza_lfmm.merge(sign_wza_kendall, on ='block')


# In[21]:


sign_wza_lfmm.merge(sign_binomial_reg, on ='block')


# In[22]:


sign_wza_kendall.merge(sign_binomial_reg, on ='block')


# In[23]:


sign_wza_lfmm.merge(sign_wza_kendall, on ='block').merge(sign_binomial_reg,on ='block')


# In[ ]:


2_199


# In[33]:


blocks_union = sign_wza_lfmm.merge(sign_wza_kendall, on ='block', how = 'outer').merge(sign_binomial_reg, on ='block', how = 'outer')['block']

for i in blocks_union:
    print(i)
    print(dict_blocks[i][0].split('_')[0]+ ':' + dict_blocks[i][0].split('_')[1] + '-' +  dict_blocks[i][0].split('_')[-1] )


# In[25]:


sign_wza_lfmm.merge(sign_wza_kendall, on ='block', how = 'outer').merge(sign_binomial_reg, on ='block', how = 'outer')


# In[26]:


sign_wza_lfmm['model'] = 'wza_lfmm_f'
sign_wza_lfmm['gen'] = 'first_gen'

sign_wza_kendall['model'] = 'wza_kendall_f'
sign_wza_kendall['gen'] = 'first_gen'

sign_binomial_reg['model'] = 'wza_binom_reg_f'
sign_binomial_reg['gen'] = 'first_gen'

sign_blocks_union = pd.concat([sign_wza_lfmm[['block', 'model', 'gen']], 
          sign_wza_kendall[['block', 'model', 'gen']],
          sign_binomial_reg[['block', 'model', 'gen']]]).reset_index(drop=True)


# In[27]:


sign_blocks_union['model'] = sign_blocks_union.groupby('block')['model'].transform(lambda x: ','.join(x.unique().astype(str)))
sign_blocks_union['gen'] = sign_blocks_union.groupby('block')['gen'].transform(lambda x: ','.join(x.unique().astype(str)))

sign_blocks_union = sign_blocks_union.drop_duplicates()


# In[28]:


sign_blocks_union


# In[29]:


sign_blocks_union.to_csv('sign_blocks_union_first_gen_BH.csv', index=None)


# In[ ]:





# In[30]:


import matplotlib
from matplotlib import pyplot as plt

from upsetplot import generate_counts, plot


# In[31]:


#sign_linages_wza_picmin = sign_linages_wza_picmin[sign_linages_wza_picmin['n_est'] > 3]


# In[ ]:





# In[32]:


# Extract 'block' columns and convert to sets (excluding set1/picmin)
set2 = set(sign_wza_lfmm['block'])
set3 = set(sign_wza_kendall['block'])
set4 = set(sign_binomial_reg['block'])

# Calculate the unions using set operations
counts = {
    (True, False, False): len(set2),  # Only in set2
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


# In[28]:


sign_wza_kendall


# In[ ]:





# In[ ]:





# In[34]:


from matplotlib import pyplot as plt

from upsetplot import generate_counts, plot, UpSet


# In[35]:


greyscale_colors = ['#666666', '#BBBBBB', '#666666', '#BBBBBB', '#666666', '#BBBBBB']

# Define green tones for significant points (alternating dark and light greens)
green_colors = ['#2aad2a', '#208420'] # 


# In[36]:


# Create the UpSet plot using the UpSet class
upset = UpSet(data, facecolor="#666666",show_counts=True)

# Style subsets with specific size
upset.style_subsets(min_degree=3, facecolor="#2aad2a")


# Plot the upset plot
plot_result = upset.plot()

plot_result["intersections"].set_ylabel("Significant blocks")
plot_result["totals"].set_xlabel("Category size")
plt.suptitle("GEA significant blocks First Generation")

plt.savefig("gea_significant_blocks.svg",bbox_inches='tight')  # Save as PNG with 300 DPI resolution

# Show the plot
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[34]:


## permutation test 


# In[35]:


import numpy as np


# In[36]:


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


# In[37]:


# Extract 'block' columns and convert to sets
set1 = set(sign_linages_wza_picmin['block'])
set2 = set(sign_wza_lfmm['block'])
set3 = set(sign_wza_kendall['block'])


# In[ ]:


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





# In[ ]:





# In[ ]:





# In[42]:


selection_atlas = pd.read_csv('../key_files/POP_EVOLUTION_selection_atlas_haplotype_sig_counts.txt',sep='\t')


# In[43]:


selection_atlas


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




