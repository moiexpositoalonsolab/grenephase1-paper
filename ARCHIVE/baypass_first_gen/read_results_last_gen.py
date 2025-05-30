#!/usr/bin/env python
# coding: utf-8

# In[1]:


import subprocess
import pandas as pd
import os
import pickle
import seaborn as sns
#base_directory = '/central/scratch/tbellagi/gea/baypass_terminal/cmd_files/shfiles/'


# In[6]:


base_directory = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/results_final_gen/'


# In[7]:


##the bfffiles have the bayes factor for each of the covariates for each fo the markers 
bffiles = [i for i in os.listdir('results_final_gen/') if 'betai_reg.out' in i]
#xtxfiles = [i for i in os.listdir('results') if 'pi_xtx.out' in i] 


# In[8]:


## caltech 
#base_directory = '/central/scratch/tbellagi/gea/baypass_terminal/cmd_files/shfiles/'


# In[34]:


partitions = set(pd.Series(bffiles).str.split('_').str[1].astype(int))


# In[ ]:


## i creates 173 partitions 


# In[37]:


list(partitions)[-1]


# In[9]:


env_vars_names = ['bio1', 'bio17']


# In[50]:


lengths = []
for i in range(173):
    length = len(pd.read_csv(f'results_final_gen/partition_{i}_chain_2_summary_betai_reg.out', sep='\s+'))
    lengths.append(length)


# In[51]:


sum(lengths)


# In[47]:


1287906


# In[42]:


pd.read_csv('results_final_gen/partition_173_chain_3_summary_betai_reg.out', sep='\s+')


# In[43]:


pd.read_csv('results_final_gen/partition_130_chain_3_summary_betai_reg.out', sep='\s+')


# In[ ]:





# In[ ]:





# In[10]:


#base_directory = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/cmd_files/shfiles/'
##the bfffiles have the bayes factor for each of the covariates for each fo the markers 
#bffiles = [i for i in os.listdir('cmd_files/shfiles/') if 'betai_reg.out' in i]
#xtxfiles = [i for i in os.listdir('results') if 'pi_xtx.out' in i]

len(bffiles)

def extract_partitions(bffiles):
    partitions = {}
    for file_name in bffiles:
        parts = file_name.split('_')
        partition = parts[1]
        chain = parts[3]
        if partition not in partitions:
            partitions[partition] = []
        partitions[partition].append(chain)
    return partitions

file_pattern = 'partition_{partition}_chain_{chain}_summary_betai_reg.out'

def process_file(partition, chain, output_dir):
    bf_file_name = file_pattern.format(partition=partition, chain=chain)
    bf_file = pd.read_csv(base_directory + bf_file_name, sep='\s+')
    
    # Load loci and environment variable names
    partition_name = 'partition_{}'.format(partition)
    with open(f'individual_gfiles/loci_{partition_name}', 'rb') as f:
        gloci = pickle.load(f)
    gdict = {num+1: locus for num, locus in enumerate(gloci)}
    bf_file['MRK'] = bf_file['MRK'].map(gdict)
    
    with open('env_vars_names', 'rb') as f:
        env_vars_names = pickle.load(f)
    env_vars_dict = {num+1: env_var for num, env_var in enumerate(env_vars_names)}
    bf_file['COVARIABLE'] = bf_file['COVARIABLE'].map(env_vars_dict)

    # Rename columns to include chain number, except for 'COVARIABLE' and 'MRK'
    bf_file = bf_file.rename(columns={col: f"{col}_chain{chain}" for col in bf_file.columns if col not in ['COVARIABLE', 'MRK']})

    # Write to file by COVARIABLE
    for covariable, df_group in bf_file.groupby('COVARIABLE'):
        output_path = os.path.join(output_dir, f"combined_covariable_{covariable}_chain{chain}.csv")
        df_group.drop(columns='COVARIABLE').to_csv(output_path, mode='a', index=False, header=not os.path.exists(output_path))


partitions = extract_partitions(bffiles)  # Assuming bffiles is defined somewhere

all_part_collected = list(partitions.keys())

output_directory = 'output_last_gen/'
os.makedirs(output_directory, exist_ok=True)

for partition, chains in partitions.items():
    for chain in chains:
        process_file(partition, chain, output_directory)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[12]:


chain1_bio1 = pd.read_csv('output_last_gen/combined_covariable_bio1_chain1.csv')#['chain'].unique()


# In[13]:


len(chain1_bio1)


# In[14]:


chain2_bio1 = pd.read_csv('output_last_gen/combined_covariable_bio1_chain2.csv')#['chain'].unique()


# In[15]:


len(chain2_bio1)


# In[16]:


chain3_bio1 = pd.read_csv('output_last_gen/combined_covariable_bio1_chain3.csv')#['chain'].unique()


# In[17]:


len(chain3_bio1)


# In[22]:


bio1_results = chain1_bio1.merge(chain2_bio1, on = 'MRK').merge(chain3_bio1, on = 'MRK')


# In[23]:


bio1_results


# In[20]:


chain1_bio1[chain1_bio1['BF(dB)_chain1'] > 0 ]


# In[ ]:





# In[26]:


def rank_baypass(df):
    """Create columns to rank loci within and across chains."""
    import pandas
    import tqdm
    import math
    chains = ['chain1', 'chain2', 'chain3']  # Your chains are named 'chain1', 'chain2', 'chain3'

    # rank per chain based on chain BF
    for chain in tqdm.tqdm(chains, desc='chains'):
        col = f'BF(dB)_{chain}'  # Adjust column names to match your DataFrame
        colBF = df[col]  # Use the 'BF(dB)' columns as-is
        ranked = dict((locus, rank+1) for (rank, locus) in enumerate(colBF.sort_values(ascending=False).index))
        df[f'rank_BF_{chain}'] = [ranked[locus] for locus in df.index]  # Adding the rank column for each chain

    # mean BF across chains
    print('getting mean BF')
    bfcols = [f'BF(dB)_{chain}' for chain in chains]  # List of BF columns across chains
    df['mean_BF(dB)'] = df[bfcols].mean(axis=1)  # Mean of BF values across chains

    # rank of mean BF
    print('ranking mean BF')
    ranked = dict((locus, rank+1) for (rank, locus) in enumerate(df['mean_BF(dB)'].sort_values(ascending=False).index))
    df['rank_mean_BF(dB)'] = [ranked[locus] for locus in df.index]  # Rank based on the mean BF

    # bool column: BF >= 20 for >=3 chains, and BF >= 15 for >= 3 chains
    print('calculating bool column')
    df['BF(dB)_gte20_for-gte3chains'] = (df[bfcols] >= 1).sum(axis=1) >= 3
    df['BF(dB)_gte15_for-gte3chains'] = (df[bfcols] >= 1).sum(axis=1) >= 3

    return df

# Example of calling the function


# In[ ]:





# In[139]:


bio1_results[bio1_results['BF(dB)_chain3'] > 0 ]


# In[ ]:


bio1_results[bio1_results['BF(dB)_chain2'] > 0]


# In[ ]:


bio1_results[bio1_results['BF(dB)_chain1'] > 0]


# In[24]:


bio1_results[(bio1_results['BF(dB)_chain1'] > 0) & (bio1_results['BF(dB)_chain2'] > 0) & (bio1_results['BF(dB)_chain3'] > 0)]


# In[ ]:





# In[140]:


sns.histplot(bio1_results['BF(dB)_chain3'].sample(1000))


# In[ ]:





# In[ ]:


sns.histplot(bio1_results['BF(dB)_chain3'])


# In[121]:


bio1_results


# In[27]:


bio1_results_bf = rank_baypass(bio1_results)


# In[150]:


bio1_results_bf[bio1_results_bf['mean_BF(dB)'] > 0 ]


# In[ ]:





# In[28]:


bio1_results[bio1_results['BF(dB)_gte20_for-gte3chains'] == True]


# In[111]:


import math

# Calculate the number of rows corresponding to the top 1% of ranked SNPs
oneperc = math.ceil(0.01 * bio1_results.shape[0])  # bio1_results.shape[0] gives the number of rows

# Identify the columns that contain the ranking information for each chain
rankcols = [col for col in bio1_results.columns if 'rank_chain' in col]

# Create a new column for rank consistency: True if ranked in the top 1% for >= 3 chains
bio1_results['rank_consistency_top1perc_for-gte3chains'] = (bio1_results[rankcols] < oneperc).sum(axis=1) >= 3


# In[112]:


bio1_results


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[39]:


shfiles = []

file = f'baypass_postproc.sh'
text = f'''#!/bin/bash
#SBATCH --job-name=baypass_postproc
#SBATCH --time=48:00:00  
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=40gb
#SBATCH --cpus-per-task=1
#SBATCH --output=baypass_postproc.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

# Source the conda.sh script
source /home/tbellagio/miniforge3/etc/profile.d/conda.sh

# Activate the conda environment
conda activate /home/tbellagio/miniforge3/envs/pipeline_snakemake

cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/
python read_results.py 
    
    
    '''
with open(file, 'w') as o:
    o.write("%s" % text)
shfiles.append(file)


# In[40]:


subprocess.run(["sbatch", shfiles[0]], check=True)


# In[ ]:





# In[9]:


bio1_results = pd.read_csv('output/combined_covariable_bio1.csv')


# In[11]:


len(bio1_results)/3


# In[12]:


bio1_results


# In[16]:


pd.read_csv('output/old/combined_covariable_bio1.csv', nrows=10)


# In[ ]:





# In[46]:


ex = ex.drop_duplicates()


# In[47]:


ex


# In[69]:





# In[71]:


bf_file = pd.read_csv('results/' + bffiles[0], sep = '\s+')#.to_csv('chat.csv') #['COVARIABLE']


# In[72]:


gdict = dict((num+1, locus) for num,locus in enumerate(gloci))


# In[75]:


bf_file['MRK'] = bf_file['MRK'].map(gdict)


# In[76]:


bf_file


# In[46]:


partitions


# In[ ]:





# In[ ]:





# In[ ]:


def label_snps(resfile):
    """Label baypass resfiles with SNP ID using pkl files saved in 002_kickoff."""
    import os
    import pandas
    import tqdm
    # find the gfile associated with resfile
    prefix = os.path.basename(resfile).split("_chain")[0]
    chain = "chain_%s" % os.path.basename(resfile).split(prefix)[1].split("_")[2]
    
    # get the pkl file that has the ordered list of SNP IDs for the gfile
    if 'new' in resfile:
        newgfile_dir = os.path.dirname(resfile).replace("new_results", "new_individual_gfiles")
        pklfile = os.path.join(newgfile_dir, f'{prefix}_noheaderidx.pkl')
    else:
        pklfile = os.path.join(gfile_dir, f'{prefix}_noheaderidx.pkl')
    gloci = pklload(pklfile)  # read in rows with SNPIDs, but only one column
    
    # read in the resfile
    print('reading table ...')
    df = pandas.read_table(resfile, delim_whitespace=True)
    
    # map MRK column (SNP index) to locus name
    print('getting gdict ...')
    gdict = dict((num+1, locus) for num,locus in enumerate(gloci))
    
    # map covariable column to environment ID
    print('getting edict ...')
    edict = dict((num+1, env) for num,env in enumerate(efile.index))
    
    # split each covariable (env) in to its own file, label index with SNP ID
    df['env'] = df['COVARIABLE'].map(edict)
    envdfs = {chain:{}}
    for env in tqdm.tqdm(uni(df['env']), desc='annotating snps'):
        envdfs[chain][env] = df[df['env']==env].copy()
        envdfs[chain][env].index = envdfs[chain][env]['MRK'].map(gdict)
        envdfs[chain][env].index.names = ['']

    return envdfs

efile = pd.read_table(op.join(baydir, 'coastal_efile_std_HEADERIDX.txt'), index_col=0, usecols=[0,1])
dview['efile'] = efile
dview['uni'] = uni
dview['gfile_dir'] = gfile_dir
dview['pklload'] = pklload


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[21]:


import subprocess

def grep_search(directory, search_text):
    try:
        # Command to execute grep
        command = ['grep', '-rl', search_text, directory]
        # Running the command
        result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # Printing the result
        if result.stdout:
            print("Files containing the search term:")
            print(result.stdout)
        else:
            print("No files found with the search term.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Usage
grep_search('cmd_files/shfiles/', 'partition_1680_chain_3')


# In[22]:


# Open the file
with open('cmd_files/shfiles/' + 'batch_102.sh', 'r') as file:
    content = file.read()


# In[23]:


content


# In[25]:


pd.read_csv('individual_gfiles/' + 'partition_1680.txt',header=None)


# In[28]:





# In[30]:


len(data_loaded)


# In[31]:


data_loaded


# In[ ]:




