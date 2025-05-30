#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import os 
import allel
import random
import subprocess
import pickle
from os import path as op
import numpy as np
from dask import delayed
import dask.dataframe as dd
from dask import delayed, compute


# In[3]:


key_files = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/key_files/'


# In[ ]:





# In[3]:


pwd -P


# In[4]:


flowers = pd.read_csv(key_files + 'merged_sample_table.csv')

#path_meixi = '/carnegie/data/Shared/Labs/Moi/Everyone/meixilin'

num_flowers_map = flowers.set_index('sample_name')['total_flower_counts'].to_dict()

#path_all_af_indexed = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/merged_hapFIRE_allele_frequency_indexed.csv'

#Choose SNPs for GEA structure correction in baypass

#os.listdir('/home/tbellagio/safedata/meixilin/grenenet/metadata/data')

#pd.read_csv('/home/tbellagio/safedata/ath_evo/grenephase1-data/merged_tables/merged_hapFIRE_allele_frequency_LDpruned.txt', sep = '\t')

#all_samples_af = os.listdir('/carnegie/nobackup/scratch/xwu/grenet/hapFIRE_frequencies/samples/snp_frequency/')

## the size of this file is 3235481 * 745  
#path_all_af = '/carnegie/nobackup/scratch/xwu/grenet/merged_frequency/merged_hapFIRE_allele_frequency.csv'

## path to ld prunned allele freq
#path_ldp_af = '/carnegie/nobackup/scratch/xwu/grenet/merged_frequency/merged_hapFIRE_allele_frequency_LDpruned.txt'

first_gen = pd.read_csv(key_files + 'generation_1_sample_names.txt',header=None)[0]
#final_gen.append('0')


# In[ ]:





# In[5]:


#allele_freq_file = key_files + 'allele_freq_maf05_mincount05_firstgensamples.csv'


# In[ ]:





# In[ ]:





# In[6]:


allele_freq_ldp = key_files + 'merged_hapFIRE_allele_frequency_LDpruned.txt'


# In[7]:


# Read the file with only the specified columns in this case the final gen (355 columns)
df_ldp = dd.read_csv(allele_freq_ldp, sep = '\t', usecols = first_gen)


# In[8]:


df_shape = df_ldp.shape


# In[9]:


df_shape


# In[10]:


# The number of rows needs to be computed; the number of columns is immediate
num_rows = df_shape[0].compute()  # This computes the actual number of rows
num_columns = df_shape[1] 


# In[11]:


num_rows


# In[12]:


num_columns


# In[13]:


## get the vcf file frothe chromosomes and positions 

## ld pruned vcf file
ld_prunned_vcf_file = key_files + 'greneNet_final_v1.1_LDpruned.recode.vcf'
ld_prunned_vcf = allel.read_vcf(ld_prunned_vcf_file)


# In[14]:


ld_prunned_chrom = ld_prunned_vcf['variants/CHROM']
ld_prunned_pos = ld_prunned_vcf['variants/POS']


# In[15]:


## convert from dask to common pandas df
df_ldp = df_ldp.compute() 


# In[16]:


flowers = flowers[flowers['sample_name'].isin(first_gen)]


# In[17]:


flowers.head(2)


# In[18]:


## generate the allale counts
num_flowers_map = flowers.set_index('sample_name')['total_flower_counts'].to_dict()

allele_counts = {}
for i in df_ldp.columns:
    num_flowers = num_flowers_map[i]
    allele_counts[i + '_minor'] = df_ldp[i] * num_flowers * 2
    maj = 1 - df_ldp[i]
    allele_counts[i + '_major'] = maj * num_flowers * 2

allele_counts = pd.concat(allele_counts,axis=1)


# In[20]:


allele_counts = allele_counts.round().astype(int)


# In[21]:


#allele_counts.to_csv('allele_counts_minor_major_LDP_bayapass.csv',index=None)


# In[ ]:





# In[22]:


## generate the pool sizes 
#pool_sizes = {}
#for i,j in zip(allele_counts.columns[:], allele_counts.columns[1:]):
#    if i.strip('_minor') == j.strip('_major'):
#        pool_sizes[i.strip('_minor')] = (allele_counts[i][0] + allele_counts[j][0]) * 2 ## pool size = 2*#individuals

#pool_sizes = pd.concat(pool_sizes,axis=1).drop_duplicates()


# In[ ]:





# In[25]:


path_analysis_neutral = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/neutral_runs/'


# In[26]:


path_analysis_neutral + 'allele_counts_1stgen_wheader.txt'


# In[27]:


# save, no need to name as 'neutral'
allele_counts.to_csv(path_analysis_neutral + 'allele_counts_1stgen_wheader.txt', sep='\t', index=False)
allele_counts.to_csv(path_analysis_neutral + 'allele_counts_1stgen_nheader.txt', sep='\t', index=False, header=False)


# In[28]:


def makedir(directory: str) -> str:
    """If directory doesn't exist, create it.

    Return directory.
    """
#     if not op.exists(directory):
    os.makedirs(directory, exist_ok=True)
    
    return directory


# In[29]:


### it took 11 hs to run for ~13 k snps and 60 samples, 
## now im running 13 k snps and almost 300 samples so i will give it more time 


# In[ ]:


-nthreads 8


# In[30]:


# create sbatch files to run baypass omega estimation
# 
gfile = path_analysis_neutral + 'allele_counts_1stgen_nheader.txt'

## create a dir 
resultsdir = makedir(path_analysis_neutral + 'results/')
shdir = makedir(path_analysis_neutral + 'shfiles/')

# create sbatch files to submit on cedar server
shfiles = []
for i in range(5):
    seed = random.randint(1,100000000)
    file = shdir + 'chain_%s.sh' % str(i+1)
    cmd = f'/home/tbellagio/bin/baypass -gfile {gfile} -seed {seed} -nthreads 4 -print_omega_samples -outprefix chain_{i+1}'
    print(cmd,'\n')
    text = f'''#!/bin/bash
#SBATCH --job-name=chain_{i+1}
#SBATCH --time=2-00:00:00  # Time limit set to 4 days
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --cpus-per-task=4
#SBATCH --output=chain_{i+1}_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

module load python/3.11_conda
conda activate /home/tbellagio/miniforge3/envs/run_baypass
export LD_LIBRARY_PATH="/home/tbellagio/miniforge3/envs/run_baypass/lib:$LD_LIBRARY_PATH"
cd /home/tbellagio/scratch/gea_grene-net/baypass_terminal/neutral_runs
{cmd}


'''
    with open(file, 'w') as o:
        o.write("%s" % text)
    shfiles.append(file)


# In[31]:


## now run the shfiles
for shfile in shfiles:
    # Submit each sbatch script to the SLURM scheduler
    subprocess.run(["sbatch", shfile], check=True)


# In[ ]:





# In[ ]:


rsync -av --include '*/' --exclude 'results/' /home/tbellagio/scratch/gea_grene-net/baypass_terminal tbellg@hpc.brc.berkeley.edu:/global/scratch/users/tbellg/gea_grene-net/


# In[ ]:


### create files for savio 


# In[32]:


# create sbatch files to run baypass omega estimation
path_analysis = '/global/scratch/users/tbellg/gea_grene-net/baypass_terminal/'
path_analysis_calc = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/'
# 
gfile = path_analysis + 'allele_counts_3rdgen_nheader.txt'
#poolsizefile = path_analysis + 'pool_sizes_3rdgen_nheader.txt'


# In[23]:


# create sbatch files to run baypass omega estimation
path_analysis_calc = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/'
# 
gfile = path_analysis + 'allele_counts_3rdgen_nheader.txt'
poolsizefile = path_analysis + 'pool_sizes_3rdgen_nheader.txt'

## create a dir 
resultsdir = makedir(path_analysis_calc + 'neutral_runs/results/')
shdir = makedir(path_analysis_calc + 'neutral_runs/shfiles_caltech/')

# create sbatch files to submit on cedar server
shfiles = []
for i in range(5):
    seed = random.randint(1,100000000)
    file = shdir + 'chain_%s.sh' % str(i+1)
    cmd = f'/home/tbellagi/bin/baypass -gfile {gfile} -poolsizefile {poolsizefile} \
-nthreads 8 -seed {seed} -print_omega_samples -outprefix chain_{i+1}'
    print(cmd,'\n')
    text = f'''#!/bin/bash
#SBATCH --job-name=chain_{i+1}
#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=40gb
#SBATCH --cpus-per-task=10
#SBATCH --output=chain_{i+1}_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

cd /home/tbellagi/bigscratch/gea/baypass_terminal
{cmd}


'''
    with open(file, 'w') as o:
        o.write("%s" % text)
    shfiles.append(file)


# In[ ]:





# In[20]:


## now run the shfiles
for shfile in shfiles:
    # Submit each sbatch script to the SLURM scheduler
    subprocess.run(["sbatch", shfile], check=True)


# In[ ]:





# In[122]:


### create env file 


# In[9]:


path_analysis = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/'


# In[10]:


#finalgen_samples = pd.read_csv(key_files + 'final_gen.csv')['sample_name']


# In[11]:


#finalgen_samples = pd.read_csv('../final_gen.csv')['sample_name']
first_gen_samples = pd.read_csv('../key_files/generation_1_sample_names.txt',header=None)[0]

samples = first_gen_samples.to_list()

#clim_sites_during_exp = pd.read_csv('/carnegie/nobackup/scratch/tbellagio/grene/data/bioclimvars_experimental_sites_era5.csv')
clim_sites_during_exp = pd.read_csv('../key_files/bioclimvars_sites_era5_year_2018.csv')

sites_af = pd.Series(samples).str.split('_').str[0].astype(int)

sites_af.name = 'site'

env = sites_af.reset_index().merge(clim_sites_during_exp).drop(['index'],axis=1)


# In[12]:


sites = env[['site', 'bio1', 'bio17']]


# In[ ]:





# In[13]:


#pca_clim_sites_during_exp = pd.read_csv('/carnegie/nobackup/scratch/tbellagio/gea_grene-net/pca_climate_var_experimental_sites.csv').drop('Unnamed: 0',axis=1)


# In[14]:


#clim = pd.merge(clim_sites_during_exp, pca_clim_sites_during_exp)


# In[15]:


## for now im gonna run the analysis on bio1, bio12 and  and the 2 first components of the pca 


# In[16]:


sites = sites.set_index('site')


# In[17]:


## save the environmental var names 
import pickle

# Some data to pickle
env_vars_names = sites.columns.tolist()
# Pickle the data into a file
with open(path_analysis + 'env_vars_names', 'wb') as f:
    pickle.dump(env_vars_names, f)


# In[18]:


from sklearn.preprocessing import StandardScaler


# In[19]:


sites


# In[21]:


# Standardize the data
scaler = StandardScaler()
sites_scaled = scaler.fit_transform(sites)


# In[22]:


env = pd.DataFrame(data=sites_scaled, index = sites.index).T


# In[23]:


env


# In[46]:


# save
# save, no need to name as 'neutral'
env.to_csv(path_analysis + 'env_firstgen_wheader.txt', sep='\t', index=False)
env.to_csv(path_analysis + 'env_firstgen_nheader.txt', sep='\t', index=False, header=False)


# In[ ]:





# In[ ]:





# In[ ]:





# In[145]:


### now create the gfile fro the whole run 


# In[30]:


path_analysis = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/'


# In[31]:


gfiles_path = os.path.join(path_analysis, 'individual_gfiles')
os.makedirs(gfiles_path, exist_ok=True)


# In[35]:


def generate_allele_counts(partition, num_flowers_map, first_gen):
    ## adjust the order so that all files have the same ########################
    partition = partition[first_gen]
    ## adjust the order ########################
    
    ## check the order 
    #if (partition.columns == first_gen).all():
    #    print('correct order 1')
    #else:
    #    print('wrong order 1')
    ## process all the columns but the index
    index_to_save = partition['0']
    partition = partition.drop('0',axis=1)
    
    allele_counts = {}
    #print(partition.columns)
    for column in partition.columns:
        num_flowers = num_flowers_map.get(column, 0)  # Default to 0 if not found
        allele_counts[column + '_minor'] = partition[column] * num_flowers * 2
        allele_counts[column + '_major'] = (1 - partition[column]) * num_flowers * 2
    # Explicitly convert to DataFrame
    result_df = pd.DataFrame(allele_counts).round().astype(int)
    result_df['positions'] = index_to_save
    
    ## check the order 

    #if (result_df.columns == correct_order).all():
    #    print('correct order')
    #else:
    #    print('wrong order')
    return result_df


# In[ ]:





# In[36]:


def save_partition(partition, partition_index, gfiles_path):
    print(len(partition))
    filename = f"{gfiles_path}/partition_{partition_index}.txt"
    ## save the index and the data separetly 
    index_to_save = partition['positions']
    partition = partition.drop('positions',axis=1)
    # Save DataFrame with index
    partition.to_csv(filename, index=False, header=False)  # Include the index when saving
    
    ## just in case save also de columns names 
    column_names = partition.columns  # Extract the index list
    with open(gfiles_path + f'/column_names_partition_{partition_index}', 'wb') as f:
        pickle.dump(column_names, f)
    
    loci_list = index_to_save.tolist()  # Extract the index list
    with open(gfiles_path + f'/loci_partition_{partition_index}', 'wb') as f:
        pickle.dump(loci_list, f)
        
    return filename


# In[ ]:





# In[12]:


first_gen = pd.read_csv(key_files + 'generation_1_sample_names.txt',header=None)[0]
allele_freq_file = key_files + 'allele_freq_maf05_mincount05_firstgensamples.csv'


# In[ ]:





# In[ ]:


#af = dd.read_csv(key_files + 'merged_hapFIRE_allele_frequency.txt', sep = '\t', usecols = first_gen)

#af = af.compute()

#af = af.reset_index(drop=True)

#snp_dict = pd.read_csv(key_files + 'var_pos_grenenet.csv')

#mask = snp_dict['total_alleles05filter'].notna()

#af = af[mask]

#af.to_csv(key_files + 'allele_freq_maf05_mincount05_firstgensamples.csv', index=None)


# In[ ]:





# In[146]:


1242385/174


# In[36]:


gfile_ddf = dd.read_csv(allele_freq_file)
af = gfile_ddf.compute()


# In[13]:


snps_dict = pd.read_csv(key_files +  'var_pos_grenenet.csv')


# In[14]:


snps_dict = snps_dict[snps_dict['total_alleles05filter'].notna()].reset_index(drop=True)


# In[15]:


af = af.reset_index(drop=True)


# In[16]:


af


# In[17]:


af[0] = snps_dict['id']


# In[18]:


af.to_csv(key_files + 'allele_freq_maf05_mincount05_firstgensamples_windexsnpid.csv',index=None)


# In[30]:


af


# In[19]:


allele_freq_file_w_index = key_files + 'allele_freq_maf05_mincount05_firstgensamples_windexsnpid.csv'


# In[31]:


# Load the data
## make 40MB partitions so that aprox 3500 snps
gfile_ddf = dd.read_csv(allele_freq_file_w_index, blocksize='25MB')
# Now repartition if necessary
#gfile_ddf = gfile_ddf.repartition(npartitions=900)


# In[32]:


first_gen = pd.read_csv(key_files + 'generation_1_sample_names.txt',header=None)[0]
first_gen = first_gen.tolist()
first_gen.append('0')


# In[ ]:





# In[37]:


# dont need it for now 
processed_partitions = gfile_ddf.map_partitions(generate_allele_counts, num_flowers_map, first_gen)

save_tasks = [delayed(save_partition)(processed_partitions.get_partition(i), i, gfiles_path)
          for i in range(processed_partitions.npartitions)]


# In[38]:


#save_tasks = [delayed(save_partition)(gfile_ddf.get_partition(i), i, gfiles_path)
#          for i in range(gfile_ddf.npartitions)]


# In[39]:


len(save_tasks)


# In[120]:


140 * 5259


# In[117]:


## used for checking 
#save_tasks = [delayed(save_partition)(processed_partitions.get_partition(i), i, gfiles_path)
#          for i in range(5)]


# In[ ]:





# In[40]:


# Compute all tasks
results = compute(*save_tasks)


# In[165]:


## after this there should be 1700 * 3 files in the folder including gfiles, snpspositions and column names 


# In[ ]:





# In[ ]:





# In[32]:


data_all = []
for i in range(gfile_ddf.npartitions):
    pickle_file_path = gfiles_path + f'/loci_partition_{i}'
    with open(pickle_file_path, 'rb') as file:
        data = pickle.load(file)
        data_all.append(data)


# In[25]:


from itertools import chain

# Flatten the list of lists
flat_list = list(chain.from_iterable(data_all))

# Convert the flattened list to a set and compare lengths to check for duplicates
unique_data = set(flat_list)
if len(unique_data) != len(flat_list):
    print("There are repeated values")
else:
    print("All values are unique")


# In[ ]:





# In[26]:


gfile_dir = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/individual_gfiles/'

path_analysis = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/'

# get gfiles
gfiles = [gfile_dir + i for i in os.listdir(gfile_dir) if '.txt' in i]


# In[27]:


partition_214 = pd.read_csv(gfiles[0],header=None)


# In[28]:


part = gfiles[0].split('/')[-1]


# In[33]:


pickle_file_path = gfiles_path + f'/column_names_' + part.replace('.txt','')
with open(pickle_file_path, 'rb') as file:
    data214_columns = pickle.load(file)


# In[34]:


pickle_file_path = gfiles_path + f'/loci_partition_0'
with open(pickle_file_path, 'rb') as file:
    data214 = pickle.load(file)


# In[35]:


partition_214.columns = data214_columns


# In[36]:


minor_cols = [i for i in partition_214.columns if 'minor' in i]


# In[37]:


partition_214 = partition_214[minor_cols]


# In[ ]:





# In[38]:


pos_column = pd.read_csv(path_analysis + 'index_to_add.csv')


# In[39]:


## making sure all the positions are present
are_equal = set(pos_column['0']) == unique_data
print(are_equal)


# In[40]:


partition_214 = pd.read_csv(path_analysis + 'individual_gfiles/partition_0.txt',header=None)


# In[41]:


#os.listdir(path_analysis + 'individual_gfiles/')


# In[42]:


pickle_file_path = gfiles_path + f'/loci_partition_0'
with open(pickle_file_path, 'rb') as file:
    data214 = pickle.load(file)


# In[43]:


pickle_file_path = gfiles_path + f'/column_names_partition_0'
with open(pickle_file_path, 'rb') as file:
    data214_columns = pickle.load(file)


# In[44]:


partition_214.columns = data214_columns


# In[45]:


partition_214.index = data214


# In[46]:


partition_214


# In[47]:


minor_cols = [i for i in partition_214.columns if 'minor' in i]


# In[48]:


partition_214 = partition_214[minor_cols]


# In[ ]:





# In[ ]:





# In[49]:


partition_214.sum(axis=1)[partition_214.sum(axis=1) == 0]


# In[50]:


rs = (partition_214 != 0).sum(axis=1)


# In[51]:


rs[rs == 0 ]


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[52]:


partition_214['pos'] = data214


# In[53]:


check = pd.read_csv(allele_counts, nrows = 5493)


# In[101]:


#check = check[final_gen]


# In[102]:


#check.columns = columns_for.columns 


# In[103]:


#check = generate_allele_counts(check, num_flowers_map, final_gen)


# In[96]:


check = check.set_index('positions')


# In[90]:


partition_214 = partition_214.set_index('pos')


# In[105]:


check


# In[106]:


partition_214


# In[62]:


(partition_214 == check).all().all()


# In[ ]:





# In[ ]:





# In[ ]:





# In[125]:


from typing import Optional, Union


# In[126]:


def read(file: str, lines=True, ignore_blank=False) -> Union[str, list]:
    """Read lines from a file.

    Return a list of lines, or one large string
    """
    with open(file, "r") as o:
        text = o.read()

    if lines is True:
        text = text.split("\n")
        if ignore_blank is True:
            text = [line for line in text if line != ""]

    return text


# In[ ]:





# In[ ]:





# In[18]:


## check omegas for convergence 


# In[5]:


path_analysis = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/neutral_runs/'


# In[6]:


files = os.listdir(path_analysis)


# In[7]:


result_files = [i for i in files if 'mat_omega' in i]


# In[8]:


result_files


# In[9]:


filename = 'chain_3_mat_omega.out'
with open(path_analysis + filename, 'r') as file:
    content = file.read()


# In[10]:


def read(filename):
    """Read the file and return a list of its lines."""
    with open(filename, 'r') as file:
        return file.readlines()


# In[11]:


# get a list of values by vectorizing the matrix for each chain
matrices = {}
for m in result_files:
    chain = op.basename(m).split("_mat")[0]
    text = read(path_analysis + m)
    rows = []
    for line in text:
        rows.extend([float(val) for val in line.split()])
    matrices[chain] = rows
len(matrices[chain])


# In[12]:


np.sqrt(len(matrices[chain]))  # num pops


# In[13]:


from scipy.stats import pearsonr


# In[14]:


# look at pairwise correlations among chains
for i,chaini in enumerate(matrices):
    for j,chainj in enumerate(matrices):
        if i < j:
            print(chaini,chainj,pearsonr(matrices[chaini], matrices[chainj]))


# In[15]:


# plot pairwise chains agains 1:1 line, print slope and intercept from data
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms

def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')

for i,chaini in enumerate(matrices):
    for j,chainj in enumerate(matrices):
        if i < j:
            slope, intercept = np.polyfit(matrices[chaini], matrices[chainj], 1)
            print(chaini, chainj, slope, intercept)
            fig, ax = plt.subplots(figsize=(5,5))
            ax.scatter(matrices[chaini], matrices[chainj])
            ax.set_xlabel(chaini)
            ax.set_ylabel(chainj)
            line = mlines.Line2D([0, 1], [0, 1], color='red')
            transform = ax.transAxes
            line.set_transform(transform)
            ax.add_line(line)           
            
            plt.show()


# In[16]:


result_files.sort()


# In[17]:


result_files


# In[18]:


# make a matrix with the average across matrices
mats = {}
for m in result_files:
    chain = op.basename(m).split("_mat")[0]
    mats[chain] = pd.read_table(path_analysis + m, delim_whitespace=True, header=None)
mats[chain].head()


# In[19]:


mats.keys()


# In[20]:


# sum across each matrix (to use in average calc below)
for chain,m in mats.items():
    print(chain)
    if chain == 'chain_1':
        summatrix = np.array(m)
    else:
        summatrix = summatrix + np.array(m)
summatrix


# In[21]:


# calculate average netural matrix
avg = summatrix / len(mats)
avg.shape


# In[22]:


# convert to dataframe
avgmat = pd.DataFrame(avg)
print(avgmat.shape)
avgmat.head()


# In[23]:


# make sure avgmat is highly correlated with other chains
avgmatlst = []
for row in avgmat.index:
    for col in avgmat.columns:
        avgmatlst.append(avgmat.loc[row, col])
print(len(avgmatlst), len(matrices[chain]))

for chain,lst in matrices.items():
    print(chain, pearsonr(matrices[chain], avgmatlst))


# In[24]:


path_analysis


# In[25]:


# save
# save, no need to name as 'neutral'
avgmat.to_csv(path_analysis + 'omegaavg_nheader_first_gen.txt', sep='\t', index=False, header=False)


# In[26]:


path_analysis = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/'
avgmat.to_csv(path_analysis + 'omegaavg_nheader_first_gen.txt', sep='\t', index=False, header=False)


# In[27]:


## check order populations and assign it to the omega matrix 


# In[ ]:





# In[29]:


from_ac = pd.Series(pd.read_csv(path_analysis + 'neutral_runs/allele_counts_1stgen_wheader.txt', sep='\t').columns.str.replace('_minor', '').str.replace('_major', '')).drop_duplicates().reset_index(drop=True)


# In[30]:


from_ac


# In[31]:


#(from_ps == from_ac).all()


# In[32]:


avgmat.columns = from_ac
avgmat.index = from_ac


# In[33]:


avgmat.to_csv(path_analysis + 'omegaavg_wheader_first_gen.txt', sep='\t', index=True, header = True)


# In[ ]:





# In[ ]:





# In[30]:


### prepare everything to run baypass


# In[5]:


gfile_dir = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/individual_gfiles'

path_analysis = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/'

# get gfiles
gfiles = [i for i in os.listdir(gfile_dir) if '.txt' in i]
len(gfiles)


# In[ ]:





# In[6]:


env = pd.read_csv(path_analysis + 'env_firstgen_nheader.txt', sep='\t')


# In[7]:


env


# In[44]:


#pool_sizes = pd.read_csv(path_analysis + 'pool_sizes_3rdgen_nheader.txt', sep='\t')


# In[8]:


get_ipython().system('pwd')


# In[46]:


gfile = gfiles[0]


# In[47]:


gfile


# In[48]:


# tbellg@hpc.brc.berkeley.edu:/global/scratch/users/tbellg/gea_grene-net/


# In[49]:


path_analysis_savio = '/global/scratch/users/tbellg/gea_grene-net/baypass_terminal/'


# In[50]:


path_analysis_caltech = '/central/scratch/tbellagi/gea/baypass_terminal/'


# In[51]:


# -poolsizefile {path_analysis_savio}pool_sizes_3rdgen_nheader.txt \
# -auxmodel \
# -nthreads 8 \


# In[52]:


# create baypass commands for each gfile for 5 chains each
cmds = []
for gfile in gfiles:
    for chain in ['chain_1', 'chain_2', 'chain_3']:
        seed = random.randint(1, 100000)
        bname = gfile.replace('.txt', '')
        outprefix = f"{bname}_{chain}"
        cmd = f'baypass \
-gfile {path_analysis_caltech}individual_gfiles/{gfile} \
-efile {path_analysis_caltech}env_firstgen_nheader.txt \
-omegafile {path_analysis_caltech}omegaavg_nheader_first_gen.txt \
-outprefix {outprefix} \
-seed {seed} \
-pilotlength 1000 \
-nval 50000 \
-nthreads 8 \
-npilot 20 '
        cmds.append(cmd)
len(cmds)


# In[53]:


678 / 3


# In[54]:


cmds[0].split(' -')


# In[55]:


path_analysis


# In[56]:


cmddir = path_analysis + 'cmd_files/'


# In[57]:


# create catfiles for baypass commands - for requests of a single node with 48 CPUs
rundir = cmddir + 'run_01/'
os.makedirs(rundir,exist_ok=True)
catfiles = []
tocat = []
for i,cmd in enumerate(cmds):
    tocat.append(cmd)
    if len(tocat) == 8 or (i+1) == len(cmds):
        fill = str(len(catfiles)).zfill(4)
        file_path = rundir + f'catfile_{fill}.txt'
        with open(file_path, 'w') as o:
            o.write('\n'.join(tocat))
        tocat = []
        catfiles.append(file_path)
len(catfiles)


# In[58]:


resdir = os.path.join(path_analysis, 'cmd_files/results/')
shdir = os.path.join(path_analysis, 'cmd_files/shfiles/')

# Create the directories if they don't exist
os.makedirs(resdir, exist_ok=True)
os.makedirs(shdir, exist_ok=True)


# In[59]:


def read_node_list(filename):
    with open(filename, 'r') as file:
        nodes = [line.strip() for line in file if line.strip()]
    return nodes


# In[60]:


## replace catfiles path to the caltech one
catfiles = pd.Series(catfiles).str.replace('/carnegie/nobackup/scratch/tbellagio/gea_grene-net', '/central/scratch/tbellagi/gea')


# In[61]:


catfiles[0]


# In[ ]:





# In[ ]:


cat /central/scratch/tbellagi/gea/baypass_terminal/cmd_files/run_01/catfile_0000.txt | parallel -j 64
>parallel_output_%j.txt 2>parallel_error_%j.txt


# In[ ]:


parallel -j 64 >parallel_output_%j.txt 2>parallel_error_%j.txt


# In[ ]:


parallel -j 64 "awk '{...}' > output_{#}.txt"


# In[ ]:


seq 10 | parallel awk \''{...}'\' file{}.txt ">" file{}.out


# In[ ]:


seq 10 | parallel " awk command > file{}.out "


# In[ ]:


module load parallel/20220522-gcc-11.3.1-ki3yyer


# In[ ]:


cat red_catfile.txt | parallel -j 2 >parallel_output_%j.txt 2>parallel_error_%j.txt


# In[ ]:


cat red_catfile.txt | parallel -j 2 " awk command > file{#}.out "


# In[ ]:


cat red_catfile.txt | parallel -j 2 awk \''{...}'\' file_8t{}.txt ">" file_8t{}.out


# In[ ]:


cat red_catfile.txt | parallel -j 2 "{} > output_8t_{#}.out 2> output_8t_{#}.err"


# In[ ]:


cat red_catfile.txt | parallel -j 2 "{} > output_8t_{%}.out 2> output_8t_{%}.err"


# In[185]:


catfiles[0].split("_")[-1].split(".")[0]


# In[ ]:


parallel -j 8 "{} > output_8t_{%}.out 2> output_8t_{%}.err"


# In[ ]:


cat red_catfile2.txt | parallel -j 2 "{} > output_8t_{%}.out 2> output_8t_{%}.err"


# In[62]:


######### FOR CALTECH ###################
######### FOR CALTECH ###################
######### FOR CALTECH ###################
## #SBATCH --cores-per-socket=32 this will make sure i only get the 64 cores nodes
node_names_64 = read_node_list('not_node_names_64.txt')
nodelist = ','.join(node_names_64)

shfiles = []
for catfile in catfiles:
    num = catfile.split("_")[-1].split(".")[0]
    shfile = op.join(shdir, op.basename(f"batch_{num}.sh"))
    with open(shfile, 'w') as o:
        text = f'''#!/bin/bash
#SBATCH --job-name=baypass_cmd_{num}
#SBATCH --partition=expansion
#SBATCH --sockets-per-node=2
#SBATCH --cores-per-socket=32
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=25
#SBATCH --cpus-per-task=8
#SBATCH --output=batch_{num}_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

cd /central/scratch/tbellagi/gea/baypass_terminal/results
module load parallel/20220522-gcc-11.3.1-ki3yyer

cat {catfile} | parallel -j 8 "{{}} > output_catfile{num}_{{%}}.out 2> output_catfile{num}_{{%}}.err"

'''
        o.write("%s" % text)
    shfiles.append(shfile)
len(shfiles)

######### FOR CALTECH ###################
######### FOR CALTECH ###################
######### FOR CALTECH ###################


# In[ ]:





# In[ ]:


rsync -av --dry-run \
--include 'individual_gfiles/***' \
--include 'cmd_files/***' \
--include 'env_firstgen_nheader.txt' \
--include 'omegaavg_nheader_first_gen.txt' \
--include '*/' \
--exclude '*' \
--prune-empty-dirs \
/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/ \
tbellagi@login.hpc.caltech.edu:/central/scratch/tbellagi/gea/baypass_terminal/


# In[ ]:


rsync -av \
--include 'individual_gfiles/***' \
--include 'cmd_files/***' \
--include 'env_firstgen_nheader.txt' \
--include 'omegaavg_nheader_first_gen.txt' \
--include '*/' \
--exclude '*' \
--prune-empty-dirs \
tbellagio@calc.carnegiescience.edu:/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/ \
/Users/tatiana/Documents/grenenet/gea/baypass_terminal/


# In[ ]:


rsync -av \
--include 'individual_gfiles/***' \
--include 'cmd_files/***' \
--include 'env_firstgen_nheader.txt' \
--include 'omegaavg_nheader_first_gen.txt' \
--include '*/' \
--exclude '*' \
--prune-empty-dirs \
/Users/tatiana/Documents/grenenet/gea/baypass_terminal/ \
tbellagi@login.hpc.caltech.edu:/central/scratch/tbellagi/gea/baypass_terminal/


# In[ ]:





# In[ ]:





# In[ ]:


['baypass',
 'gfile /central/scratch/tbellagi/gea/baypass_terminal/individual_gfiles/partition_86.txt',
 'efile /central/scratch/tbellagi/gea/baypass_terminal/env_3rdgen_nheader.txt',
 'omegafile /central/scratch/tbellagi/gea/baypass_terminal/omegaavg_nheader.txt',
 'outprefix partition_86_chain_1',
 'seed 20642',
 'pilotlength 1000',
 'nval 50000',
 'nthreads 8',
 'npilot 20 ']


# In[ ]:


cat {catfile} | parallel -j 8 "{{}} > {{}} > output_catfile_{num}_{{#}}.out 2> output_catfile_{num}_{{#}}.err"


# In[ ]:


carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/cmd_files/run_01/catfile_0000.txt()


# In[95]:


catfile


# In[ ]:


# create slurm jobs, cat catfiles at GNU parallel
resdir = os.path.join(path_analysis, 'cmd_files/results/')
shdir = os.path.join(path_analysis, 'cmd_files/shfiles/')

# Create the directories if they don't exist
os.makedirs(resdir, exist_ok=True)
os.makedirs(shdir, exist_ok=True)

shfiles = []
for catfile in nb(catfiles):
    num = catfile.split("_")[-1].split(".")[0]
    shfile = op.join(shdir, op.basename(f"batch_{num}.sh"))
    with open(shfile, 'w') as o:
        text = f'''#!/bin/bash
#SBATCH --job-name=baypass_cmd_{num}
#SBATCH --account=fc_moilab
#SBATCH --partition=savio2_knl
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=25
#SBATCH --cpus-per-task=64
#SBATCH --output=batch_{num}_%j.out
#SBATCH --mail-user=tbellagio@berkeley.edu
#SBATCH --mail-type=FAIL

export PATH="${{PATH}}:/project/def-saitken/programs/baypass_2.2/sources"

# Initialize Mamba
source /global/home/users/tbellg/miniforge3/etc/profile.d/conda.sh

conda activate /global/home/users/tbellg/miniforge3/envs/baypass
cd /global/scratch/users/tbellg/gea_grene-net/baypass_terminal/results

# Add the conda environment's library path to LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/global/home/users/tbellg/miniforge3/envs/baypass/lib:$LD_LIBRARY_PATH

cat {catfile} | parallel -j 64 --progress

'''
        o.write("%s" % text)
    shfiles.append(shfile)
len(shfiles)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


# create slurm jobs, cat catfiles at GNU parallel
resdir = os.path.join(path_analysis, 'cmd_files/results/')
shdir = os.path.join(path_analysis, 'cmd_files/shfiles/')

# Create the directories if they don't exist
os.makedirs(resdir, exist_ok=True)
os.makedirs(shdir, exist_ok=True)

shfiles = []
for catfile in nb(catfiles):
    num = catfile.split("_")[-1].split(".")[0]
    shfile = op.join(shdir, op.basename(f"batch_{num}.sh"))
    with open(shfile, 'w') as o:
        text = f'''#!/bin/bash
#SBATCH --job-name=baypass_cmd_{num}
#SBATCH --partition=?
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=25
#SBATCH --cpus-per-task=?
#SBATCH --output=batch_{num}_%j.out
#SBATCH --mail-user=tbellagio@berkeley.edu
#SBATCH --mail-type=FAIL

conda activate /global/home/users/tbellg/miniforge3/envs/baypass
cd /home/tbellagi/bigscratch/gea/baypass_terminal

cat {catfile} | parallel -j 64 --progress

'''
        o.write("%s" % text)
    shfiles.append(shfile)
len(shfiles)


# In[ ]:





# In[ ]:





# In[ ]:





# In[222]:


path_analysis + 'cmd_files/shfiles'


# In[223]:


path_analysis + 'cmd_files/results'


# In[224]:


#SBATCH --qos=savio_long


# In[225]:


# create slurm jobs
# Define the directories based on path_analysis
resdir = os.path.join(path_analysis, 'cmd_files/results/')
shdir = os.path.join(path_analysis, 'cmd_files/shfiles/')

# Create the directories if they don't exist
os.makedirs(resdir, exist_ok=True)
os.makedirs(shdir, exist_ok=True)

shfiles = []
for num, cmd in enumerate(cmds):
    shfile = shdir + f"batch_{num}.sh"
    with open(shfile, 'w') as o:
        text = f'''#!/bin/bash

#SBATCH --job-name=baypass_cmd_{num}
#SBATCH --account=fc_moilab
#SBATCH --partition=savio2_knl
#SBATCH --time=3-00:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cores=1
#SBATCH --mem-per-cpu=25M
#SBATCH --cpus-per-task=8
#SBATCH --output=batch_{num}_%j.out
#SBATCH --mail-user=tbellagio@berkeley.edu
#SBATCH --mail-type=FAIL

# Initialize Mamba
source /global/home/users/tbellg/miniforge3/etc/profile.d/conda.sh

conda activate /global/home/users/tbellg/miniforge3/envs/baypass
cd /global/scratch/users/tbellg/gea_grene-net/baypass_terminal/results

# Add the conda environment's library path to LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/global/home/users/tbellg/miniforge3/envs/baypass/lib:$LD_LIBRARY_PATH

{cmd}
'''
        o.write("%s" % text)
        os.chmod(shfile, 0o755)  # Sets the file to be executable by user, readable and executable by group and others

    shfiles.append(shfile)

len(shfiles)


# In[226]:


shfiles[0]


# In[ ]:





# In[ ]:





# In[ ]:





# In[96]:


new = pd.read_csv('/home/tbellagio/scratch/gea_grene-net/baypass_terminal/individual_gfiles/partition_249.txt',header=None)


# In[97]:


new


# In[ ]:





# In[93]:


new.to_csv('/home/tbellagio/scratch/gea_grene-net/baypass_terminal/individual_gfiles/partition_248.txt', header=None,
          index=None)


# In[ ]:





# In[89]:


11960 / 2


# In[ ]:


# create sbatch files to run baypass omega estimation
path_analysis = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/'
# 
gfile = path_analysis + 'allele_counts_3rdgen_nheader.txt'
poolsizefile = path_analysis + 'pool_sizes_3rdgen_nheader.txt'

## create a dir 
shdir = makedir(path_analysis + 'neutral_runs/shfiles')

# create sbatch files to submit on cedar server
shfiles = []
for i in range(5):
    seed = random.randint(1,100000000)
    file = shdir + 'chain_%s.sh' % str(i+1)
    cmd = f'baypass -gfile {gfile} -poolsizefile {poolsizefile} \
-nthreads 8 -seed {seed} -print_omega_samples -outprefix chain_{i+1}'
    print(cmd,'\n')
    text = f'''#!/bin/bash
#SBATCH --job-name=chain_{i}
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=25
#SBATCH --cpus-per-task=8
#SBATCH --output=chain_{i}_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

export PATH="${{PATH}}:/home/tbellagio/bin"

cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/neutral_runs

{cmd}


'''
    with open(file, 'w') as o:
        o.write("%s" % text)
    shfiles.append(file)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[214]:


og = pd.read_csv(path_ldp_af, sep = '\t', usecols=final_gen)


# In[231]:


og = generate_allele_counts(og, num_flowers_map)


# In[ ]:





# In[232]:


og.loc[data,:]


# In[229]:


pickle_file_path = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/individual_gfiles/loci_partition_1.txt'
with open(pickle_file_path, 'rb') as file:
    data = pickle.load(file)


# In[233]:


pd.read_csv('/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/individual_gfiles/partition_1.txt')


# In[212]:


pd.read_csv('/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/individual_gfiles/partition_1.txt')


# In[ ]:





# In[197]:


import dask.dataframe as dd
from dask import delayed, compute

from dask.distributed import Client


# In[172]:


from dask.distributed import Client

client = Client()  # Starts a local Dask client

# Apply function across partitions
_ = df.map_partitions(calculate_and_save_allele_counts, 
                          meta=object,  # You might need to adjust this based on your function's return type
                          num_flowers_map=num_flowers_map).compute()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


callset_from_vcf


# In[ ]:


snps = parallel_read(f,
                     lview=lview, dview=dview,
                     verbose=False)


# In[ ]:





# In[29]:


# 13985 * 745 
af = pd.read_csv('/home/tbellagio/safedata/ath_evo/grenephase1-data/merged_tables/merged_hapFIRE_allele_frequency_LDpruned.txt', sep = '\t')


# In[40]:


af


# In[39]:


len(af.columns)


# In[21]:


af.columns.str.split('_').str[0].unique()


# In[27]:


af.columns.str.split('_').str[1]


# In[24]:


af.columns.str.split('_').str[2].unique()


# In[ ]:


safedata/ath_evo/grenephase1-data/merged_tables/merged_hapFIRE_allele_frequency.csv


# In[ ]:





# In[9]:


1503 * 1500


# In[ ]:




