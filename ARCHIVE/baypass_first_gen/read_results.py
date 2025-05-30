import pandas as pd
import os
import pickle

base_directory = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/results/'

##the bfffiles have the bayes factor for each of the covariates for each fo the markers 
bffiles = [i for i in os.listdir(base_directory) if 'betai_reg.out' in i]
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
        output_path = os.path.join(output_dir, f"combined_covariable_{covariable}.csv")
        df_group.drop(columns='COVARIABLE').to_csv(output_path, mode='a', index=False, header=not os.path.exists(output_path))


partitions = extract_partitions(bffiles)  # Assuming bffiles is defined somewhere

all_part_collected = list(partitions.keys())

output_directory = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/output/'
os.makedirs(output_directory, exist_ok=True)

for partition, chains in partitions.items():
    for chain in chains:
        process_file(partition, chain, output_directory)