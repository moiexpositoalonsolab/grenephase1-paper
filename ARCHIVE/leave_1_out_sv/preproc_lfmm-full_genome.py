import sys
import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import dask.dataframe as dd
import sys

#path_ldp_af = '/carnegie/nobackup/scratch/xwu/grenet/merged_frequency/merged_hapFIRE_allele_frequency_LDpruned.txt'
path_all_af_indexed = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/merged_hapFIRE_allele_frequency_indexed.csv'

# Parse command line arguments
number = int(sys.argv[1])
train_samples = sys.argv[2].strip('"').split()

# File and directory setup
directory_name = f"split_number_{number}"
os.makedirs(directory_name, exist_ok=True)  # Creates the directory if it doesn't exist


# Here you can add code to use train_samples and test_samples within the newly created directory.
# For example, you might want to save some files related to the train/test data in this directory.
allele_freq = dd.read_csv(path_all_af_indexed, sep = ',', usecols=train_samples)

allele_freq = allele_freq.compute() 

#n_components = len(allele_freq.columns)

#scaler = StandardScaler()
#scaler.fit(allele_freq)
#scaled = scaler.fit_transform(allele_freq)

#scaled_allele_freq = pd.DataFrame(scaled, columns=allele_freq.columns)

# Perform PCA
#pca = PCA(n_components=n_components)  # Set the desired number of components
#pca.fit(scaled_allele_freq)

# Get the transformed data (projected onto the principal components)
#transformed_data = pca.transform(scaled_allele_freq)

#explain_var_ratio = pca.explained_variance_ratio_

#cumulative_variance = np.cumsum(explain_var_ratio)

# Set the desired cumulative explained variance threshold
#threshold = 0.96

# Find the number of components that exceed the threshold
#num_components = np.sum(cumulative_variance <= threshold) + 1

num_components = 2
# Open the file in write mode and write the integer with a newline character at the end
with open(f'split_number_{number}/num_components_full_genome.txt', 'w') as file:
    file.write(f"{num_components}\n")

clim_sites_during_exp = pd.read_csv('/carnegie/nobackup/scratch/tbellagio/grene/data/bioclimvars_experimental_sites_era5.csv')

sites_af = pd.Series(train_samples).str.split('_').str[0].astype(int)

sites_af.name = 'site'

env = sites_af.reset_index().merge(clim_sites_during_exp).drop(['index'],axis=1)

## for now only work with bio1 and bio12 
env = env[['site', 'bio1', 'bio12']]

##scale it
for i in env.columns[1:]:
    env[i] = (env[i] - np.mean(env[i])) / np.std(env[i])

## delete the sites column
env = env.drop('site',axis=1)

env.to_csv(f'split_number_{number}/environment_lea_full_genome.csv', sep = ',', index=False)