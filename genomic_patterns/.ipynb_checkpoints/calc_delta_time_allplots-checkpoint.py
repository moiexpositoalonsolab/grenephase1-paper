import pandas as pd
import statsmodels.api as sm
import json
from statsmodels.formula.api import ols
import argparse
import logging

# Set up logging

# Set up the argument parser
parser = argparse.ArgumentParser(description='Process some site and plot identifiers.')

# Add arguments to the parser
parser.add_argument('site', type=str, help='The site identifier to process')

# Parse the command line arguments
args = parser.parse_args()

# Extract the site and plot values from the arguments
site = args.site
print(site)


logging.basicConfig(filename=f'log_site_{site}.log', level=logging.DEBUG)

logging.debug(f'Starting script for site: {site}')



flowers = pd.read_csv('../key_files/merged_sample_table.csv')
num_flowers_map = flowers.set_index('sample_name')['total_flower_counts'].to_dict()
## ge thte mask for the 0.05maf filtering 
var_pos = pd.read_csv('../key_files/var_pos_grenenet.csv')
mask = var_pos['maf05filter'].notna()


## get the samples i need 
merged_hapFIRE_allele_frequency = pd.read_csv('../key_files/merged_hapFIRE_allele_frequency.txt',nrows=1, sep = '\t')
merged_hapFIRE_allele_frequency = merged_hapFIRE_allele_frequency.T
merged_hapFIRE_allele_frequency = merged_hapFIRE_allele_frequency.reset_index()
site_plot_samples = [i for i in merged_hapFIRE_allele_frequency['index'] if i.startswith(str(site) + '_')]
print(site_plot_samples)
## import the dataset 
logging.debug(f'Filtering completed for site {site_plot_samples}')

merged_hapFIRE_allele_frequency = pd.read_csv('../key_files/merged_hapFIRE_allele_frequency.txt', sep = '\t', usecols = site_plot_samples)
## filter

logging.debug(f'Filtering completed for site {site}, shape: {merged_hapFIRE_allele_frequency.shape}')

merged_hapFIRE_allele_frequency = merged_hapFIRE_allele_frequency[mask]
# reset index after filtering 
merged_hapFIRE_allele_frequency = merged_hapFIRE_allele_frequency.reset_index(drop=True)


## generate the allale counts
num_flowers_map = flowers.set_index('sample_name')['total_flower_counts'].to_dict()

allele_counts = {}
for i in merged_hapFIRE_allele_frequency.columns:
    num_flowers = num_flowers_map[i]
    allele_counts[i] = merged_hapFIRE_allele_frequency[i] * num_flowers * 2

allele_counts = pd.concat(allele_counts,axis=1)
allele_counts = allele_counts.T
allele_counts = allele_counts.reset_index()


allele_counts['gen'] = allele_counts['index'].str.split('_').str[1].astype(int)
allele_counts['plot'] = allele_counts['index'].str.split('_').str[2].astype(int)
allele_counts = allele_counts.drop('index',axis=1)
allele_counts = allele_counts.melt(id_vars=['gen', 'plot'])
allele_counts.columns = ['gen', 'plot', 'snp', 'count']
allele_counts = allele_counts.drop('plot',axis=1)
print(len(allele_counts))

logging.debug(f'Filtering completed for site {site}, shape: {allele_counts.shape}')


# Open the file in append mode to add to existing allele_counts without overwriting
with open(f'results_site_{site}.jsonl', 'a') as file:
    buffer = []
    for snp, group in allele_counts.groupby('snp'):
        if len(group) > 1:  # Ensure there's enough data for regression
            model = ols('count ~ gen', data=group).fit()
            slope = model.params['gen'].item()  # Convert numpy float64 to Python float
            p_value = model.pvalues['gen'].item()  # Convert numpy float64 to Python float
            
            # Prepare the result as a JSON object
            result = {
                "snp": str(snp),
                "slope": slope,
                "p_value": p_value,
            }
            
            # Append result to buffer
            buffer.append(result)
            
            # Check if buffer has reached the chunk size of 100
            if len(buffer) == 100:
                # Write all buffered items to the file as JSON Lines
                for item in buffer:
                    file.write(json.dumps(item) + '\n')
                # Clear the buffer
                buffer.clear()
    # Write all results to the file as JSON Lines
    if buffer:
        for item in buffer:
            file.write(json.dumps(item) + '\n')




