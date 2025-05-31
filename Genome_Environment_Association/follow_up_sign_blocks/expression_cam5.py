#!/usr/bin/env python
# coding: utf-8

"""
Analysis of CAM5 Isoform Expression and Ecotype Frequency Dynamics

This script analyzes the expression levels of two CAM5 gene isoforms (AT2G27030.1 and AT2G27030.3)
in relation to climate (annual mean temperature, bio1). It performs linear regression to assess
these relationships. Additionally, the script visualizes the change in frequency of different
ecotypes across generations for selected experimental sites using Plotly.

Inputs:
- ../key_files/transcptomic_data_meta.csv: Metadata for transcriptomic samples.
- ../key_files/transcriptomics_data_100t.csv: Transcriptomic expression data.
- ../key_files/founder_ecotype_frequency.txt: List of founder ecotypes.
- ../key_files/merged_ecotype_frequency.txt: Merged ecotype frequency data.
- ../key_files/bioclimvars_experimental_sites_era5.csv: ERA5 climate data for experimental sites.

Outputs:
- Interactive plots displaying the relationship between CAM5 isoform expression and bio1.
- Interactive plots showing ecotype frequency change across generations for selected sites/plots.
- A pandas DataFrame (results_df) containing linear regression results for CAM5 expression vs. bio1.
"""

# --- Import Libraries ---
import pandas as pd
import seaborn as sns
import numpy as np
import statsmodels.api as sm
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px

# --- Load Data ---
# Load transcriptomic metadata and data
meta = pd.read_csv('../key_files/transcptomic_data_meta.csv')
t = pd.read_csv('../key_files/transcriptomics_data_100t.csv')
t = t.drop('Unnamed: 0', axis=1)

# Load founder ecotype list
grenenet_ecotypes = pd.read_csv('../key_files/founder_ecotype_frequency.txt', sep = '\t', header=None)[0]

# Load climate data (annual mean temperature)
bio1 = meta.set_index('index')['worldClim_bio.1']

# Load merged ecotype frequency data
ecotype_freq_raw = pd.read_csv('../key_files/merged_ecotype_frequency.txt', sep = '\t')

# --- CAM5 Expression Analysis ---
# Identify CAM5 isoforms and extract expression data
cam5_isoforms = [i for i in t.columns if 'AT2G27030' in i]
cam5_expresion = t.set_index('id')[cam5_isoforms]

# Combine CAM5 expression with bio1 climate data
cam5_expresion_bio1 = pd.concat([cam5_expresion, bio1], axis=1).dropna()

# Perform linear regression for each CAM5 isoform against bio1
print("Performing linear regression for CAM5 isoform expression vs. bio1...")

x1 = cam5_expresion_bio1['worldClim_bio.1']
y1 = cam5_expresion_bio1['AT2G27030.1']
X1 = sm.add_constant(x1)
model1 = sm.OLS(y1, X1, missing='drop').fit()

x2 = cam5_expresion_bio1['worldClim_bio.1']
y2 = cam5_expresion_bio1['AT2G27030.3']
X2 = sm.add_constant(x2)
model2 = sm.OLS(y2, X2, missing='drop').fit()

# Store regression results in a DataFrame
results = {
    'Model': ['AT2G27030.1', 'AT2G27030.3'],
    'Slope': [model1.params['worldClim_bio.1'], model2.params['worldClim_bio.1']],
    'p-value': [model1.pvalues['worldClim_bio.1'], model2.pvalues['worldClim_bio.1']],
    'R-squared': [model1.rsquared, model2.rsquared],
    'Number of samples': [int(model1.nobs), int(model2.nobs)] # Ensure nobs is int for DataFrame
}
results_df = pd.DataFrame(results)

print("Linear regression results:")
print(results_df)

# Filter for Grenenet ecotypes and reset index for plotting
cam5_expresion_bio1_grenenet = cam5_expresion_bio1.reset_index()
cam5_expresion_bio1_grenenet.columns = ['ecotype', 'AT2G27030.1', 'AT2G27030.3', 'bio1']
cam5_expresion_bio1_grenenet = cam5_expresion_bio1_grenenet[cam5_expresion_bio1_grenenet['ecotype'].isin(grenenet_ecotypes)]

# Generate regression plots for Grenenet ecotypes
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
sns.regplot(data=cam5_expresion_bio1_grenenet, x='bio1', y='AT2G27030.1')
plt.title('AT2G27030.1 Expression vs. bio1 (Grenenet Ecotypes)')
plt.xlabel('bio1 (Annual Mean Temperature)')
plt.ylabel('Expression Level')

plt.subplot(1, 2, 2)
sns.regplot(data=cam5_expresion_bio1_grenenet, x='bio1', y='AT2G27030.3')
plt.title('AT2G27030.3 Expression vs. bio1 (Grenenet Ecotypes)')
plt.xlabel('bio1 (Annual Mean Temperature)')
plt.ylabel('Expression Level')

plt.tight_layout()
plt.show() # Display interactive plots

# --- Ecotype Frequency Dynamics Visualization ---
print("Visualizing ecotype frequency dynamics...")

# Define hot sites for filtering (example: site 13)
hot_sites_str = ['13'] # List of site numbers as strings

# Select columns corresponding to samples from hot sites in the raw ecotype frequency data
hot_sites_samples = [col for col in ecotype_freq_raw.columns if col.startswith(tuple(hot_sites_str))]
ecotype_freq = ecotype_freq_raw[hot_sites_samples].T.reset_index()

# Extract generation, plot, and site from sample names
ecotype_freq['generation'] = ecotype_freq['index'].str.split('_').str[1].astype(int)
ecotype_freq['plot'] = ecotype_freq['index'].str.split('_').str[2]
ecotype_freq['site'] = ecotype_freq['index'].str.split('_').str[0]

# Melt the DataFrame to long format for plotting
ecotype_freq = ecotype_freq.drop(columns=['index']).melt(id_vars = ['plot','site', 'generation'])
ecotype_freq.columns = ['plot','site', 'generation', 'ecotype', 'freq']

# Define unique sites, plots, and ecotypes for plotting
unique_sites_plot = hot_sites_str
unique_plots = ecotype_freq['plot'].unique()
unique_ecotypes = ecotype_freq['ecotype'].unique()

# Get and extend color palette for ecotypes
base_colors = px.colors.qualitative.Plotly
extended_colors = [base_colors[i % len(base_colors)] for i in range(len(unique_ecotypes))]
color_mapping = {ecotype: extended_colors[i] for i, ecotype in enumerate(unique_ecotypes)}

# Create subplots for each site and plot combination
fig = make_subplots(
    rows=len(unique_sites_plot), 
    cols=len(unique_plots), 
    shared_xaxes=True, 
    shared_yaxes=True, 
    subplot_titles=[f"Site {site} - Plot {plot}" for site in unique_sites_plot for plot in unique_plots]
)

# Create a dictionary to track whether the legend for each ecotype has been shown
legend_shown = {ecotype: False for ecotype in unique_ecotypes}

# Iterate over each site and plot to create traces for each ecotype's frequency change
for row_idx, site in enumerate(unique_sites_plot):
    for col_idx, plot in enumerate(unique_plots):
        # Filter data for the current site and plot
        plot_data = ecotype_freq[(ecotype_freq['site'] == site) & (ecotype_freq['plot'] == plot)].sort_values(by='generation')
   
        # Iterate over each unique ecotype in the current filtered data
        for ecotype in unique_ecotypes:
            # Filter data for the current ecotype
            ecotype_data = plot_data[plot_data['ecotype'] == ecotype]
            
            if not ecotype_data.empty:
                # Only show the legend once per ecotype
                show_legend = not legend_shown[ecotype]
                legend_shown[ecotype] = True

                # Add line trace for the current ecotype's frequency across generations
                fig.add_trace(
                    go.Scatter(
                        x=ecotype_data['generation'], 
                        y=ecotype_data['freq'], 
                        mode='lines+markers',
                        name=f'Ecotype {ecotype}', 
                        line=dict(color=color_mapping[ecotype]),
                        legendgroup=ecotype,
                        showlegend=show_legend
                    ),
                    row=row_idx + 1, col=col_idx + 1
                )

# Update layout for ecotype frequency plot
fig.update_layout(
    height=500 * len(unique_sites_plot), 
    title_text="Change in Ecotype Frequency Across Generations for Each Site and Plot",
    xaxis_title="Generation",
    yaxis_title="Frequency",
    legend_title="Ecotype"
)

fig.show() # Display interactive plots
