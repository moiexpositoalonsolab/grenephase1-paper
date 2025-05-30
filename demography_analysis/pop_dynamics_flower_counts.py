#!/usr/bin/env python
# coding: utf-8

"""
Population Dynamics Analysis - Flower Counts
This script analyzes the number of flowers collected across three generations at different experimental sites.
It correlates flower counts with environmental variables (bio1 - annual mean temperature) and generates
visualizations showing temporal trends in flower production.

Input:
- SURVIVAL_total_flowers_collected.csv: Raw flower count data
- bioclimvars_experimental_sites_era5.csv: Climate data
- generation_1_parallelism.txt: Parallelism data

Output:
- Various plots showing flower count dynamics across generations
- Statistical analysis of trends and correlations
"""

# --- Import Libraries ---
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from scipy import stats

# Set matplotlib parameters for PDF output
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

# --- Load and Prepare Data ---
# Load flower count data
survival = pd.read_csv('../../grenenet_cleaning/grenenet/SURVIVAL_total_flowers_collected.csv')

# Load climate data (annual mean temperature - bio1)
climate = pd.read_csv('../key_files/bioclimvars_experimental_sites_era5.csv')[['site', 'bio1']]

# Load parallelism data
parallelism = pd.read_csv('../key_files/generation_1_parallelism.txt', sep='\t')

# --- Data Cleaning and Preparation ---
# Select relevant columns and create site_plot identifier
survival = survival[['site', 'plot', '1_flowerstotal','2_flowerstotal','3_flowerstotal']]
survival['site_plot'] = survival['site'].astype(str) + '_' + survival['plot'].astype(str)
survival = survival[['site_plot', '1_flowerstotal','2_flowerstotal','3_flowerstotal']]
survival.columns = ['site_plot', 1, 2, 3]

# Convert to long format and clean data
df = survival.set_index('site_plot')
df_cleaned = df.dropna(thresh=2)  # Keep rows with at least 2 non-NaN values
df_cleaned['site'] = df_cleaned.index.str.split('_').str[0]

# Prepare climate data
climate['site'] = climate['site'].astype(str)
grouped_sites_with_bio1 = pd.merge(df_cleaned[['site']], climate, on='site').drop_duplicates()
grouped_sites_with_bio1_sorted = grouped_sites_with_bio1.sort_values(by='bio1')

# --- Analysis Functions ---
def fit_poly(group_long):
    """
    Fit a second-degree polynomial regression to the data.
    
    Args:
        group_long: DataFrame with 'Generation' and 'Population' columns
        
    Returns:
        model: Fitted polynomial regression model
        x_fine: X values for plotting
        y_fine: Predicted Y values for plotting
    """
    x = group_long['Generation'].values.reshape(-1, 1)
    y = group_long['Population'].values
    
    poly = PolynomialFeatures(degree=2)
    x_poly = poly.fit_transform(x)
    
    model = LinearRegression()
    model.fit(x_poly, y)
    
    x_fine = np.linspace(min(x), max(x), 100).reshape(-1, 1)
    x_fine_poly = poly.transform(x_fine)
    y_fine = model.predict(x_fine_poly)
    
    return model, x_fine, y_fine

def fit_lr(group_long):
    """
    Fit a linear regression to the data.
    
    Args:
        group_long: DataFrame with 'Generation' and 'Population' columns
        
    Returns:
        model: Fitted linear regression model
        x_fine: X values for plotting
        y_fine: Predicted Y values for plotting
    """
    x = group_long['Generation'].values.reshape(-1, 1)
    y = group_long['Population'].values
    
    model = LinearRegression().fit(x, y)
    x_fine = np.linspace(min(x), max(x), 100).reshape(-1, 1)
    y_fine = model.predict(x_fine)
    
    return model, x_fine, y_fine

# --- Visualization Setup ---
# Set up color palette for temperature-based visualization
colors = sns.color_palette('coolwarm', n_colors=27)

# Set up plot parameters
s_size = 60
linewidth = 5
alpha_scatter = 0.5
n_rows = 1
n_cols = 27
plt.rcParams['axes.axisbelow'] = True

# Create figure and axes
fig, axes = plt.subplots(n_rows, n_cols, figsize=(30, 4), sharex=True, sharey=True)
axes = axes.flatten()

# Initialize dictionaries for storing results
site_min_in_last_gen = {}
coeff_sites = {}

# --- Main Analysis Loop ---
for i, site in enumerate(grouped_sites_with_bio1_sorted['site']):
    group = df_cleaned[df_cleaned['site'] == site]
    group_long = group.melt(var_name='Generation', value_name='Population', value_vars=[1, 2, 3], ignore_index=False)
    group_long = group_long.dropna()
    group_long['Generation'] = pd.to_numeric(group_long['Generation'])
    group_long['Population'] = pd.to_numeric(group_long['Population'])

    ax = axes[i] if len(axes) > 1 else axes

    if group_long['Generation'].nunique() == 3:
        model, x_fine, y_fine = fit_poly(group_long)
        intercept = model.intercept_
        coefficients = model.coef_
        coeff_sites[site] = [intercept, coefficients]
        
        if (coefficients[1] < 0) & (coefficients[2] > 0):
            # Evolutionary rescue pattern
            sns.scatterplot(data=group_long, x='Generation', y='Population', ax=ax, 
                          color=colors[i], alpha=0.7, s=s_size, edgecolor=None)
            ax.plot(x_fine, y_fine, color=colors[i], linewidth=linewidth)
        else:
            sns.scatterplot(data=group_long, x='Generation', y='Population', ax=ax, 
                          color=colors[i], alpha=0.2, s=s_size, edgecolor=None)
            ax.plot(x_fine, y_fine, color=colors[i], linewidth=linewidth, alpha=0.2)

    elif group_long['Generation'].nunique() < 3:
        model, x_fine, y_fine = fit_lr(group_long)
        intercept = model.intercept_
        coefficients = model.coef_
        coeff_sites[site] = [intercept, coefficients]
        
        if coefficients[0] > 0:
            # Positive trend
            sns.scatterplot(data=group_long, x='Generation', y='Population', ax=ax, 
                          color=colors[i], alpha=0.2, s=s_size, edgecolor=None)
            ax.plot(x_fine, y_fine, color=colors[i], linewidth=linewidth, alpha=0.2)
        else:
            # Negative trend
            sns.scatterplot(data=group_long, x='Generation', y='Population', ax=ax, 
                          color=colors[i], alpha=0.2, s=s_size, edgecolor=None)
            ax.plot(x_fine, y_fine, color=colors[i], linewidth=linewidth, alpha=0.2)

    # Clean up plot appearance
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    # Add labels and grid
    bio1_value = int(grouped_sites_with_bio1_sorted.loc[grouped_sites_with_bio1_sorted["site"] == site, "bio1"].values[0])
    ax.set_title(f'Site {site}')
    ax.set_xlabel('Generation')
    ax.set_ylabel('Population Size')
    ax.grid(True, color='lightgrey', alpha=0.7)
    ax.tick_params(axis='x', which='both', bottom=False)
    ax.tick_params(axis='y', which='both', left=False)
    ax.set_xlabel('')

# Add common x-axis label
fig.supxlabel('Generations')

# Adjust layout
plt.tight_layout()
plt.subplots_adjust(wspace=0.1, hspace=0.1)

# Save plots
plt.savefig('pop_dynamics_highlighted_er.png')
plt.savefig('pop_dynamics_highlighted_er.pdf')

# --- Process Results ---
# Combine site data
combined_sites_data = {site: [intercept] + coefficients.tolist() 
                      for site, (intercept, coefficients) in coeff_sites.items()}

# Clean up coefficients
for site, values in combined_sites_data.items():
    if len(values) == 4 and values[1] == 0:
        combined_sites_data[site] = [values[0]] + values[2:]
    elif len(values) == 2:
        combined_sites_data[site] = values + [0]

# Create results DataFrame
coeff_sitesdf = pd.DataFrame(combined_sites_data).T.reset_index()
coeff_sitesdf.columns = ['site', 'intercept', 'linear_term', 'cuadratic_term']
coeff_sitesdf = coeff_sitesdf.merge(climate)

# Process parallelism data
parallelism['site'] = parallelism['site'].astype(str)
parallelism_snps = parallelism[parallelism['source'] == 'snp'][['site', 'mean']]
parallelism_snps.columns = ['site', 'p_snps']
parallelism_ecotypes = parallelism[parallelism['source'] == 'ecotype'][['site', 'mean']]
parallelism_ecotypes.columns = ['site', 'p_ecotypes']

# Merge all data
coeff_sitesdf = coeff_sitesdf.merge(parallelism_snps).merge(parallelism_ecotypes)
coeff_sitesdf = coeff_sitesdf[coeff_sitesdf['site'] != '57']
coeff_sitesdf = coeff_sitesdf[coeff_sitesdf['cuadratic_term'] != 0]

# --- Create Final Plots ---
# Plot 1: Quadratic term vs Repeatability
plt.figure(figsize=(9, 6))
x = coeff_sitesdf['p_ecotypes']
y = coeff_sitesdf['cuadratic_term']
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

sns.regplot(data=coeff_sitesdf, x='p_ecotypes', y='cuadratic_term', 
           scatter=False, line_kws={'color': 'grey'})

sns.scatterplot(data=coeff_sitesdf, x='p_ecotypes', y='cuadratic_term', 
               hue='bio1', palette='coolwarm', edgecolor=None, s=150)

# Add site labels
for i in range(coeff_sitesdf.shape[0]):
    site = coeff_sitesdf['site'].iloc[i]
    x_value = coeff_sitesdf['p_ecotypes'].iloc[i]
    y_value = coeff_sitesdf['cuadratic_term'].iloc[i]
    plt.text(x_value, y_value, site, fontsize=12, ha='right', color='dimgray')

# Add statistics
plt.text(min(x), max(y), f'R²: {r_value**2:.2f} pvalue: {p_value:.4f}', 
         fontsize=12, color='dimgray')

# Style plot
plt.ylabel('Quadratic term')
plt.xlabel('Repeatability 1st gen')
plt.rc('font', family='sans-serif', size=12, weight='normal')
plt.rc('axes', titlesize=12, labelsize=12)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)

# Clean up plot appearance
for spine in plt.gca().spines.values():
    spine.set_visible(False)
plt.tick_params(axis='both', colors='#4D4D4D')
plt.grid(True, color='lightgrey', alpha=0.7)

# Save final plot
plt.tight_layout()
plt.savefig('parallelism_ecotypes_cuadratic_term_flower_counts.png')
plt.savefig('parallelism_ecotypes_cuadratic_term_flower_counts.pdf')

plt.show()
