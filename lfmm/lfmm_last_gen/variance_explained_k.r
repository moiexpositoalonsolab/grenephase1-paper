#!/usr/bin/env Rscript

# ============================================================================
# Variance Explained Analysis for LFMM Model (Last Generation)
# ============================================================================
#
# This script analyzes how much variance is explained by different numbers of 
# components (K) in the LFMM model for the last generation data to determine 
# the optimal K value for the subsequent analysis.
#
# Required Input Files:
# 1. ../data/delta_p_maf05_mincount05_lastgensamples.csv
#    - Allele frequency data for the last generation
# 2. env_train_files/env_train_bio1.csv
#    - Environmental training data for bio1 (used as reference)
#
# Output Files:
# 1. lfmm_fullresults_all_k/variance_explained_k.csv
#    - CSV file containing variance explained for different K values
# 2. lfmm_fullresults_all_k/variance_explained_k.pdf
#    - PDF plot showing variance explained vs K
#

# --- Import Required Libraries ---
library("lfmm")
library("dplyr")
library("qvalue")
library('data.table')
library("ggplot2")

# --- Set Up Paths and Parameters ---
analysis_dir = '.'
results_dir = 'lfmm_fullresults_all_k'
delta_p_file = '../data/delta_p_maf05_mincount05_lastgensamples.csv'

# --- Load and Prepare Data ---
# Read environmental training data for bio1
env_train_file <- sprintf('%s/env_train_files/env_train_bio1.csv', analysis_dir)
env_train <- read.csv(env_train_file, header = TRUE, colClasses = "numeric")
print(dim(env_train))

# Read allele frequency data
allele_freq <- fread(delta_p_file, sep = ',', header = TRUE, check.names = FALSE)
colnames(allele_freq) <- NULL
allele_freq <- t(as.matrix(allele_freq))
print(dim(allele_freq))

# --- Define K Values to Test ---
# Test K values from 1 to 10
k_values <- 1:10

# --- Calculate Variance Explained for Each K ---
# Initialize vector to store variance explained values
variance_explained <- numeric(length(k_values))

# Calculate variance explained for each K
for (i in seq_along(k_values)) {
    k <- k_values[i]
    print(paste("Testing K =", k))
    
    # Run LFMM model
    mod.lfmm <- lfmm_ridge(Y = allele_freq, X = env_train, K = k)
    
    # Calculate variance explained
    variance_explained[i] <- mod.lfmm$var_explained
}

# --- Save Results ---
# Create data frame with results
results_df <- data.frame(K = k_values, Variance_Explained = variance_explained)

# Save to CSV
write.csv(results_df, 
          file = paste0(results_dir, '/variance_explained_k.csv'),
          row.names = FALSE)

# --- Create and Save Plot ---
# Create plot of variance explained vs K
p <- ggplot(results_df, aes(x = K, y = Variance_Explained)) +
    geom_line() +
    geom_point() +
    theme_minimal() +
    labs(x = "Number of Components (K)",
         y = "Variance Explained",
         title = "Variance Explained by Different K Values (Last Generation)")

# Save plot
ggsave(paste0(results_dir, '/variance_explained_k.pdf'), p, width = 10, height = 6)



