#!/usr/bin/env Rscript

# ============================================================================
# LFMM Analysis for Multiple Bioclimatic Variables (Last Generation)
# ============================================================================
#
# This script performs LFMM (Latent Factor Mixed Model) analysis for a given
# bioclimatic variable using a specified number of components (K) for the last
# generation data.
#
# Required Input Files:
# 1. ../data/delta_p_maf05_mincount05_lastgensamples.csv
#    - Allele frequency data for the last generation
# 2. env_train_files/env_train_*.csv
#    - Environmental training data for each bioclimatic variable
#
# Output Files:
# 1. lfmm_fullresults_all_k/w_calibration_pvalue_full_genome_*.csv
#    - Calibrated p-values
# 2. lfmm_fullresults_all_k/effect_sizes_simple_full_genome_*.csv
#    - Effect sizes
#

# --- Import Required Libraries ---
library("lfmm")
library("dplyr")
library("qvalue")
library('data.table')

# --- Get Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
samples <- args[1]  # Space-separated list of samples
bio_name <- args[2]  # Bioclimatic variable name
k <- as.integer(args[3])  # Number of components

# Clean up samples argument
samples <- gsub("\\[|\\]|'", "", samples)
samples_vector <- strsplit(samples, " ")[[1]]
samples_vector <- trimws(samples_vector)

# --- Set Up Paths ---
analysis_dir = '.'
results_dir = 'lfmm_fullresults_all_k'
delta_p_file = '../data/delta_p_maf05_mincount05_lastgensamples.csv'

# --- Load Data ---
# Load environmental training data
env_train_file <- sprintf('%s/env_train_files/env_train_%s.csv', analysis_dir, bio_name)
env_train <- read.csv(env_train_file, header = TRUE, colClasses = "numeric")
print(dim(env_train))

# Load allele frequency data
allele_freq <- fread(delta_p_file, select = samples_vector, sep = ',', header = TRUE, check.names = FALSE)
colnames(allele_freq) <- NULL
allele_freq <- t(as.matrix(allele_freq))
print(dim(allele_freq))

# --- Set Up Output Files ---
output_file_w_calibration_pvalue <- paste0(results_dir, '/w_calibration_pvalue_full_genome_', bio_name, '_k', k, '.csv')
output_file_wo_calibration_pvalue <- paste0(results_dir, '/wo_calibration_pvalue_full_genome_', bio_name, '_k', k, '.csv')
output_file_effect_sizes_simple <- paste0(results_dir, '/effect_sizes_simple_full_genome_', bio_name, '_k', k, '.csv')

# --- Run LFMM Analysis ---
# Train the model
mod.lfmm <- lfmm_ridge(Y = allele_freq, X = env_train, K = k)

# Perform association testing
pv <- lfmm_test(Y = allele_freq, X = env_train, lfmm = mod.lfmm, calibrate = "gif")

# --- Save Results ---
# Save p-values
p_values <- pv$pvalue
write.csv(p_values, file = output_file_wo_calibration_pvalue)

# Save calibrated p-values
calibrated_pvalues <- pv$calibrated.pvalue
write.csv(calibrated_pvalues, file = output_file_w_calibration_pvalue)

# Save effect sizes
beta_i <- pv$B
write.csv(beta_i, file = output_file_effect_sizes_simple)

# Beta coefficients for the biovariable
#biovar_beta = matrix(pv$B[, bio_name])

# Indices for significant p-values
#indices_biovar = which(pv$pvalue[, bio_name] < 1.5453657571674064e-08)
#write.csv(indices_biovar, file = output_indices_biovar)