library("lfmm")
library("dplyr")
library("qvalue")
library('data.table')

# Get the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# args[2] will be the space-separated train samples
samples <- args[1]

# Remove the outer brackets and quotes
samples <- gsub("\\[|\\]|'", "", samples)

# Split the string into a vector of strings
samples_vector <- strsplit(samples, " ")[[1]]

# Trim any leading or trailing whitespace from each element
samples_vector <- trimws(samples_vector)

print(samples_vector)
# args[1] will be the split number (assuming it is the first argument after the script name)
# args[2] will be the biovariable name, e.g., 'bio12'
bio_name <- args[2]  # capture the biovariable name

k <- as.integer(args[3])
print(k)
# /global/scratch/users/tbellg/gea_grene-net/lfmm_full
# Dynamic paths and other initialization
analysis_dir = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/lfmm_full'

results_dir = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/lfmm_full/lfmm_fullresults_all_k'
delta_p_file = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/key_files/delta_p_maf05_mincount05_firstgensamples.csv'

# Load training environmental variables dynamically based on the bio_name
env_train_file <- sprintf('%s/env_train_files/env_train_%s.csv', analysis_dir, bio_name)
env_train <- read.csv(env_train_file, header = TRUE, colClasses = "numeric")
print(dim(env_train))
# Load allele frequencies
allele_freq <- fread(delta_p_file, select = samples_vector, sep = ',', header = TRUE, check.names = FALSE)
colnames(allele_freq) <- NULL
allele_freq <- t(as.matrix(allele_freq))

output_file_w_calibration_pvalue <- paste0(results_dir, '/w_calibration_pvalue_full_genome_', bio_name, '_k', k, '.csv')
output_file_wo_calibration_pvalue <- paste0(results_dir, '/wo_calibration_pvalue_full_genome_', bio_name, '_k', k, '.csv')
output_file_effect_sizes_simple <- paste0(results_dir, '/effect_sizes_simple_full_genome_', bio_name, '_k', k, '.csv')

#output_indices_biovar = paste0(results_dir,'/snps_indices_', bio_name, '.csv')


print(dim(allele_freq))
print(dim(env_train))
# Model training
mod.lfmm <- lfmm_ridge(Y = allele_freq, X = env_train, K = k)

# Association testing
pv <- lfmm_test(Y = allele_freq, X = env_train, lfmm = mod.lfmm, calibrate = "gif")

## save the p values
p_values <- pv$pvalue
write.csv(p_values, file = output_file_wo_calibration_pvalue)

## save the p values
calibrated_pvalues <- pv$calibrated.pvalue
write.csv(calibrated_pvalues, file = output_file_w_calibration_pvalue)

## save the betas
beta_i <- pv$B
write.csv(beta_i, file = output_file_effect_sizes_simple)

# Beta coefficients for the biovariable
#biovar_beta = matrix(pv$B[, bio_name])

# Indices for significant p-values
#indices_biovar = which(pv$pvalue[, bio_name] < 1.5453657571674064e-08)
#write.csv(indices_biovar, file = output_indices_biovar)