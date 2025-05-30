# Get the command line argumentsdelta_p
args <- commandArgs(trailingOnly = TRUE)

# args[1] will be the split number (assuming it is the first argument after the script name)
split_number <- as.integer(args[1])

# args[2] will be the space-separated train samples
train_samples <- args[2]
print(train_samples)

# Remove brackets, single quotes, and newlines
clean_string <- gsub("\\[|\\]|'|\n", "", train_samples)

# Split the cleaned string into a vector
train_samples_vector <- strsplit(clean_string, "\\s+")[[1]]

# Print the result
print(train_samples_vector)

library("lfmm")
library("dplyr")
library("qvalue")
library('data.table')

#path_ldp_af = '/carnegie/nobackup/scratch/xwu/grenet/merged_frequency/merged_hapFIRE_allele_frequency_LDpruned.txt'
#path_all_af_indexed = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/merged_hapFIRE_allele_frequency_indexed.csv'
#delta_p_normed_file = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/leave_1_out_sv/delta_p_normed.csv'
delta_p_file = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/leave_1_out/delta_p_mean.csv'


results_dir = sprintf('/carnegie/nobackup/scratch/tbellagio/gea_grene-net/leave_1_out/split_number_%d', split_number)

#k <- as.integer(readLines(sprintf('%s/num_components_full_genome.txt', results_dir)))
k=2


#get the environemntal varialbes to trian the model on 
env_train <- read.csv(sprintf('%s/env_training_sites.csv', results_dir), header = TRUE, colClasses = "numeric")
env_train <- env_train[, c('bio1', 'bio12')]
rownames(env_train) <- NULL
colnames(env_train) <- NULL
env_train <- as.matrix(env_train)
print(dim(env_train))
## get the environmental variables fro testing 
env_test <- read.csv(sprintf('%s/env_testing_sites.csv', results_dir), header = TRUE, colClasses = "numeric")
bio1 = env_test[, 'bio1']
bio1 <- as.matrix(bio1)
bio12 = env_test[, 'bio12']
bio12 <- as.matrix(bio12)

## get the allele freq 
allele_freq <- fread(delta_p_file, select = train_samples_vector, sep = ',', header = TRUE, check.names = FALSE)
print(dim(allele_freq))
colnames(allele_freq) <- NULL
allele_freq <- t(as.matrix(allele_freq))


### create output files 
output_file_wo_calibration_pvalue <- paste0(results_dir, "/delta_p_mean/wo_calibration_pvalue_full_genome.csv")
output_file_w_calibration_pvalue <- paste0(results_dir, "/delta_p_mean/w_calibration_pvalue_full_genome.csv")

output_file_effect_sizes_simple <- paste0(results_dir, "/delta_p_mean/effect_sizes_simple_full_genome.csv")

output_file_w_calibration_qvalue = paste0(results_dir, "/delta_p_mean/w_calibration_qvalue_full_genome.csv")
output_file_wo_calibration_qvalue = paste0(results_dir, "/delta_p_mean/wo_calibration_qvalue_full_genome.csv")

output_af_bio1 = paste0(results_dir,'/delta_p_mean/predicted_allele_freq_bio1.csv')
output_af_bio12 = paste0(results_dir,'/delta_p_mean/predicted_allele_freq_bio12.csv')

output_indices_bio1 = paste0(results_dir,'/delta_p_mean/snps_indices_bio1.csv')
output_indices_bio12 = paste0(results_dir,'/delta_p_mean/snps_indices_bio12.csv')

## train teh model
mod.lfmm <- lfmm_ridge(Y = allele_freq,
                       X = env_train,
                       K = k)

## Perform association testing using the fitted model:
pv <- lfmm_test(Y = allele_freq,
                X = env_train,
                lfmm = mod.lfmm,
                calibrate = "gif")


## save the p values 
p_values <- pv$pvalue
write.csv(p_values, file = output_file_wo_calibration_pvalue)

## save the calibrated p values  
calibrated_pvalues <- pv$calibrated.pvalue
write.csv(calibrated_pvalues, file = output_file_w_calibration_pvalue)

## save the betas 
beta_i <- pv$B
write.csv(beta_i, file = output_file_effect_sizes_simple)

## save the q values without calibration 
q_qvalues_wo_calibration <- qvalue(p_values)$qvalue
write.csv(q_qvalues_wo_calibration, file = output_file_wo_calibration_qvalue)

## save the q values with calibration 
qvalues_w_calibration <- qvalue(calibrated_pvalues)$qvalue
write.csv(qvalues_w_calibration, file = output_file_w_calibration_qvalue)

## predict 
## save the af mean to then reescale , use as intersection 
af_means <- colMeans(allele_freq)

## get hte betas for the bio 1 adn bio12
bio1_beta = matrix(pv$B[,1])
bio12_beta = matrix(pv$B[,2])

## get the indices for the singifican p values 
indices_bio1 = which(pv$pvalue[,1] < 1.5453657571674064e-08)
indices_bio12 = which(pv$pvalue[,2] < 1.5453657571674064e-08)

write.csv(indices_bio1, file = output_indices_bio1)
write.csv(indices_bio12, file = output_indices_bio12)

## adding the mean allele freq as an intercept
## in r thisis how you do it 
pred_af_bio1 = bio1 %*% t(bio1_beta) #+ rep(af_means, each = nrow(bio1))
pred_af_bio12 = bio12 %*% t(bio12_beta) # + rep(af_means, each = nrow(bio12))

pred_af_filt_bio1 = pred_af_bio1[,indices_bio1]
pred_af_filt_bio12 = pred_af_bio12[,indices_bio12]

write.csv(pred_af_filt_bio1, file = output_af_bio1)
write.csv(pred_af_filt_bio12, file = output_af_bio12)

## this has to be done once per env var 
#bio1 = as.matrix(env[,1])
#bio12 = as.matrix(env[,2])

#prs_effect_size_bio1 <- effect_size(Y = allele_freq, X = bio1, lfmm.object = mod.lfmm) 
#write.csv(prs_effect_size_bio1, file = output_file_effect_sizes_prs_bio1)

#prs_effect_size_bio12 <- effect_size(Y = allele_freq, X = bio12, lfmm.object = mod.lfmm) 
#write.csv(prs_effect_size_bio12, file = output_file_effect_sizes_prs_bio12)




