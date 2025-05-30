# Get the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# args[1] will be the split number (assuming it is the first argument after the script name)
split_number <- as.integer(args[1])

# args[2] will be the space-separated train samples
train_samples <- args[2]

# Remove the outer brackets and quotes
train_samples <- gsub("\\[|\\]|'", "", train_samples)

# Split the string into a vector of strings
train_samples_vector <- strsplit(train_samples, ",")[[1]]

# Trim any leading or trailing whitespace from each element
train_samples_vector <- trimws(train_samples_vector)

print(train_samples_vector)
# If the identifiers are numeric, convert them to numeric
# This step is necessary if you need to use them as numeric indices

library("LEA")
library("dplyr")
library("qvalue")
library('data.table')

#path_ldp_af = '/carnegie/nobackup/scratch/xwu/grenet/merged_frequency/merged_hapFIRE_allele_frequency_LDpruned.txt'
path_all_af_indexed = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/merged_hapFIRE_allele_frequency_indexed.csv'

results_dir = sprintf('/carnegie/nobackup/scratch/tbellagio/gea_grene-net/leave_1_out/split_number_%d', split_number)

k <- as.integer(readLines(sprintf('%s/num_components_full_genome.txt', results_dir)))
env <- read.csv(sprintf('%s/environment_lea_full_genome.csv', results_dir), header = TRUE, colClasses = "numeric")

rownames(env) <- NULL
colnames(env) <- NULL
env <- as.matrix(env)

output_file_pvalue <- paste0(results_dir, "/pvalues_full_genome.csv")
output_file_zscores <- paste0(results_dir, "/zscores_full_genome.csv")
output_file_gif <- paste0(results_dir, "/gif_full_genome.csv")
output_file_qvalue = paste0(results_dir, "/qvalues_full_genome.csv")
output_file_notadj_pvalues <- paste0(results_dir, "/notadj_pvalues_full_genome.csv")
output_file_notadj_qvalues <- paste0(results_dir, "/notadj_qvalues_full_genome.csv")


#allele_freq <- read.table(path_all_af_indexed, sep = ',', header = TRUE, check.names = FALSE)
#allele_freq <- allele_freq[, train_samples_vector]  # Use dynamic columns
allele_freq <- fread(path_all_af_indexed, select = train_samples_vector, sep = ',', header = TRUE, check.names = FALSE)

colnames(allele_freq) <- NULL
allele_freq <- t(as.matrix(allele_freq))

mod2 <- lfmm2(input = allele_freq, env = env, K = k)
pv <- lfmm2.test(object = mod2, input = allele_freq, env = env)

p_values <- pv$pvalues
zscores <- pv$zscores
gif <- pv$gif

write.csv(p_values, file = output_file_pvalue)
write.csv(zscores, file = output_file_zscores)
write.csv(gif, file = output_file_gif)

##calculate qvalues 
## gif has already been conducted 
#p_values[, 1] <- NULL 
qvalues <- qvalue(p_values)$qvalue
write.csv(qvalues, file = output_file_qvalue)

# Calculate p-values using the chi-squared distribution
# pchisq gives the CDF from 0 to 'q', so use 'lower.tail = FALSE' to get 1 - CDF
notadj_pvalues <- pchisq((zscores^2), df = 1, lower.tail = FALSE)
write.csv(notadj_pvalues, file = output_file_notadj_pvalues)

notadj_qvalues <- qvalue(notadj_pvalues)$qvalue
write.csv(notadj_qvalues, file = output_file_notadj_qvalues)






