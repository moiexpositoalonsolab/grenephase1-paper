.libPaths("/home/tbellagio/miniforge3/envs/r-environment")

library(devtools) # Needed to install the PicMin package from Tom's GitHub install_github("TBooker/PicMin", force = T)

.libPaths("/home/tbellagio/miniforge3/envs/r-environment/lib/R/library")

library(PicMin)
library(tidyverse)
library(poolr)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
numReps <- as.numeric(args[1])

# Create an empty list to store the data frames
data_list <- list()

# Loop through the 12 CSV files and read them into the list, keeping only 'gene' and 'Z_pVal'
for (run in 0:11) {
  file_name <- paste0("wza_run", run, "_results.csv")
  temp_df <- read.csv(file_name)
  
  # Keep only 'gene' and 'Z_pVal', rename 'Z_pVal' to indicate the run
  temp_df <- temp_df[, c("gene", "Z_pVal")]
  colnames(temp_df)[2] <- paste0("Z_pVal_run", run)
  
  data_list[[run + 1]] <- temp_df
}

# Merge all data frames by the 'gene' column
all_lins_p <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), data_list)

# Set the 'gene' column as row names
row.names(all_lins_p) <- all_lins_p$gene

# Remove the 'gene' column since it's now the row names
all_lins_p$gene <- NULL

nLins = 12
n = 12 # corresponds to the number of lineages present (i.e. no missing data)
# Run 10,000 replicate simulations of this situation and build the correlation matrix for the order statistics from them
emp_p_null_dat <- t(replicate(40000, PicMin:::GenerateNullData(1.0, n, 0.5, 3,
10000)))
# Calculate the order statistics' p-values for each simulation
emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, PicMin:::orderStatsPValues))
# Take a look at the p-values - aren't they nice?# Use those p-values to construct the correlation matrix

null_pMax_cor_unscaled <- cor(emp_p_null_dat_unscaled) 
lins_p_5 <-  as.matrix(all_lins_p[rowSums(is.na(all_lins_p)) == nLins-n,])


# Make some containers for the PicMin results
resulting_p <- rep(-1,
             nrow(lins_p_5))
resulting_n <- rep(-1,
             nrow(lins_p_5))
                     
# For each of the lines in the dataframe, perform PicMin
for (i in seq(nrow(lins_p_5)) ){
    test_result <- PicMin:::PicMin(na.omit(lins_p_5[i,]),
                                   null_pMax_cor_unscaled,                                     
                                   numReps = numReps)
    # Store the p-value
    resulting_p[i] <- test_result$p
    resulting_n[i] <- test_result$config_est
}

lins_p_5_result = data.frame(numLin = n ,
                                p = resulting_p,
                                q = p.adjust(resulting_p, method = "fdr"),
                                n_est = resulting_n,
                                locus = row.names(lins_p_5) )
picMin_results <- cbind( lins_p_5_result,
                         read.csv(text=row.names(lins_p_5), 
                                  header=FALSE, 
                                  sep = "_",
                                  col.names=c('redundan','scaffold','start')) )



# Dynamically modify the output file name based on numReps
output_file <- paste0("picMin_results_deltap", numReps, "rep.csv")

# Write the results to the output file
write.csv(picMin_results, output_file, row.names = FALSE)

# Inform the user where the file was saved
cat("Results saved to:", output_file, "\n")
