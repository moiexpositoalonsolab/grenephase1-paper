.libPaths("/home/tbellagio/miniforge3/envs/r-environment")

packageVersion("fastmap")

library(devtools) # Needed to install the PicMin package from Tom's GitHub install_github("TBooker/PicMin", force = T)

.libPaths("/home/tbellagio/miniforge3/envs/r-environment/lib/R/library")

library(PicMin)
library(tidyverse)
library(poolr)

snp_names <- read.csv("../key_files/var_pos_grenenet.csv", header = TRUE, sep = ",")

snp_names <- snp_names[snp_names$total_alleles05filter != "", ]

lfmm_bio1 <- read.csv("lfmm/lfmm_bio1_kk20_results.csv", header = TRUE, sep = ",")

lfmm_bio1 <- data.frame(name = snp_names$id, lfmm = lfmm_bio1$p_value)

lmm_bio1 <- read.csv("lmm/results_lmm_bio1.csv", header = TRUE, sep = ",")

lmm_bio1 <- data.frame(name = snp_names$id, lmm = lmm_bio1$beta_p)

kendall_bio1 <- read.csv("kendall_tau/kendall_corr_bio1.csv", header = TRUE, sep = ",")

kendall_bio1 <- data.frame(name = snp_names$id, kendall = kendall_bio1$K_tau_p)

#put all data frames into list 
df_list <- list(lfmm_bio1,
                lmm_bio1, 
               kendall_bio1)

#merge all data frames in list - use the 'window' variable to merge 
all_lins <- df_list %>% reduce(full_join, by='name')

# remove the column named "window"
all_lins_p <- all_lins[ , !(names(all_lins) %in% c("name"))]

# Use the "window" column as row.names 
rownames(all_lins_p) <- all_lins$name

nLins = 3
n = 3 # corresponds to the number of lineages present (i.e. no missing data)
# Run 10,000 replicate simulations of this situation and build the correlation matrix for the order statistics from them
emp_p_null_dat <- t(replicate(40000, PicMin:::GenerateNullData(1.0, n, 0.5, 3,
10000)))
# Calculate the order statistics' p-values for each simulation
emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, PicMin:::orderStatsPValues))
# Take a look at the p-values - aren't they nice?# Use those p-values to construct the correlation matrix

## becasue The emp_p_null_dat_unscaled is created using the orderStatsPValues() function, which calculates order statistics. 
# When calculating order statistics, you typically remove the highest value (since it's redundant for analysis). Therefore, for 5 species, this function outputs the 4 smallest p-values (excluding the largest one).
## teh nwe cannot run this with 2 sets of pvalues 

null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled ) 

# Select the loci that have daa for exactly 7 lineages
lins_p_5 <-  as.matrix(all_lins_p[rowSums(is.na(all_lins_p)) == nLins-n,])

# Make some containers for the PicMin results
resulting_p <- rep(-1,
             nrow(lins_p_5))
resulting_n <- rep(-1,
             nrow(lins_p_5))
numReps = 100000 ## 100x larger than before

print('about to un big run')
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

write.csv(picMin_results, "picMin_results_lfmm_lmm_kendalltau_100000rep.csv", row.names = FALSE)