### Goal: to explore frequency trajectories in genes of possible importance
# The idea is to unpack the LFMM or others into something that is 
# visually clear

################################################################################
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
library(topr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

#####*********************************************************************######
### FUNCTIONS ####


read_merged_allele_subset<-function(headfile="data-intermediate/merged_hapFIRE_allele_frequency-header.csv",
                                    subsetfile="data-intermediate/merged_hapFIRE_allele_frequency-AT5G10140-FLC.csv",
                                    snpfile="data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-AT5G10140-FLC.csv",
                                    frequencyp0file="data/average_seedmix_p0.txt",
                                    worldclimsitesdatafile="grene/data/worldclim_sitesdata.rda",
                                    sampledatafile="grene/data/samples_data.rda",
                                    addclimate=TRUE,
                                    addflowers=FALSE
){
  
  
  # Read header
  ahead<-read.csv(headfile,header = F)
  # Read SNP locations
  alocations<-read.table(snpfile,header = F)
  # Read the allele frequencies
  a<-read.csv(subsetfile,header = F)
  # Add columns needed
  colnames(a)<-ahead
  a$snp=alocations$V2
  a$chr=alocations$V1
  a$pos=alocations$V4
  
  # Pivot to add other things  
  a<-
    a %>% 
    # head %>%
    pivot_longer(
      cols = -c(snp, chr, pos),
      names_to = c("site", "year", "rep"),
      # names_pattern = "X(\\d+)_(\\d+)_(\\d+)",
      names_pattern = "(\\d+)_(\\d+)_(\\d+)",
      values_to = "freq"
    )
  
  # Merge with starting frequencies
  p0<-read.table(frequencyp0file, header=F)
  p0$snp<-paste0(p0$V1,"_",p0$V2)
  p0<-p0 %>% dplyr::rename(chr=V1,pos=V2,startfreq=V3)
  a<-
    p0 %>% 
    dplyr::select(startfreq, snp) %>% 
    merge(.,a,by='snp') 
  # add the year zero
  azero<-dplyr::filter(a, year==1) %>% 
    mutate(year=0, freq=startfreq)
  a<-rbind(a,azero)
  # Merge with climate
  if(addclimate){
    # Read climate
    load(worldclimsitesdatafile)
    a<-
      worldclim_sitesdata %>% 
      dplyr::select(site,starts_with("bio")) %>% 
      merge(.,a,by='site')
  }
  # merge with flower numbers
  if(addflowers){
    # Read sample size
    load("grene/data/samples_data.rda")
    samples_tmp<-samples_data %>% 
      dplyr::filter(year<2021) %>% 
      dplyr::mutate(year= year-2017) %>% 
      dplyr::rename(flowers=flowerscollected) %>%
      dplyr::rename(rep=plot) %>%
      dplyr::select(site, year, rep, flowers) %>% 
      group_by(site,year,rep) %>% 
      summarise(flowers=sum(flowers))
    # dplyr::select(site,plot,year, date,month,day,flowers) %>% 
    # dplyr::mutate(julian_day = yday(ymd(date)))
    a<-
      a %>% 
      merge(.,
            samples_tmp, by.x=c("site","year","rep"), 
            by.y=c("site","year","rep"), all.x=T
      )
  }
  # END
  return(a)  
}


#####*********************************************************************######
#### Get intermedaite datasets ####

# Intermediate datasets
p0<-read.table("data/average_seedmix_p0_LDpruned.txt", header=F)
p0<-read.table("data/average_seedmix_p0.txt", header=F)
p0$snp<-paste0(p0$V1,"_",p0$V2)
p0<-p0 %>% dplyr::rename(chr=V1,pos=V2,startfreq=V3)


# Site cliamte
load("grene/data/worldclim_sitesdata.rda")

# TAIR
tair<-read.csv("data-external/TAIR10_parsedgenes.csv",header=T)
head(tair)

# Annotated SNPs
load("data-intermediate/alleles-TAIR10_parsedgenes.rda")
annotated_snps %>% dim
annotated_snps %>% head

#####*********************************************************************######
##### GET GENE SUBSETS #####
###################################

#####*********************************************************************######
##### Read cam5 #####
cam<-read_merged_allele_subset(
  subsetfile = "data-intermediate/merged_hapFIRE_allele_frequency-AT2G27030-CAM5.csv",
  headfile="data-intermediate/merged_hapFIRE_allele_frequency-header.csv",
  snpfile="data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-AT2G27030-CAM5.csv",
  addflowers = T
)

head(cam)
tail(cam)
dim(cam)

# get LMs just for corroboration
camres<-
  cam %>% 
  dplyr::group_by(snp) %>% 
  dplyr::summarize(r=cor.test(freq-startfreq,bio1)$estimate,
                   p=cor.test(freq-startfreq,bio1)$p.value) %>% 
  merge(.,cam,by="snp")
# minimanhattan  
camres %>% 
  ggplot(.)+
  geom_point(aes(y=-log10(p), x=pos))+
  geom_vline( xintercept = 11531967)+
  geom_vline( xintercept = 11534358)+
  xlim(c(11531967-500,11534358+500))+
  theme_minimal()

# create top snps and neutral
camrestop<-
  camres %>% 
  dplyr::filter(-log10(p)>14) %>% 
  dplyr::filter(startfreq<0.5) 

#### estimate selection in site 4 
siteofinterest=46 # Wuetzburg
siteofinterest=4 # Cadiz
siteofinterest=4 # Cadiz
sres<-
  camrestop %>% dplyr::filter(site==siteofinterest, p==min(p))

# Get the dataset for calculations
d <- sres

# Subset data for each generation
p0 <- d %>% dplyr::filter(year == 0) %>% dplyr::select(freq, rep, year, flowers)
p1 <- d %>% dplyr::filter(year == 1) %>% dplyr::select(freq, rep, year, flowers)
p2 <- d %>% dplyr::filter(year == 2) %>% dplyr::select(freq, rep, year, flowers)
p3 <- d %>% dplyr::filter(year == 3) %>% dplyr::select(freq, rep, year, flowers)

#####*********************************************************************######
##### Grid search and likelihood ######
# 
# 
# # Subset dataset
# d <- sres
# 
# # Selection coefficient to evaluate
# s_values <- seq(0, 0.1, length.out = 100) # Range of s values
# s_values <- seq(0, 10, length.out = 100) # Range of s values
# s_values <- seq(-10, 10, length.out = 100) # Range of s values
# iterations <- 10 # Number of sampling iterations per s value
# 
# # Subset data for each generation
# p0 <- d %>% dplyr::filter(year == 0) %>% dplyr::select(freq, rep, year, flowers)
# p1 <- d %>% dplyr::filter(year == 1) %>% dplyr::select(freq, rep, year, flowers)
# p2 <- d %>% dplyr::filter(year == 2) %>% dplyr::select(freq, rep, year, flowers)
# p3 <- d %>% dplyr::filter(year == 3) %>% dplyr::select(freq, rep, year, flowers)
# 
# # Function to compute log-likelihood across multiple generations
# compute_likelihood_across_generations <- function(pstart, generations, flowers, scoef) {
#   log_likelihood <- 0
#   
#   # Iterate through each generation
#   for (gen in seq_along(generations)) {
#     # p_adjusted <- p0 + s
#     p_adjusted <- (pstart * (1 + scoef)) / (pstart * (1 + scoef) + (1 - pstart)) # THIS IS THE POPGEN DEFINITION
#     
#     if (any(p_adjusted <= 0 | p_adjusted >= 1)) {
#       message("Adjusted probabilities (p0 + s) must be between 0 and 1. re-scaling")
#       p_adjusted<-ifelse(p_adjusted <= 0, 0,p_adjusted)
#       p_adjusted<-ifelse(p_adjusted >=1, 1, p_adjusted)
#     }
#     # Convert observed frequencies to counts
#     k <- round(generations[[gen]] * flowers[[gen]])
#     
#     # Compute log-likelihood for this generation
#     # log_likelihood <- log_likelihood + 
#     #   sum(k * log(p_adjusted) + (flowers[[gen]] - k) * log(1 - p_adjusted))
#     # message("k=", k,", flowers=",flowers, ", p_adj=",p_adjusted)
#     log_likelihood <- sum(dbinom(k, size = flowers[[gen]], prob = p_adjusted, log = TRUE))
#     
#     # Update p0 for the next generation
#     pstart <- generations[[gen]]
#   }
#   
#   return(log_likelihood)
# }
# 
# # Perform sampling and likelihood calculations for each s value
# results <- data.frame(
#   s = numeric(0),
#   iteration = numeric(0),
#   log_likelihood = numeric(0)
# )
# 
# # Iteratively compute likelihoods
# for (s in s_values) {
#   message("Likelihood with s ", s)
#   for (i in 1:iterations) {
#     message("Iteration ", i)
#     
#     # Compute log-likelihood across generations
#     # ll <- compute_likelihood_across_generations(
#     #   pstart=p0$freq, 
#     #   generations=list(p1$freq, p2$freq, p3$freq), 
#     #   flowers=list(p1$flowers, p2$flowers, p3$flowers), 
#     #   scoef=s
#     # )
#     # Debug with one generation
#     ll <- compute_likelihood_across_generations(
#       pstart=p0$freq, 
#       generations=list(p1$freq), 
#       flowers=list(p1$flowers), 
#       scoef=s
#     )
#     
#     # Store results
#     results <- rbind(results, 
#                      data.frame(s = s, iteration = i, log_likelihood = ll))
#   }
# }
# 
# # Summarize results to compute posterior probabilities
# posterior <- results %>%
#   group_by(s) %>%
#   summarize(mean_log_likelihood = mean(log_likelihood), .groups = "drop") %>%
#   mutate(posterior_prob = exp(mean_log_likelihood - max(mean_log_likelihood))) %>%
#   mutate(posterior_prob = posterior_prob / sum(posterior_prob)) # Normalize
# 
# # Plot posterior probabilities
# selection_plot<-
# ggplot(posterior, aes(x = s, y = posterior_prob)) +
#   geom_line() +
#   theme_minimal() +
#   labs(
#       # title = "Posterior Probability of Selection Coefficient (s)",
#        x = "Selection coefficient (s)",
#        y = "Posterior probability")+
#   geom_vline(xintercept = 0,lty="dotted")
# selection_plot

# ggsave(file = paste0("figs/fig-CAM5-selection-coefficient-site",siteofinterest,".pdf"),
#        height=4, width = 5,
#        selection_plot
#         )
# ggsave(file = paste0("figs/fig-CAM5-selection-coefficient-site",siteofinterest,".png"),
#        height=4, width = 5,
#        selection_plot
# )

# # Function to calculate posterior mean and 95% HPDI
# calculate_posterior_summary <- function(data) {
#   # Filter out rows with -Inf log-likelihoods
#   valid_data <- data %>%
#     filter(log_likelihood > -Inf)
#   
#   # Check if there are any valid rows left
#   if (nrow(valid_data) == 0) {
#     stop("All log-likelihood values are -Inf. Cannot compute posterior summary.")
#   }
#   
#   # Calculate weights from log-likelihoods
#   weights <- exp(valid_data$log_likelihood - max(valid_data$log_likelihood))
#   
#   # Normalize weights
#   weights <- weights / sum(weights)
#   
#   # Calculate posterior mean
#   posterior_mean <- sum(valid_data$s * weights)
#   
#   # Calculate 95% HPDI
#   sorted_data <- valid_data %>%
#     arrange(s)
#   cum_weights <- cumsum(weights[order(valid_data$s)])
#   
#   # Find indices for the 95% interval
#   lower_index <- which.min(abs(cum_weights - 0.025))
#   upper_index <- which.min(abs(cum_weights - 0.975))
#   
#   hpdi <- c(sorted_data$s[lower_index], sorted_data$s[upper_index])
#   
#   # Return results
#   list(
#     posterior_mean = posterior_mean,
#     hpdi = hpdi
#   )
# }
# calculate_posterior_summary(results)
 
# 
# p0=0.1
# s=0.5
# (p0 * (1 + s)) / (p0 * (1 + s) + (1 - p0))


#####*********************************************************************######
#####* New likelihood approach ####

# Function to compute the likelihood (pseudo-posterior) of s for binomial sampling
compute_likelihood_distribution <- function(p0, k, n, s_grid = seq(-0.05, 1, length.out = 100)) {
  # Initialize matrix to store likelihoods
  likelihood_matrix <- matrix(0, nrow = length(s_grid), ncol = length(p0))
  
  for (i in seq_along(p0)) {
    # Data for this observation
    p <- p0[i]
    successes <- k[i]
    trials <- n[i]
    
    # Calculate likelihood for each s in the grid
    for (j in seq_along(s_grid)) {
      s <- s_grid[j]
      p_adjusted <- (p * (1 + s)) / (p * (1 + s) + (1 - p))
      
      if (any(p_adjusted <= 0 | p_adjusted >= 1)) {
        message("Adjusted probabilities (p1) must be between 0 and 1. re-scaling")
        p_adjusted<-ifelse(p_adjusted <= 0, 0,p_adjusted)
        p_adjusted<-ifelse(p_adjusted >=1, 1, p_adjusted)
      }
      
      # Check for valid probabilities
      if (p_adjusted > 0 && p_adjusted < 1) {
        # Compute binomial likelihood
        likelihood_matrix[j, i] <- dbinom(successes, size = trials, prob = p_adjusted, log = TRUE)
      } else {
        likelihood_matrix[j, i] <- -Inf # Invalid s produces -Inf likelihood
      }
    }
  }
  
  # Combine likelihoods across observations by summing log-likelihoods
  total_log_likelihood <- rowSums(likelihood_matrix)
  
  # Convert to pseudo-posterior by exponentiating and normalizing
  pseudo_posterior <- exp(total_log_likelihood - max(total_log_likelihood))
  pseudo_posterior <- pseudo_posterior / sum(pseudo_posterior)
  
  return(data.frame(s = s_grid, pseudo_posterior = pseudo_posterior))
}


# Example data for three generations
pstart <- p0$freq

# Observed frequencies for generations
p1_freq <- p1$freq
p2_freq <- p2$freq
p3_freq <- p3$freq


# Observed counts for each generation
k1 <- round(p1_freq * p1$flowers)
k2 <- round(p2_freq * p2$flowers)
k3 <- round(p3_freq * p3$flowers)

# Flowers
flowers1<-p1$flowers
flowers2<-p2$flowers
flowers3<-p3$flowers

# Combine data into lists
generations <- list(k1, k2, k3)
flowers_list <- list(p1$flowers, p2$flowers, p3$flowers)

# Compute pseudo-posterior distribution
pseudo_posterior1 <- compute_likelihood_distribution(pstart, k1, flowers1)
pseudo_posterior2 <- compute_likelihood_distribution(pstart, k2, flowers2)
pseudo_posterior3 <- compute_likelihood_distribution(pstart, k3, flowers3)

# Plot the pseudo-posterior distribution
library(ggplot2)
ggplot() +
  geom_line(data=pseudo_posterior1, aes(x = s, y = pseudo_posterior)) +
  geom_line(data=pseudo_posterior2, aes(x = s, y = pseudo_posterior)) +
  geom_line(data=pseudo_posterior3, aes(x = s, y = pseudo_posterior)) +
  theme_minimal() +
  labs(
       y = "Probability of selection coefficient (s)",
       x = "Selection coefficient (s)"
       )
