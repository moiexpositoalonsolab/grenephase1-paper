### Goal: to explore frequency trajectories in genes of possible importance
# The idea is to unpack the LFMM or others into something that is 
# visually clear

# cleaned version

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
siteofinterest=43 # Lesbos
siteofinterest=5 # Madrid
siteofinterest=c(4,5,32,43) # Cadiz
siteofinterest=46 # Wuetzburg
siteofinterest=4 # Cadiz
siteofinterest=24 #  Brixen im Thale, Austria
sres<-
  camrestop %>% dplyr::filter(site %in% siteofinterest, p==min(p))

# Get the dataset for calculations
d <- sres


#####*********************************************************************######

# Subset data for each generation
p0 <- d %>% dplyr::filter(year == 0) %>% dplyr::select(freq, rep, year, flowers)
p1 <- d %>% dplyr::filter(year == 1) %>% dplyr::select(freq, rep, year, flowers)
p2 <- d %>% dplyr::filter(year == 2) %>% dplyr::select(freq, rep, year, flowers)
p3 <- d %>% dplyr::filter(year == 3) %>% dplyr::select(freq, rep, year, flowers)

#####*********************************************************************######
##### Grid search and likelihood ######

# Function to compute the likelihood (pseudo-posterior) of s for binomial sampling
compute_likelihood_distribution <- function(p0, k, n, s_grid = seq(-1, 1, length.out = 100)) {
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


# Function to compute p(t) after t generations under selection
  predict_p_t <- function(p0, s, t) {
    p_t <- p0
    for (i in 1:t) {
      p_t <- (p_t * (1 + s)) / (p_t * (1 + s) + (1 - p_t))
    }
    return(p_t)
  }
  
# Function to compute the likelihood (pseudo-posterior) of s for binomial sampling
compute_likelihood_distribution_t <- function(pstart, k, n, s_grid = seq(-1, 1, length.out = 100), mytime=1) {

  # Initialize matrix to store likelihoods
  likelihood_matrix <- matrix(0, nrow = length(s_grid), ncol = length(pstart))
  
  i=1
  for (i in seq_along(p0)) {
    # Data for this observation
    p <- pstart[i]
    successes <- k[i]
    trials <- n[i]
    
    # Calculate likelihood for each s in the grid
    j=1 # debug
    for (j in seq_along(s_grid)) {
      s <- s_grid[j]
      # p_adjusted <- (p * (1 + s)) / (p * (1 + s) + (1 - p))
      p_adjusted <- predict_p_t(p, s, mytime) # Use t-1 steps to project to t
      
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
pseudo_posterior1 <- compute_likelihood_distribution_t(pstart=pstart, k=k1, n=flowers1 ,s_grid = seq(-1, 2, length.out = 100),mytime=1)
pseudo_posterior2 <- compute_likelihood_distribution_t(pstart, k2, flowers2, s_grid = seq(-1, 2, length.out = 100),mytime=2)
pseudo_posterior3 <- compute_likelihood_distribution_t(pstart, k3, flowers3, s_grid = seq(-1, 2, length.out = 100),mytime=3)

pseudo_posterior1

# Plot the pseudo-posterior distribution
library(ggplot2)
selectionplot<-
ggplot() +
  geom_line(data=pseudo_posterior1, aes(x = s, y = pseudo_posterior),color= "#74C476", lwd=2) +
  # geom_line(data=pseudo_posterior2, aes(x = s, y = pseudo_posterior),color= "#238B45") +
  # geom_line(data=pseudo_posterior3, aes(x = s, y = pseudo_posterior),color= "#00441B") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(-1,+2,by=0.25))+
  labs(
    y = "Probability (data | s)",
    x = "Selection coefficient (s)"
  )
selectionplot

save_plot(paste0("figs/fig-CAM5-selection-coefficient-site-",siteofinterest,".pdf"),
          selectionplot,base_height = 6,base_width = 6)
save_plot(paste0("figs/fig-CAM5-selection-coefficient-site-",siteofinterest,".png"),
          selectionplot,base_height = 6,base_width = 6)

# Get summary statistics
pseudo_posterior<-pseudo_posterior1

# Normalize the pseudo-posterior
pseudo_posterior <- pseudo_posterior %>%
  mutate(pseudo_posterior = pseudo_posterior / sum(pseudo_posterior))

# Weighted mean
weighted_mean <- weighted.mean(pseudo_posterior$s, pseudo_posterior$pseudo_posterior)

# Weighted standard deviation
weights <- pseudo_posterior$pseudo_posterior
weighted_var <- sum(weights * (pseudo_posterior$s - weighted_mean)^2) / sum(weights)
weighted_sd <- sqrt(weighted_var)

# Effective sample size
effective_sample_size <- (sum(weights)^2) / sum(weights^2)

# Weighted standard error
weighted_se <- weighted_sd / sqrt(effective_sample_size)

# 95% Confidence Interval
lower_bound <- weighted_mean - 1.96 * weighted_se 
upper_bound <- weighted_mean + 1.96 * weighted_se

sink(paste0("tables/test-selection-coef-CAM5-site-",siteofinterest,".txt"))
# Print results
cat("Weighted Mean:", weighted_mean, "\n")
cat("Weighted SE:", weighted_se, "\n")
cat("95% CI:", (lower_bound ), "to", (upper_bound ), "\n")
sink()

mean(p1$freq-p0$freq)
mean(p3$freq-p0$freq)
mean(p3$freq)
mean(p2$freq)
mean(p1$freq)
mean(p0$freq)
