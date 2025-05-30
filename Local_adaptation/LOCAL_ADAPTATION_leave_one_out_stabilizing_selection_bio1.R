rm(list=ls())

library(corrplot)
library(MCMCglmm)
library(r2glmm)
library(genio)
library(MASS)
library(grDevices)
library(dplyr)
library(patchwork)

## local adaptation

## This script run the leave-one-out analysis on all sites to predict the rest site


prefix <- "LOCAL_ADAPTATION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

## load the stabilizing selection table (era5)
mydata <- read.delim("data-intermediate/LOCAL_ADAPTATION_stabilizing_selection_data_2024Jun18_collectionsite_era5_allyears_pt1_pt0.txt",header = T) %>%
  filter(!site==33)
mydata$ecotype <- factor(mydata$ecotype,levels = unique(mydata$ecotype))
mydata$site <- factor(mydata$site,levels = unique(mydata$site))
unique_sites <- levels(mydata$site)

### First build stabilizing selection model with the climate of collection site with only Bio1 ### 

## generation 1
df1 <- mydata %>%
  filter(generation==1) %>% 
  mutate(log_p_ratio = log(pt1/pt0)) %>% 
  filter(!is.infinite(log_p_ratio) & !is.na(log_p_ratio)) %>%
  mutate(bio1_diff_sq = (bio1_ecotype - bio1_site)^2)

## the univariate stabilizing selection model with only bio1

all_posterior <- as.data.frame(matrix(nrow=493,ncol=30))

for(i in 1:30){
  
  df <- df1 %>%
    filter(site != unique_sites[i])
  
  m1.1 <- MCMCglmm(log_p_ratio ~ 1 + bio1_diff_sq, random = ~ site + us(1+bio1_diff_sq):ecotype, data = df,
                   family = "gaussian", verbose = T,nitt = 30000,burnin = 3000, pr = TRUE, saveX = TRUE, saveZ = TRUE,thin = 5)
  
  posterior_m1.1 <- apply(m1.1$Sol,2,mean)
  all_posterior[,i] <- posterior_m1.1
  print(i)
}

write.table(all_posterior,"data-intermediate/stabilizing_selection_leaveoneout_posterior_allruns_bio1_ecotypeorigin.txt",append = F,quote = F,sep = "\t",row.names = F)


## Secondly, build stabilizing selection model with the climate of collection site with Bio1 and Bio18 ### 

## generation 1
df1 <- mydata %>%
  filter(generation==1) %>% 
  mutate(log_p_ratio = log(pt1/pt0)) %>% 
  filter(!is.infinite(log_p_ratio) & !is.na(log_p_ratio)) %>%
  mutate(bio1_diff_sq = (bio1_ecotype - bio1_site)^2) %>%
  mutate(bio18_diff_sq = (bio18_ecotype / 10 - bio18_site / 10)^2) %>%
  mutate(bio1_bio18_diff = ( bio1_ecotype - bio1_site ) * (bio18_ecotype / 10 - bio18_site / 10)) 

## the univariate stabilizing selection model with bio1

all_posterior <- as.data.frame(matrix(nrow=957,ncol=30))

for(i in 1:30){
  
  df <- df1 %>%
    filter(site != unique_sites[i])
  
  m1.2 <- MCMCglmm(log_p_ratio ~ 1 + bio1_diff_sq + bio18_diff_sq + bio1_bio18_diff, random = ~ site + ecotype + us(bio1_diff_sq + bio18_diff_sq + bio1_bio18_diff):ecotype, 
                   data = df,family = "gaussian", verbose = T,nitt = 30000,burnin = 3000, pr = TRUE, saveX = TRUE, saveZ = TRUE,thin = 5)
  
  posterior_m1.2 <- apply(m1.2$Sol,2,mean)
  all_posterior[,i] <- posterior_m1.2
  print(i)
}

write.table(all_posterior,"data-intermediate/stabilizing_selection_leaveoneout_posterior_allruns_bio1_bio18_ecotypeorigin.txt",append = F,quote = F,sep = "\t",row.names = F)

















# all_vs_posterior <- as.data.frame(matrix(nrow=231,ncol=30))
# rownames(all_vs_posterior) <- levels(mydata$ecotype)
# all_wmax_posterior <- as.data.frame(matrix(nrow=231,ncol=30))
# rownames(all_wmax_posterior) <- levels(mydata$ecotype)
# all_site_posterior <- as.data.frame(matrix(nrow=30,ncol=29))
# rownames(all_site_posterior) <- levels(mydata$site)
# 

# for(i in 1:30){
#   all_vs_posterior[,i] <- -1/(all_posterior[2,i] + all_posterior[263:493,i])
#   #if Vs < 0, then we fix it with the fix effect of z1_sq
#   neg_index <- which(all_vs_posterior[,i]<0)
#   all_vs_posterior[neg_index,i] <- -1/(all_posterior[2,i])
#   all_wmax_posterior[,i] <- exp(all_posterior[33:263,i])
#   all_site_posterior[i,] <-  exp(all_posterior[3:32,i])
# }
# 
# stabilizing_selection_parameters <- as.data.frame(matrix(nrow=231,ncol=5))
# colnames(stabilizing_selection_parameters) <- c("ecotype","wmax_mean","wmax_sd","vs_mean","vs_sd")
# stabilizing_selection_parameters$ecotype <- unique(mydata$ecotype)
# stabilizing_selection_parameters$wmax_mean <- apply(all_wmax_posterior,1,mean)
# stabilizing_selection_parameters$wmax_sd <- apply(all_wmax_posterior,1,sd)
# stabilizing_selection_parameters$vs_mean <- apply(all_vs_posterior,1,mean)
# stabilizing_selection_parameters$vs_sd <- apply(all_vs_posterior,1,sd)
# 
# write.table(stabilizing_selection_parameters,"stabilizing_selection_wmax_vs.txt",sep="\t",append = F,quote = F,row.names = F)
# write.table(all_vs_posterior,"stabilizing_selection_wmax_estimation_allruns.txt",sep="\t",append = F,quote = F,row.names = F)
# write.table(all_wmax_posterior,"stabilizing_selection_vs_estimation_allruns.txt",sep="\t",append = F,quote = F,row.names = F)
# write.table(all_site_posterior,"stabilizing_selection_site_estimation_allruns.txt",sep="\t",append = F,quote = F,row.names = F)


## load the PRS-calculated stabilizing selection table
mydata <- read.delim("data-intermediate/LOCAL_ADAPTATION_stabilizing_selection_data_2024Jun18_collectionsite_BSLMM_PRS_allyears_pt1_pt0.txt",header = T) %>%
  filter(!site==33)
mydata$ecotype <- factor(mydata$ecotype,levels = unique(mydata$ecotype))
mydata$site <- factor(mydata$site,levels = unique(mydata$site))
unique_sites <- levels(mydata$site)

### First build stabilizing selection model with the climate of collection site ### 

## generation 1
df1 <- mydata %>%
  filter(generation==1) %>% 
  mutate(log_p_ratio = log(pt1/pt0)) %>% 
  filter(!is.infinite(log_p_ratio) & !is.na(log_p_ratio)) %>%
  mutate(bio1_diff_sq = (bio1_ecotype_prs - bio1_site)^2)

## the univariate stabilizing selection model with only bio1

all_posterior <- as.data.frame(matrix(nrow=493,ncol=30))

for(i in 1:30){
  
  df <- df1 %>%
    filter(site != unique_sites[i])
  
  m1.1 <- MCMCglmm(log_p_ratio ~ 1 + bio1_diff_sq, random = ~ site + us(1+bio1_diff_sq):ecotype, data = df,
                   family = "gaussian", verbose = T,nitt = 30000,burnin = 3000, pr = TRUE, saveX = TRUE, saveZ = TRUE,thin = 5)
  
  posterior_m1.1 <- apply(m1.1$Sol,2,mean)
  all_posterior[,i] <- posterior_m1.1
  print(i)
}

write.table(all_posterior,"data-intermediate/stabilizing_selection_leaveoneout_posterior_allruns_ecotypeprs.txt",append = F,quote = F,sep = "\t",row.names = F)

