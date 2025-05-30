rm(list=ls())

library(ggplot2)
library(boot)
library(dplyr)

## local adaptation

## This script will use the leave-one-out results to estimate parameters (V, Wmax, and Wavg) on Bio1 only
## It will use bootstrap to estimate the mean values of different parameters
## We will use the the climate of collection for ecotype climate

prefix <- "LOCAL_ADAPTATION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

## load the dataframe
mydata <- read.delim("data-intermediate/LOCAL_ADAPTATION_stabilizing_selection_data_2024Jun18_collectionsite_era5_allyears_pt1_pt0.txt",header = T)
mydata$ecotype <- factor(mydata$ecotype,levels = unique(mydata$ecotype))
mydata$site <- factor(mydata$site,levels = unique(mydata$site))
unique_sites <- levels(mydata$site)

## load stabilizing selection model posteriors (bio1 and bio18 model)

posterior <- read.delim("data-intermediate/stabilizing_selection_leaveoneout_posterior_allruns_bio1_ecotypeorigin.txt")

wmax_pos <- matrix(nrow=231,ncol=30)
wavg_pos <- matrix(nrow=30,ncol=30)
v_bio1 <- matrix(nrow=231,ncol=30)

for(i in 1:ncol(posterior)){
  wmax_pos[,i] <- posterior[32:262,i]
  v_bio1[,i] <- posterior[2,i] + posterior[263:493,i]
  wavg_pos[,i] <- append(posterior[3:31,i],NA,after = i-1)
}

## use bootstrap to estimate the mean/median values of the parameters 
samplemean <- function(data,i){return(mean(data[i],na.rm=T))}


bio1_fix <- mean(boot(as.numeric(posterior[2,]),samplemean,R=1000)$t)



wmax_mean <- c()
for(i in 1:231){
  b <- boot(wmax_pos[i,],samplemean,R=1000)
  wmax_mean <- c(wmax_mean, mean(b$t))
}

wavg_mean <- c()
for(i in 1:30){
  b <- boot(wavg_pos[i,],samplemean,R=1000)
  wavg_mean <- c(wavg_mean, mean(b$t))
}

v_bio1_mean <- c()
for(i in 1:231){
  b <- boot(v_bio1[i,],samplemean,R=1000)
  v_bio1_mean <- c(v_bio1_mean, mean(b$t))
}

## now estimate the true Wmax, Wavg and V values (in the exponential scale)

Wmax <- exp(wmax_mean)
Wavg <- exp(-wavg_mean)
Vs <- - 1 / v_bio1_mean
index <- which(Vs<0)
table(Vs<0)
Vs[Vs<0] <- -1 / bio1_fix

site_wavg <- tibble(site = unique_sites[which(unique_sites!=33)],
                    wavg = Wavg)
write.table(site_wavg,paste0("data-intermediate/",prefix,"site_wavg.txt"),append = F,quote = F,sep = "\t",col.names = T,row.names = F)

## read ecotypes 

ecotypes <- read.delim("data/founder_ecotype_frequency.txt",header=F)

stabilizing_parameters <- as.data.frame(matrix(ncol=3,nrow=231))
colnames(stabilizing_parameters) <- c("ecotype","Wmax","Vs")
stabilizing_parameters$ecotype <- ecotypes$V1
stabilizing_parameters$Wmax <- Wmax
stabilizing_parameters$Vs <- Vs

write.table(stabilizing_parameters,paste0("data-intermediate/",prefix,"ecotype_Wmax_Vs.txt"),append = F,quote = F,sep = "\t",col.names = T,row.names = F)

cor.test(Wmax,Vs,method = "spearman")

ggplot(stabilizing_parameters,aes(x=Wmax,y=Vs))+
  geom_point(size=4,pch=16)+
  annotate(geom = "text",label= "Spearman Rho -0.3433",x = 2,y=3000,size=5)+
  theme_minimal()
ggsave(file = paste0("figs/",prefix,"Wmax_Vs_relationship"),device = "pdf",width = 5,height = 5)


## investigate the random site intercept (random effect)

generation_1_parallelism <- read.delim("data-intermediate/generation_1_parallelism.txt") %>%
  filter(source =="snp") %>%
  mutate(site_effect = Wavg)
unique_sites <- unique(generation_1_parallelism$site)


cor.test(generation_1_parallelism$mean,generation_1_parallelism$site_effect,method = "spearman")

ggplot(generation_1_parallelism,aes(x=mean,y=site_effect))+
  geom_point(cex=4)+
  annotate("text",x=0.5,y=5,label = "spearman's correlation rho: -0.7463849",size=5)+
  xlab("SNP parallelism (r)") + 
  ylab("random experiemtnal site effect")+
  theme_minimal()
ggsave(file = paste0("figs/",prefix,"site_effect_snp_parallelism.pdf"),device = "pdf",width = 5,height = 5)


