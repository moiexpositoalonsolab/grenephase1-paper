
rm(list=ls())

library(ggplot2)
library(dplyr)
library(patchwork)
library(SNPRelate)
library(SeqArray)
library(vcfR)

## POPULATION EVOLUTION

## This script will analyze the real experimental delta frequency (snp and ecotype) across 3 generations
## It will calculate delta frequency in two ways: 1) all compare to p0 p1-p0, p2-p0 2) per generation change p1-p0, p2-p1


prefix <- "POP_EVOLUTION"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

#####  load meta information, SNP and ecotype frequencies ######

meta <- read.delim("data/merged_sample_table.csv",sep=",",stringsAsFactors = F)
terminal_index <- read.delim("data/site_plot_terminal_index.txt")

## SNP frequencies
ld_p0 <- read.delim("data/average_seedmix_p0_LDpruned.txt",header=F)
ld_p <- read.delim("data/merged_hapFIRE_allele_frequency_LDpruned.txt",header=T,check.names = F)

## ecotype frequencies
ecotype_p0 <- read.delim("data/founder_ecotype_frequency.txt",header=F)
ecotype_p <- read.delim("data/merged_ecotype_frequency.txt",header=T,check.names = F)


##### calculate the per generation delta frequency: ALL compare to p0 
#### we are using plots that have all 3 generations of data

## generation 3
index_gen_3 <- which(terminal_index$generation1==1 &terminal_index$generation2==2 & terminal_index$generation3==3)


## generation 1 frequency change
gen_1_snp_deltap_avg <- c()
gen_1_snp_deltap_var <- c()
gen_1_ecotype_deltap_avg <- c()
for(i in index_gen_3){
  snp_p0 <- ld_p0$V3
  snp_p1 <- ld_p[,terminal_index$index1[i]]
  eco_p0 <- ecotype_p0$V2
  eco_p1 <- ecotype_p[,terminal_index$index1[i]]
  gen_1_snp_deltap_avg <- c(gen_1_snp_deltap_avg,mean(snp_p1-snp_p0))
  gen_1_snp_deltap_var <- c(gen_1_snp_deltap_var,var(snp_p1-snp_p0))
  gen_1_ecotype_deltap_avg <- c(gen_1_ecotype_deltap_avg,mean(eco_p1 - eco_p0 ))
}
mean(gen_1_snp_deltap_avg)
mean(gen_1_snp_deltap_var)
mean(gen_1_ecotype_deltap_avg)


## generation 2 frequency change
gen_2_snp_deltap_avg <- c()
gen_2_snp_deltap_var <- c()
gen_2_ecotype_deltap_avg <- c()
for(i in index_gen_3){
  snp_p0 <- ld_p0$V3
  snp_p1 <- ld_p[,terminal_index$index2[i]]
  eco_p0 <- ecotype_p0$V2
  eco_p1 <- ecotype_p[,terminal_index$index2[i]]
  gen_2_snp_deltap_avg <- c(gen_2_snp_deltap_avg,mean(snp_p1-snp_p0))
  gen_2_snp_deltap_var <- c(gen_2_snp_deltap_var,var(snp_p1-snp_p0))
  gen_2_ecotype_deltap_avg <- c(gen_2_ecotype_deltap_avg,mean(eco_p1 - eco_p0 ))
}
mean(gen_2_snp_deltap_avg)
mean(gen_2_snp_deltap_var)
mean(gen_2_ecotype_deltap_avg)



## generation 3 frequency change
gen_3_snp_deltap_avg <- c()
gen_3_snp_deltap_var <- c()
gen_3_ecotype_deltap_avg <- c()
for(i in index_gen_3){
  snp_p0 <- ld_p0$V3
  snp_p1 <- ld_p[,terminal_index$index3[i]]
  eco_p0 <- ecotype_p0$V2
  eco_p1 <- ecotype_p[,terminal_index$index3[i]]
  gen_3_snp_deltap_avg <- c(gen_3_snp_deltap_avg,mean(snp_p1-snp_p0))
  gen_3_snp_deltap_var <- c(gen_3_snp_deltap_var,var(snp_p1-snp_p0))
  gen_3_ecotype_deltap_avg <- c(gen_3_ecotype_deltap_avg,mean(eco_p1 - eco_p0 ))
}
mean(gen_3_snp_deltap_avg)
mean(gen_3_snp_deltap_var)
mean(gen_3_ecotype_deltap_avg)

snp_delta_1 <- as.data.frame(matrix(nrow=3*142,ncol=2))
colnames(snp_delta_1) <- c("delta_frequency","generation")
ecotype_delta <- as.data.frame(matrix(nrow=3*142,ncol=2))
colnames(ecotype_delta) <- c("delta_frequency","generation")

snp_delta_1$delta_frequency <- c(gen_1_snp_deltap_avg,gen_2_snp_deltap_avg,gen_3_snp_deltap_avg)
snp_delta_1$generation <- c(rep(1,142),rep(2,142),rep(3,142))

ecotype_delta$delta_frequency <- c(gen_1_ecotype_deltap_avg,gen_2_ecotype_deltap_avg,gen_3_ecotype_deltap_avg)
ecotype_delta$generation <- c(rep(1,142),rep(2,142),rep(3,142))

means_1 <- snp_delta_1 %>%
  group_by(generation) %>%
  dplyr::summarize(mean_delta = mean(delta_frequency))

p1 <- ggplot(snp_delta_1,aes(x=delta_frequency))+
  geom_histogram(alpha = 0.9,position = 'identity',bins = 30,fill="grey40")+
  geom_density(alpha=0.8,fill="grey80")+
  facet_wrap(~generation)+
  geom_vline(xintercept = 0.0009,linetype="dashed", color = "black")+
  geom_text(data = means_1, aes(x = -0.03, y = 80, label = paste("Mean =", round(mean_delta, 5))), vjust = -1) +
  ylab("Frequency")+
  xlim(c(-0.05,0.05))+
  ylim(c(0,100))+
  theme_minimal()

p1 <- ggplot(snp_delta_1,aes(x=delta_frequency))+
  geom_histogram(alpha = 0.9,position = 'identity',bins = 30,fill="grey40")+
  geom_density(alpha=0.8,fill="grey80")+
  facet_wrap(~generation)+
  geom_vline(xintercept = 0.0009,linetype="dashed", color = "black")+
  geom_text(data = means_1, aes(x = -0.03, y = 80, label = paste("Mean =", round(mean_delta, 5))), vjust = -1) +
  ylab("Frequency")+
  xlim(c(-0.05,0.05))+
  ylim(c(0,100))+
  theme_minimal()


p1 <- ggplot(snp_delta_1,aes(x=delta_frequency))+
  #geom_histogram(alpha = 0.9,position = 'identity',bins = 30,fill="grey40")+
  geom_density(alpha=0.8,fill="grey80")+
  facet_wrap(~generation)+
  geom_vline(xintercept = 0.0009,linetype="dashed", color = "black")+
  geom_text(data = means_1, aes(x = -0.03, y = 80, label = paste("Mean =", round(mean_delta, 5))), vjust = -1) +
  ylab("Frequency")+
  xlim(c(-0.05,0.05))+
  ylim(c(0,100))+
  theme_minimal()


p1

##### calculate the per generation delta frequency: per generation change p1-p0, p2-p1 

# 
# logit <- function(p){return(log(p/(1-p)))}
# arcsine <- function(p){return(asin(sqrt(p)))}
# normalize <- function(p1,p0) {return((p1-p0)/p0/(1-p0))}

## generation 1 frequency change
gen_1_snp_deltap_avg <- c()
gen_1_ecotype_deltap_avg <- c()
for(i in index_gen_1){
  snp_p0 <- ld_p0$V3
  snp_p1 <- ld_p[,terminal_index$index1[i]]
  eco_p0 <- ecotype_p0$V2
  eco_p1 <- ecotype_p[,terminal_index$index1[i]]
  #gen_1_snp_deltap_avg <- c(gen_1_snp_deltap_avg,mean(arcsine(snp_p1)-arcsine(snp_p0)))
  gen_1_snp_deltap_avg <- c(gen_1_snp_deltap_avg,mean(snp_p1-snp_p0))
  gen_1_ecotype_deltap_avg <- c(gen_1_ecotype_deltap_avg,mean(eco_p1 - eco_p0 ))
}
mean(gen_1_snp_deltap_avg)
mean(gen_1_ecotype_deltap_avg)

hist(gen_1_snp_deltap_avg)
## generation 2 frequency change
gen_2_snp_deltap_avg <- c()
gen_2_ecotype_deltap_avg <- c()
for(i in index_gen_2){
  snp_p0 <- ld_p[,terminal_index$index1[i]]
  snp_p1 <- ld_p[,terminal_index$index2[i]]
  eco_p0 <- ecotype_p[,terminal_index$index1[i]]
  eco_p1 <- ecotype_p[,terminal_index$index2[i]]
  gen_2_snp_deltap_avg <- c(gen_2_snp_deltap_avg,mean(snp_p1-snp_p0))
  #gen_2_snp_deltap_avg <- c(gen_2_snp_deltap_avg,mean(arcsine(snp_p1)-arcsine(snp_p0)))
  gen_2_ecotype_deltap_avg <- c(gen_2_ecotype_deltap_avg,mean(eco_p1 - eco_p0 ))
}
mean(gen_2_snp_deltap_avg)
mean(gen_2_ecotype_deltap_avg)
hist(gen_2_snp_deltap_avg)



## generation 3 frequency change
gen_3_snp_deltap_avg <- c()
gen_3_ecotype_deltap_avg <- c()
for(i in index_gen_3){
  snp_p0 <- ld_p[,terminal_index$index2[i]]
  snp_p1 <- ld_p[,terminal_index$index3[i]]
  eco_p0 <- ecotype_p[,terminal_index$index2[i]]
  eco_p1 <- ecotype_p[,terminal_index$index3[i]]
  gen_3_snp_deltap_avg <- c(gen_3_snp_deltap_avg,mean(snp_p1-snp_p0))
  #gen_3_snp_deltap_avg <- c(gen_3_snp_deltap_avg,mean(arcsine(snp_p1)-arcsine(snp_p0)))
  gen_3_ecotype_deltap_avg <- c(gen_3_ecotype_deltap_avg,mean(eco_p1 - eco_p0 ))
}
mean(gen_3_snp_deltap_avg)
mean(gen_3_ecotype_deltap_avg)
hist(gen_3_snp_deltap_avg)



snp_delta_2 <- as.data.frame(matrix(nrow=326+203+142,ncol=2))
colnames(snp_delta_2) <- c("delta_frequency","generation")

snp_delta_2$delta_frequency <- c(gen_1_snp_deltap_avg,gen_2_snp_deltap_avg,gen_3_snp_deltap_avg)
snp_delta_2$generation <- c(rep(1,326),rep(2,203),rep(3,142))

ecotype_delta$delta_frequency <- c(gen_1_ecotype_deltap_avg,gen_2_ecotype_deltap_avg,gen_3_ecotype_deltap_avg)
ecotype_delta$generation <- c(rep(1,326),rep(2,203),rep(3,142))

means_2 <- snp_delta_2 %>%
  group_by(generation) %>%
  dplyr::summarize(mean_delta = mean(delta_frequency))


p2 <- ggplot(snp_delta_2,aes(x=delta_frequency))+
  #geom_histogram(alpha = 0.9,position = 'identity',bins = 30,fill="grey40")+
  geom_density(alpha=0.8,fill="grey80")+
  facet_wrap(~generation)+
  geom_vline(xintercept = 0,linetype="dashed", color = "black")+
  geom_text(data = means_2, aes(x = -0.03, y = 60, label = paste("Mean =", round(mean_delta, 5))), vjust = -1) +
  ylab("Density")+
  xlim(c(-0.05,0.05))+
 # ylim(c(0,100))+
  theme_minimal()


p2
combined <- p1/p2

ggsave(plot = combined,filename =  paste0("figs/",prefix,"real_grene_net_snp_delta_frequency.pdf"),device = "pdf",width = 10,height = 8)

ggsave(plot = p2,filename =  paste0("figs/",prefix,"real_grene_net_snp_delta_frequency_per_generation_density.pdf"),device = "pdf",width = 8,height = 2)


## Now we are going to simulate SNP frequencies for each generation using the previous generation frequency 
## and known sampled flowers

genofile <- seqOpen("data-intermediate/greneNet_final_v1.1_LDpruned.gds")
# get the dosage of reference allele
gt <- t(seqGetData(genofile, "$dosage"))

dim(gt)
seqClose(genofile)

## simulation based on random sampling from the founder population
p_drift_1 <- as.data.frame(matrix(ncol=1000,nrow=13985))
p_drift_2 <- as.data.frame(matrix(ncol=1000,nrow=13985))
p_drift_3 <- as.data.frame(matrix(ncol=1000,nrow=13985))
N1 <-  sample(x=meta$total_flower_counts[meta$generation==1],size=1000,replace = T)
N2 <-  sample(x=meta$total_flower_counts[meta$generation==2],size=1000,replace = T)
N3 <-  sample(x=meta$total_flower_counts[meta$generation==3],size=1000,replace = T)

for(i in 1:1000){
  ## gen 1 sampling
  index_1 <- sample(ncol(gt),N1[i],replace = T,prob = ecotype_p0$V2)
  gt_ <- gt[,index_1]
  if(N1[i]==1){
    p1 <- 1 - gt_/2
  }else{
    p1 <- apply(gt_,1,function(x) 1 - sum(x)/N1[i]/2)
  }
  p_drift_1[,i] <- p1
  
  ## gen 2 sampling
  ecotype_p1 <- apply(ecotype_p[,which(meta$generation==1)],1,mean)
  index_2 <- sample(ncol(gt),N2[i],replace = T,prob = ecotype_p1)
  gt_ <- gt[,index_2]
  if(N2[i]==1){
    p2 <- 1 - gt_/2
  }else{
    p2 <- apply(gt_,1,function(x) 1 - sum(x)/N2[i]/2)
  }
  p_drift_2[,i] <- p2
  
  
  ## gen 3 sampling
  ecotype_p2 <- apply(ecotype_p[,which(meta$generation==2)],1,mean)
  
  index_3 <- sample(ncol(gt),N3[i],replace = T,prob = ecotype_p2)
  gt_ <- gt[,index_3]
  if(N3[i]==1){
    p3 <- 1 - gt_/2
  }else{
    p3 <- apply(gt_,1,function(x) 1 - sum(x)/N3[i]/2)
  }
  p_drift_3[,i] <- p3
  
  print(i)
}


## generation 1 frequency change
gen_1_snp_deltap_avg_sim <- c()
for(i in 1:1000){
  snp_p0 <- ld_p0$V3
  snp_p1 <- p_drift_1[,i]
  gen_1_snp_deltap_avg_sim <- c(gen_1_snp_deltap_avg_sim,mean(snp_p1-snp_p0))
}
mean(gen_1_snp_deltap_avg_sim)

## generation 2 frequency change
gen_2_snp_deltap_avg_sim <- c()
for(i in 1:1000){
  snp_p0 <- p_drift_1[,i]
  snp_p1 <- p_drift_2[,i]
  gen_2_snp_deltap_avg_sim <- c(gen_2_snp_deltap_avg_sim,mean(snp_p1-snp_p0))
}
mean(gen_2_snp_deltap_avg_sim)

## generation 3 frequency change
gen_3_snp_deltap_avg_sim <- c()
for(i in 1:1000){
  snp_p0 <- p_drift_2[,i]
  snp_p1 <- p_drift_3[,i]
  gen_3_snp_deltap_avg_sim <- c(gen_3_snp_deltap_avg_sim,mean(snp_p1-snp_p0))
}
mean(gen_3_snp_deltap_avg_sim)

ks.test(x=gen_3_snp_deltap_avg_sim,y = gen_3_snp_deltap_avg)
#p-value = 0.3581
ks.test(x=gen_2_snp_deltap_avg_sim,y = gen_2_snp_deltap_avg)
#p-value = 0.008316
ks.test(y=gen_1_snp_deltap_avg_sim,x = gen_1_snp_deltap_avg)
#p-value = 0.01048


snp_delta_sim <- as.data.frame(matrix(nrow=1000*3,ncol=2))
colnames(snp_delta_sim) <- c("delta_frequency","generation")

snp_delta_sim$delta_frequency <- c(gen_1_snp_deltap_avg_sim,gen_2_snp_deltap_avg_sim,gen_3_snp_deltap_avg_sim)

snp_delta_sim$generation <- c(rep(1,1000),rep(2,1000),rep(3,1000))

means_sim <- snp_delta_sim %>%
  group_by(generation) %>%
  dplyr::summarize(mean_value = mean(delta_frequency))

p_sim <- ggplot(snp_delta_sim,aes(x=delta_frequency))+
  #geom_histogram(alpha = 0.9,position = 'identity',bins = 50,fill="grey40")+
  geom_density(alpha=0.8,fill="grey80")+
  geom_vline(data = means_sim,aes(xintercept = mean_value),linetype="dashed", color = "black")+
  geom_text(data = means_sim, aes(x = -0.05, y = 100, label = paste("Mean =", round(mean_value, 5))), vjust = -1) +
  facet_wrap(~ generation) +
  ylab("Frequency")+
  theme_minimal()

p_sim
ggsave(plot = p_sim,filename =  paste0("figs/",prefix,"simulated_grene_net_snp_delta_frequency_per_generation.pdf"),device = "pdf",width = 8,height = 4)

