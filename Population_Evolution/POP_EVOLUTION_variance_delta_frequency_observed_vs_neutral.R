rm(list=ls())

library(ggplot2)
library(dplyr)
library(patchwork)
library(SNPRelate)
library(SeqArray)
library(vcfR)

## POPULATION EVOLUTION

## This script will analyze the real experimental delta frequency (snp and ecotype) across 3 generations
## and compare it with the neutral evolution and see if it is significantly different


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



### now the neutral simulation

## The real greneNET data to use random sampling to estimate the null distribution


genofile <- seqOpen("data-intermediate/greneNet_final_v1.1_LDpruned.gds")
# get the dosage of reference allele
gt <- t(seqGetData(genofile, "$dosage"))
dim(gt)
seqClose(genofile)
#write.table(gt,"genotype_reference_dosage.txt",col.names = F,row.names = F,append = F,quote = F,sep = "\t")
founder_ecotype_frequency <- read.delim("data/founder_ecotype_frequency.txt",header=F)

## simulation based on random sampling from the founder population
p_drift_1 <- as.data.frame(matrix(ncol=142,nrow=13985))
p_drift_2 <- as.data.frame(matrix(ncol=142,nrow=13985))
p_drift_3 <- as.data.frame(matrix(ncol=142,nrow=13985))

N1 <- c()
N2 <- c()
N3 <- c()
for(i in index_gen_3){
  N1 <- c(N1,meta$total_flower_counts[terminal_index$index1[i]])
  N2 <- c(N1,meta$total_flower_counts[terminal_index$index2[i]])
  N3 <- c(N1,meta$total_flower_counts[terminal_index$index3[i]])
}


for(i in 1:142){
  ## gen 1 sampling
  index_1 <- sample(ncol(gt),N1[i],replace = T,prob = founder_ecotype_frequency$V2)
  gt_ <- gt[,index_1]
  if(N1[i]==1){
    p1 <- 1 - gt_/2
  }else{
    p1 <- apply(gt_,1,function(x) 1 - sum(x)/N1[i]/2)
  }
  p_drift_1[,i] <- p1
  
  ## gen2 sampling
  avail_genotypes <- sort(unique(index_1))
  if(length(avail_genotypes) == 1){
    index_2 <- avail_genotypes
    gt_ <- gt[,index_2]
    p2 <- 1 - gt_/2
  } else{
    genotype_frequency <- table(index_1) / N1[i]
    index_2 <- sample(avail_genotypes,N2[i],replace = T,prob = genotype_frequency)
    gt_ <- gt[,index_2]
    if(N2[i]==1){
      p2 <- 1 - gt_/2
    }else{
      p2 <- apply(gt_,1,function(x) 1 - sum(x)/N2[i]/2)
    }
  }
  p_drift_2[,i] <- p2
  
  ## gen3 sampling
  avail_genotypes <- sort(unique(index_2))
  if(length(avail_genotypes) == 1){
    index_3 <- avail_genotypes
    gt_ <- gt[,index_3]
    p3 <- 1 - gt_/2
  } else{
    genotype_frequency <- table(index_2) / N2[i]
    index_3 <- sample(avail_genotypes,N3[i],replace = T,prob = genotype_frequency)
    gt_ <- gt[,index_3]
    if(N3[i]==1){
      p3 <- 1 - gt_/2
    }else{
      p3 <- apply(gt_,1,function(x) 1 - sum(x)/N3[i]/2)
    }
  }
  
  p_drift_3[,i] <- p3 
  
  print(i)
}


## generation 1 frequency change
sim_gen_1_snp_deltap_avg <- c()
sim_gen_1_snp_deltap_var <- c()
for(i in 1:142){
  snp_p0 <- ld_p0$V3
  snp_p1 <- p_drift_1[,i]
  sim_gen_1_snp_deltap_avg <- c(sim_gen_1_snp_deltap_avg,mean(snp_p1-snp_p0))
  sim_gen_1_snp_deltap_var <- c(sim_gen_1_snp_deltap_var,var(snp_p1-snp_p0))
}
mean(sim_gen_1_snp_deltap_avg)
sd(sim_gen_1_snp_deltap_avg)
mean(sim_gen_1_snp_deltap_var)
sd(sim_gen_1_snp_deltap_var)

## generation 2 frequency change
sim_gen_2_snp_deltap_avg <- c()
sim_gen_2_snp_deltap_var <- c()

for(i in 1:142){
  snp_p0 <- ld_p0$V3
  snp_p1 <- p_drift_2[,i]
  sim_gen_2_snp_deltap_avg <- c(sim_gen_2_snp_deltap_avg,mean(snp_p1-snp_p0))
  sim_gen_2_snp_deltap_var <- c(sim_gen_2_snp_deltap_var,var(snp_p1-snp_p0))
}
mean(sim_gen_2_snp_deltap_avg)
mean(sim_gen_2_snp_deltap_var)
sd(sim_gen_2_snp_deltap_var)


## generation 3 frequency change
sim_gen_3_snp_deltap_avg <- c()
sim_gen_3_snp_deltap_var <- c()
for(i in 1:142){
  snp_p0 <- ld_p0$V3
  snp_p1 <-p_drift_3[,i]
  sim_gen_3_snp_deltap_avg <- c(sim_gen_3_snp_deltap_avg,mean(snp_p1-snp_p0))
  sim_gen_3_snp_deltap_var <- c(sim_gen_3_snp_deltap_var,var(snp_p1-snp_p0))
}
mean(sim_gen_3_snp_deltap_avg)
mean(sim_gen_3_snp_deltap_var)
sd(sim_gen_3_snp_deltap_var)


t.test(sim_gen_1_snp_deltap_avg,gen_1_snp_deltap_avg,alternative = "less")
t.test(sim_gen_2_snp_deltap_avg,gen_2_snp_deltap_avg,alternative = "less")
t.test(sim_gen_3_snp_deltap_avg,gen_3_snp_deltap_avg,alternative = "less")

foldchange_1 <- gen_1_snp_deltap_var/sim_gen_1_snp_deltap_var
t.test(foldchange_1)
foldchange_2 <- gen_2_snp_deltap_var/sim_gen_2_snp_deltap_var
t.test(foldchange_2)
foldchange_3 <- gen_3_snp_deltap_var/sim_gen_3_snp_deltap_var
t.test(foldchange_3)





