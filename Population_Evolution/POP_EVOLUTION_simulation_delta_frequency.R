rm(list=ls())

library(ggplot2)
library(data.table)
library(SNPRelate)
library(SeqArray)
library(vcfR)
library(dplyr)
library(patchwork)

## POPULATION EVOLUTION

## This script will simulate GrENE-net under neutral / selection situations and plot the delta frequency in two ways: 
## 1) all compare to p0 p1-p0, p2-p0 2) per generation change p1-p0, p2-p1


prefix <- "POP_EVOLUTION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

## The real greneNET data to use random sampling to estimate the null distribution

meta <- read.delim("data/merged_sample_table.csv",sep=",",stringsAsFactors = F)
terminal_index <- read.delim("data/site_plot_terminal_index.txt")
index_gen_3 <- which(terminal_index$generation1==1 &terminal_index$generation2==2 & terminal_index$generation3==3)

ld_p0 <- read.delim("data/average_seedmix_p0_LDpruned.txt",header=F)
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



snp_delta <- as.data.frame(matrix(nrow=142*3,ncol=2))
colnames(snp_delta) <- c("delta_frequency","generation")

snp_delta$delta_frequency <- c(sim_gen_1_snp_deltap_avg,sim_gen_2_snp_deltap_avg,sim_gen_3_snp_deltap_avg)
snp_delta$generation <- c(rep(1,142),rep(2,142),rep(3,142))


means <- snp_delta %>%
  group_by(generation) %>%
  dplyr::summarize(mean_value = mean(delta_frequency))

p1 <- ggplot(snp_delta,aes(x=delta_frequency))+
  geom_histogram(alpha = 0.9,position = 'identity',bins = 50,fill="grey40")+
  geom_density(alpha=0.8,fill="grey80")+
  geom_vline(data = means,aes(xintercept = mean_value),linetype="dashed", color = "black")+
  geom_text(data = means, aes(x = -0.05, y = 200, label = paste("Mean =", round(mean_value, 5))), vjust = -1,cex=5) +
  facet_wrap(~ generation) +
  ylab("Frequency")+
  xlim(c(-0.1,0.1))+
  ylim(c(0,250))+
  theme_minimal(base_size = 18)

p1
ggsave(plot = p1,filename =  paste0("figs/",prefix,"simulated_grene_net_snp_delta_frequency_ALLTOFOUNDER.pdf"),device = "pdf",width = 10,height = 4)


##### calculate the per generation delta frequency: per generation change p1-p0, p2-p1 .. 

## generation 1 frequency change
gen_1_snp_deltap_avg <- c()
for(i in 1:500){
  snp_p0 <- ld_p0$V3
  snp_p1 <- p_drift_1[,i]
  gen_1_snp_deltap_avg <- c(gen_1_snp_deltap_avg,mean(snp_p1-snp_p0))
}
mean(gen_1_snp_deltap_avg)


## generation 2 frequency change
gen_2_snp_deltap_avg <- c()
for(i in 1:500){
  snp_p0 <- p_drift_1[,i]
  snp_p1 <- p_drift_2[,i]
  gen_2_snp_deltap_avg <- c(gen_2_snp_deltap_avg,mean(snp_p1-snp_p0))
}
mean(gen_2_snp_deltap_avg)
var(gen_2_snp_deltap_avg)



## generation 3 frequency change
gen_3_snp_deltap_avg <- c()
for(i in 1:500){
  snp_p0 <- p_drift_2[,i]
  snp_p1 <-p_drift_3[,i]
  gen_3_snp_deltap_avg <- c(gen_3_snp_deltap_avg,mean(snp_p1-snp_p0))
}
mean(gen_3_snp_deltap_avg)

snp_delta_2 <- as.data.frame(matrix(nrow=500*3,ncol=2))
colnames(snp_delta_2) <- c("delta_frequency","generation")

snp_delta_2$delta_frequency <- c(gen_1_snp_deltap_avg,gen_2_snp_deltap_avg,gen_3_snp_deltap_avg)
snp_delta_2$generation <- c(rep(1,500),rep(2,500),rep(3,500))

means_2 <- snp_delta_2 %>%
  group_by(generation) %>%
  dplyr::summarize(mean_value = mean(delta_frequency))

p2 <- ggplot(snp_delta_2,aes(x=delta_frequency))+
  geom_histogram(alpha = 0.9,position = 'identity',bins = 50,fill="grey40")+
  geom_density(alpha=0.8,fill="grey80")+
  geom_vline(data = means_2,aes(xintercept = mean_value),linetype="dashed", color = "black")+
  geom_text(data = means_2, aes(x = -0.05, y = 200, label = paste("Mean =", round(mean_value, 5))), vjust = -1) +
  facet_wrap(~ generation) +
  ylab("Frequency")+
  xlim(c(-0.1,0.1))+
  ylim(c(0,250))+
  theme_minimal()

combined <- p1 / p2
ggsave(plot = combined,filename =  paste0("figs/",prefix,"simulated_grene_net_snp_delta_frequency.pdf"),device = "pdf",width = 10,height = 8)
