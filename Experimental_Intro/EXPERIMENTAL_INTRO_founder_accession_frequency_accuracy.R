rm(list=ls())

library(ggplot2)
library(patchwork)
library(dplyr)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)

## EXPERIMENTAL_INTRO

## This script estimate the allele and accession frequencies of 231 founders and plot the correlation between eight replicates

prefix <- "EXPERIMENTAL_INTRO_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

df <- as.data.frame(matrix(nrow=231,ncol=8))

for(i in 1:8){
  tmp <- read.delim(paste0("data/seed_mix/s",i,"_ecotype_frequency.txt"),header=F)
  df[,i] <- tmp$V2
}

founder_freq_avg <- apply(df,1,mean)

mat <- cor(df,method = "pearson")
mean(mat[upper.tri(mat)])
quantile(mat[upper.tri(mat)],0.25)
quantile(mat[upper.tri(mat)],0.75)
sqrt(mean((founder_freq_avg-1/231)^2))
mean(abs(founder_freq_avg-1/231))


df1 <- tibble(accession = 1:231, freq = founder_freq_avg)
p1 <- ggplot(df1,aes(x=accession,y=freq))+
  geom_hline(yintercept = 1/231,colour = "red",linewidth = 1,linetype = 1)+
  geom_point(pch=16,cex=2)+
  ylim(c(0,0.1))+
  ylab("accession frequency")+
  xlab("GrENE-net founder")+
  theme_minimal(base_size = 18)
p1

ggsave(plot = p1,filename = "figs/EXPERIMENTAL_INTRO_founder_accession_frequency_observed_expexted.pdf",device = "pdf",width = 8,height = 8)
ggsave(plot = p1,filename = "figs/EXPERIMENTAL_INTRO_founder_accession_frequency_observed_expexted.png",device = "png",width = 8,height = 8,dpi=600,units="in",bg = "white")


## read the founder allele frequencies

allele <- as.data.frame(matrix(nrow=3235480,ncol=8))

for(i in 1:8){
  tmp <- read.delim(paste0("data/seed_mix/s",i,"_snp_frequency.txt"),header=F)
  allele[,i] <- tmp$V3
}

colnames(allele) <- paste0("seed_replicate_",1:8)

mat <- cor(allele,method = "pearson")


pdf("figs/EXPERIMENTAL_INTRO_seed_mix_allele_correlation.pdf",width = 8,height = 8,useDingbats = F)
corrplot(mat,diag = T,method = "ellipse",type="upper")
dev.off()

allele_af_seed <- apply(allele,1,mean)

## read the vcf file

## load the VCF file and convert it to snp matrix
seqVCF2GDS("data/greneNet_final_v1.1.recode.vcf", "data-intermediate/greneNet_final_v1.1.gds")

genofile <- seqOpen("data-intermediate/greneNet_final_v1.1.gds")
snp_matrix <- as.matrix(seqGetData(genofile, "$dosage"))
variant.id <- paste0(1,"_",seqGetData(genofile, "position"))
seqClose(genofile)

allele_af_vcf <- 1- t(snp_matrix) %*% founder_freq_avg /2

df_2 <- tibble(observed =allele_af_seed,
               expected =allele_af_vcf )

p2 <- ggplot(df_2,aes(x=observed,y=expected))+
  geom_point(pch=16,color="grey50",cex=1.5)+
  geom_abline(slope = 1,intercept = 0,color="red3",lty=2)+
  xlab("Estimated allele frequency from seed mixture")+
  ylab("Expected allele frequency from the founder VCF")+
  annotate(geom="text",label= "Pearson correlation r: 0.9998",x = 0.8,y=0.2,cex=6)+
  theme_minimal(base_size = 18)
p2
ggsave(plot = p2,filename = "figs/EXPERIMENTAL_INTRO_founder_allele_frequency_observed_expexted.pdf",device = "pdf",width = 8,height = 8)
ggsave(plot = p2,filename = "figs/EXPERIMENTAL_INTRO_founder_allele_frequency_observed_expexted.png",device = "png",width = 8,height = 8,dpi=600,units="in",bg = "white")
