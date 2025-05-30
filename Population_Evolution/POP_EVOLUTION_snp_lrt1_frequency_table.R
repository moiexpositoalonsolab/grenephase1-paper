rm(list=ls())

library(ggplot2)
library(data.table)


## POPULATION EVOLUTION

## This script will use the delta snp frequency to generate the diffusion map and PCA plots to
## describe the population change dynamics

prefix <- "POP_EVOLUTION_"
PATH <-  "~/xingwu@berkeley.edu - Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

af <- read.delim("data-intermediate/snp_allele_frequency.txt",header=F)
maf0.05_index <- af$V1 > 0.05

p0 <- fread("data/average_seedmix_p0.txt",header=F)
p0 <- as.numeric(p0$V3)
p0 <- p0[maf0.05_index]

freq <- fread("data/merged_hapFIRE_allele_frequency.csv",header=T,sep=",",stringsAsFactors = F,nThread = 10)
freq_matrix <- as.matrix(freq)
freq_matrix <- freq_matrix[maf0.05_index,]

dim(freq_matrix)

terminal_index <- read.delim("data/site_plot_terminal_index.txt")
## generation 1
index_gen_1 <- terminal_index$index1[terminal_index$generation1==1]
meta <- read.delim("data/merged_sample_table.csv",sep=",")
meta_gen1 <- meta[index_gen_1,]


freq_matrix_generation1 <- freq_matrix[,index_gen_1]
colnames(freq_matrix_generation1)

for( i in 1:ncol(freq_matrix_generation1)){
  index <- which(freq_matrix_generation1[,i]>=1)
  if(length(index > 0)){
    print(i)
    freq_matrix_generation1[index,i] <- 0.9999999999
  }
}

site <- unique(meta$site[index_gen_1])

for(i in site){
  index <-  which(meta_gen1$site==i)
  if (length(index) >= 2){
    name <- paste("data-intermediate/snp_lrt1/gen1_site_",i,"_freq.csv",sep="")
    fwrite(freq_matrix_generation1[,index],name,append =F,row.names = F,col.names = F, quote = F,sep=",",nThread = 7)
    print(i)
  }
}
}
