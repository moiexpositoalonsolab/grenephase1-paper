rm(list=ls())

library(ggplot2)
library(data.table)


## POPULATION EVOLUTION

## This script will use ecotype frequency data to generate the phenotype input for survival gwas

prefix <- "POP_EVOLUTION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

meta <- read.delim("data/merged_sample_table.csv",sep=",",stringsAsFactors = F)
terminal_index <- read.delim("data/site_plot_terminal_index.txt")

## ecotype frequencies
ecotype_p0 <- read.delim("data/founder_ecotype_frequency.txt",header=F)
ecotype_p <- read.delim("data/merged_ecotype_frequency.txt",header=T,check.names = F)

## site 4,5,32

sites <- c(4,5,32)
for(i in sites){
  site_ <- i
  ## generation 1 
  index <- which(meta$site == i & meta$generation==1)
  delta_p <- apply(ecotype_p[,index],2,function(x) (x - ecotype_p0$V2))
  mean_delta_p <- apply(delta_p,1,mean)
  name <- paste0("data-intermediate/survival_gwas/site_",i,"_generation_1_frequency_change.txt")
  write.table(mean_delta_p,name,append =F,row.names = F,col.names = F, quote = F)
  
  
  ## generation 3
  index <- which(meta$site == i & meta$generation==3)
  delta_p <- apply(ecotype_p[,index],2,function(x) (x - ecotype_p0$V2))
  mean_delta_p <- apply(delta_p,1,mean)
  name <- paste0("data-intermediate/survival_gwas/site_",i,"_generation_3_frequency_change.txt")
  write.table(mean_delta_p,name,append =F,row.names = F,col.names = F, quote = F)
  
}
