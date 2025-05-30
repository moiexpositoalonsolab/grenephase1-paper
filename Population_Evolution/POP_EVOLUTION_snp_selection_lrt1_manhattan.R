rm(list=ls())

## POPULATION EVOLUTION

## This script will use the the manhattan_plot_Xing.R to plot manhattan plots for snp selection signal 
## results (LRT1 test) in  generation1 across all sites 

prefix <- "POP_EVOLUTION_"
PATH <-  "~/xingwu@berkeley.edu - Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)
source("scripts/manhattan_plot_Xing.R")

af <- read.delim("data-intermediate/snp_allele_frequency.txt",header=F)
locations <- read.delim("data-intermediate/chromosome_locations.txt",header=F)
locations <- locations[af$V1>0.05,]
site_list <- read.delim("data-intermediate/list.txt",header=F)

for(i in 1:length(site_list$V1)){
  name <- site_list[i,]
  pvalue <- read.delim(paste("data-intermediate/snp_lrt1/gen1_",name,"_pvalue.txt",sep = ""),header=F)
  pvalue$V1[pvalue$V1==0] <- runif(sum(pvalue$V1==0),min = 10^-20,max=10^-16)
  data <- as.data.frame(matrix(nrow=nrow(pvalue),ncol=4))
  colnames(data) <- c("CHR","BP","SNP","P")
  data$CHR <- locations[,1]
  data$BP <- locations[,2]
  data$SNP <- locations[,3]
  data$P <- pvalue$V1
  sig_snps <- data$SNP[data$P<0.05/1174429/30]
  output_name <- paste("figs/snp_lrt1/","snp_lrt1_generation1_",name,".png",sep="")
  png(filename =output_name,width = 8,height = 0.7,units = "in",res = 1000)
  par(mar = c(1, 1, 1, 1))
  manhattan_new(data,suggestiveline = F,genomewideline = F,highlight = sig_snps,col = c("grey40","grey80"),cex=0.6,ylim=c(0,20+0.5),
                xlab = "",xaxt = "n",ylab=name,yaxt="n")
  dev.off()
  print(i)
}



