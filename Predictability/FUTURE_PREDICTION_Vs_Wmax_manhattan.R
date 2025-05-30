rm(list=ls())

library(dplyr)


## future prediction

## This script generates the manhattan plots for Vs and Wmax

prefix <- "FUTURE_PREDICTION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

source("scripts/manhattan_plot_Xing.R")

Wmax_co <- read.delim("data-intermediate/Wmax_covariates_formated.txt",header=F)
colnames(Wmax_co) <- c("SNP","CHR","BP","P")
Wmax_co$bon <- p.adjust(Wmax_co$P)
sig_snps <- Wmax_co$SNP[Wmax_co$bon<0.05]
png(filename ="figs/FUTURE_PREDICTION_Wmax_GWAS_manhattan.png",width = 8,height = 3,units = "in",res = 1000)
manhattan_new(Wmax_co,chr = "CHR",bp = "BP",p = "P",snp = "SNP",suggestiveline = F,genomewideline = F,highlight = sig_snps,col = c("grey40","grey80"),cex=0.6,ylim=c(0,max(-log10(Wmax_co$P))+0.5),
              xlab = "Wmax GWAS",xaxt = "n",ylab="-log10(p-value)")
dev.off()


Vs_co <- read.delim("data-intermediate/Vs_covariates_formated.txt",header=F)
colnames(Vs_co) <- c("SNP","CHR","BP","P")
Vs_co$bon <- p.adjust(Vs_co$P)
sig_snps <- Vs_co$SNP[Vs_co$bon<0.05]
png(filename ="figs/FUTURE_PREDICTION_Vs_GWAS_manhattan.png",width = 8,height = 3,units = "in",res = 1000)
manhattan_new(Vs_co,chr = "CHR",bp = "BP",p = "P",snp = "SNP",suggestiveline = F,genomewideline = F,highlight = sig_snps,col = c("grey40","grey80"),cex=0.6,ylim=c(0,max(-log10(Vs_co$P))+0.5),
              xlab = "Vs GWAS",xaxt = "n",ylab="-log10(p-value)")
dev.off()

