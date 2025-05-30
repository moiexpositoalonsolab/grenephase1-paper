rm(list=ls())

library(dplyr)

## climate prs

## This script will make the ecotype optimal climate dataframe using their PRS scores

prefix <- "CLIMATE_GWAS_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)


vcf_id <- read.delim("data-external/1001g_grenet_accessions_climate.txt",header=F)$V1
info <- read.delim("data/1001g_regmap_grenet_ecotype_info_corrected_2024May16.csv",sep=",")
grenet_index <- which(vcf_id %in%info$ecotype_id[info$GrENE.net==T])


grenet_ecotypes <- read.delim("data/founder_ecotype_frequency.txt",header=F)

## First check the order of grenet ecotype id does match
table(vcf_id[grenet_index] == grenet_ecotypes$V1)

## load the era5 data
era5 <- read.delim("data-external/bioclimvars_ecotypes_era5.csv",sep=",") %>%
  arrange(match(ecotype,vcf_id))
era5_grenet <- era5[grenet_index,]

## build the prs table
prs_bio <- as.data.frame(matrix(nrow=231,ncol = 20))
colnames(prs_bio) <- c("ecotype",paste0("PRS_bio",1:19))
prs_bio$ecotype <- vcf_id[grenet_index]


for(i in 1:19){
  prd <- read.delim(paste0("data-intermediate/climate_prs/era5_bslmm_bio",i,".prdt.txt"),header=F)
  m <- lm(era5_grenet[,i+1]~prd$V1[grenet_index])
  prs_bio[,i+1] <- m$coefficients[1]+prd$V1[grenet_index]
}

## check with the real era5 data
plot(prs_bio$PRS_bio4,era5_grenet$bio4_new)
abline(a = 0,b=1)
summary(lm(era5_grenet$bio4_new~prs_bio$PRS_bio4))

write.table(prs_bio,paste0("data-intermediate/",prefix,"bslmm_prs_ecotype_climate.txt"),append = F,quote = F,sep = "\t",row.names = F,col.names = T)




