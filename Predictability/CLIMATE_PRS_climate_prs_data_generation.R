rm(list=ls())

library(dplyr)

## local adaptation

## This script generates this script generates the climate data for climate prs using 1001g genome
## It will generate 2 verstions: era5 and worldclim

prefix <- "LOCAL_ADAPTATION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

vcf_id <- read.delim("data-external/1001g_grenet_accessions_climate.txt",header=F)$V1
info <- read.delim("data/1001g_regmap_grenet_ecotype_info_corrected_2024May16.csv",sep=",")

## first check the ids are matched
table(info$ecotype_id[info$GrENE.net==T | info$source=="1TG"] %in% vcf_id)

worldclim <- read.delim("data-external/1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv",sep=",")
worldclim <- worldclim[,2:21] %>%
  filter(ecotypeid %in% vcf_id ) %>%
  arrange(match(ecotypeid,vcf_id))

era5 <- read.delim("data-external/bioclimvars_ecotypes_era5.csv",sep=",") %>%
  arrange(match(ecotype,vcf_id))

table(era5$ecotype == vcf_id)

NA_accessions <- info$ecotype_id[info$country=="USA" | info$GrENE.net==T]

era5[era5$ecotype %in% NA_accessions,2:20] <- "NA"
worldclim[worldclim$ecotypeid %in% NA_accessions,2:20 ] <- "NA"


training_ids <- era5$ecotype[era5$bio1_new != "NA"]
testing_ids <- era5$ecotype[era5$bio1_new == "NA"]
write.table(training_ids,"data-intermediate/climate_prs/training_ecotype_ids.txt",append = F,quote = F,col.names = F,row.names = F)
write.table(testing_ids,"data-intermediate/climate_prs/testing_ecotype_ids.txt",append = F,quote = F,col.names = F,row.names = F)


## 1001g prs data for bslmm
for(i in 1:19){
  name_era5 = paste0("data-intermediate/climate_prs/era5_bslmm_bio",i,".txt")
  write.table(era5[,i+1],file = name_era5,append = F,quote = F,row.names = F,col.names = F)
  name_worldclim = paste0("data-intermediate/climate_prs/worldclim_bslmm_bio",i,".txt")
  write.table(worldclim[,i+1],file = name_worldclim,append = F,quote = F,row.names = F,col.names = F)
}


## 1001g prs data for hapfm
for(i in 1:19){
  name_era5 = paste0("data-intermediate/climate_prs/era5_hapfm_bio",i,".txt")
  write.table(era5[!era5$ecotype %in% NA_accessions,i+1],file = name_era5,append = F,quote = F,row.names = F,col.names = F)
  name_worldclim = paste0("data-intermediate/climate_prs/worldclim_hapfm_bio",i,".txt")
  write.table(worldclim[!worldclim$ecotypeid %in% NA_accessions,i+1],file = name_worldclim,append = F,quote = F,row.names = F,col.names = F)
}


## The covariates data for hapfm
covariates <- read.delim("data-intermediate/climate_prs/1001g_grenet_climate_10PC.txt",header=F)
table(covariates$V2 == vcf_id)
training_covariates <- covariates[covariates$V2 %in% training_ids,c(1,3:12)]
write.table(training_covariates,"data-intermediate/climate_prs/TRAINING_10PC.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = F)

grenet_id <- info$ecotype_id[info$GrENE.net==T]
testing_grenet_covariates <- covariates[covariates$V2 %in% grenet_id,c(1,3:12)]
write.table(testing_grenet_covariates,"data-intermediate/climate_prs/TESTING_grenet_10PC.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = F)

