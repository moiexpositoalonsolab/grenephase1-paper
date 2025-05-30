rm(list=ls())

library(dplyr)
library(ggplot2)
library(maps)

## climate prs

## This script will compare worldclim and era5 two climate data on climate prs
## The PRS methods are bslmm and hapfm

prefix <- "CLIMATE_GWAS_"
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



world_map <- map_data("world")
TG <- info[info$source == "1TG",]
lon_range <- range(TG$Longitude_corrected)
lat_range <- range(TG$Latitude_corrected) + c(-10, 10)

ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),fill = "grey90", color = "white") +
  #geom_point(data = TG[TG$GrENE.net==T,], aes(x = Longitude_corrected, y = Latitude_corrected), color = "red", size = 1) +
  geom_point(data = TG[TG$GrENE.net==F,], aes(x = Longitude_corrected, y = Latitude_corrected), color = "blue", size = 2) +
  theme_minimal() +
  coord_fixed(xlim = lon_range, ylim = lat_range) +
  labs(x = "Longitude", y = "Latitude")+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20))
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),fill = "grey90", color = "white") +
  geom_point(data = TG[TG$GrENE.net==T,], aes(x = Longitude_corrected, y = Latitude_corrected), color = "red", size = 2) +
  #geom_point(data = TG[TG$GrENE.net==F,], aes(x = Longitude_corrected, y = Latitude_corrected), color = "blue", size = 1) +
  theme_minimal() +
  coord_fixed(xlim = lon_range, ylim = lat_range) +
  labs(x = "Longitude", y = "Latitude")+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20))



table(era5$ecotype == vcf_id)
table(era5$ecotype == worldclim$ecotypeid)

grenet_index <- which(vcf_id %in%info$ecotype_id[info$GrENE.net==T])

bioclim_summary <- as.data.frame(matrix(nrow=38,ncol=6))
colnames(bioclim_summary) <- c("bioclim","source","pve","pge","n_gamma","r2")
bioclim_summary$bioclim <- rep(paste0("bio",1:19),2)
bioclim_summary$source <- c(rep("worldclim",19),rep("era5",19))

## lets first look at the estimated heritability and polygenicity from different climate data
for(i in 1:19){
  hyp <- read.delim(paste0("data-intermediate/climate_prs/worldclim_bslmm_bio",i,".hyp.txt"),sep="\t",skip = 1,header=F)[,1:6]
  colnames(hyp) <- c("h","pve","rho","pge","pi","n_gamma")
  bioclim_summary[i,3] <- mean(tail(hyp$pve,n=10000))
  bioclim_summary[i,4] <- mean(tail(hyp$pge,n=10000))
  bioclim_summary[i,5] <- round(mean(tail(hyp$n_gamma,n=10000)))
  
  hyp <- read.delim(paste0("data-intermediate/climate_prs/era5_bslmm_bio",i,".hyp.txt"),sep="\t",skip = 1,header=F)[,1:6]
  colnames(hyp) <- c("h","pve","rho","pge","pi","n_gamma")
  bioclim_summary[i+19,3] <- mean(tail(hyp$pve,n=10000))
  bioclim_summary[i+19,4] <- mean(tail(hyp$pge,n=10000))
  bioclim_summary[i+19,5] <- round(mean(tail(hyp$n_gamma,n=10000)))
}


## Now lets look the BSLMM prediction on different datasets
era5_true <- era5[grenet_index,2:20]
worldclim_true <- worldclim[grenet_index,2:20]

for(i in 1:19){
  prd <- read.delim(paste0("data-intermediate/climate_prs/worldclim_bslmm_bio",i,".prdt.txt"),header=F)
  bioclim_summary[i,6] <- summary(lm(worldclim_true[,i]~prd$V1[grenet_index]))$adj.r.squared
  
  prd <- read.delim(paste0("data-intermediate/climate_prs/era5_bslmm_bio",i,".prdt.txt"),header=F)
  bioclim_summary[i+19,6] <- summary(lm(era5_true[,i]~prd$V1[grenet_index]))$adj.r.squared
}

plot(prd$V1[grenet_index]+9.901,era5_true[,i],pch=16,cex=2,xlim =c (-5,20),ylim=c(-15,25))
abline(a=0,b=1,lty=2,col="red",lwd=2)

ggplot(bioclim_summary,aes(x=bioclim,y = r2,col=source))+
  geom_boxplot()


## Now lets look at hapfm's prediction on the era5 data

hapfm_summary <- as.data.frame(matrix(nrow = 57,ncol=3))
colnames(hapfm_summary) <- c("bioclim","r2","source")
hapfm_summary$bioclim <- c(rep(paste0("bio",1:19),2),paste0("bio",1:19))
hapfm_summary$source <- c(rep("KNN",19),rep("xmeans",19),rep("modularity",19))


## Load the covariates
C <- read.delim("data-intermediate/climate_prs/TESTING_grenet_10PC.txt",header = F)

## KNN first 

KNN_dm <- read.delim("data-intermediate/climate_prs/TESTING_1001g_grenet_climate_KNN_haplotypeDM.txt",check.names = F)
testing_ids <- read.delim("data-intermediate/climate_prs/testing_ecotype_ids.txt",header=F)

testing_grenet <-testing_ids$V1[testing_ids$V1 %in% info$ecotype_id[info$GrENE.net==T]]

grenet_dm <- KNN_dm[testing_ids$V1 %in% info$ecotype_id[info$GrENE.net==T],]

table(testing_grenet == vcf_id[grenet_index])

for(i in 1:19){
  alpha <- read.delim(paste0("data-intermediate/climate_prs/era5_bio",i,"_KNN_10PC_alpha.txt"),header=F)
  beta <- read.delim(paste0("data-intermediate/climate_prs/era5_bio",i,"_KNN_10PC_beta.txt"),header=F)
  prd <- as.matrix(C) %*% alpha$V1 + as.matrix(grenet_dm) %*% beta$V2
  hapfm_summary[i,2] <- summary(lm(era5_true[,i]~prd))$adj.r.squared
  
  alpha <- read.delim(paste0("data-intermediate/climate_prs/era5_bio",i,"_KNN_alpha.txt"),header=F)
  beta <- read.delim(paste0("data-intermediate/climate_prs/era5_bio",i,"_KNN_beta.txt"),header=F)
  prd <- alpha$V1 + as.matrix(grenet_dm) %*% beta$V2
  hapfm_summary[i,2] <- summary(lm(era5_true[,i]~prd))$adj.r.squared
}


## xmeans then

xmeans_dm <- read.delim("data-intermediate/climate_prs/TESTING_1001g_grenet_climate_xmeans_haplotypeDM.txt",check.names = F)
testing_ids <- read.delim("data-intermediate/climate_prs/testing_ecotype_ids.txt",header=F)

testing_grenet <-testing_ids$V1[testing_ids$V1 %in% info$ecotype_id[info$GrENE.net==T]]

grenet_xmeans_dm <- xmeans_dm[testing_ids$V1 %in% info$ecotype_id[info$GrENE.net==T],]

table(testing_grenet == vcf_id[grenet_index])

for(i in 1:19){
  alpha <- read.delim(paste0("data-intermediate/climate_prs/era5_bio",i,"_xmeans_alpha.txt"),header=F)
  beta <- read.delim(paste0("data-intermediate/climate_prs/era5_bio",i,"_xmeans_beta.txt"),header=F)
  prd <- alpha$V1 + as.matrix(grenet_xmeans_dm) %*% beta$V2
  hapfm_summary[i+19,2] <- summary(lm(era5_true[,i]~prd))$adj.r.squared
}


## modularity last

modularity_dm <- read.delim("data-intermediate/climate_prs/TESTING_1001g_grenet_climate_modularity_haplotypeDM.txt",check.names = F)
testing_ids <- read.delim("data-intermediate/climate_prs/testing_ecotype_ids.txt",header=F)

testing_grenet <-testing_ids$V1[testing_ids$V1 %in% info$ecotype_id[info$GrENE.net==T]]

grenet_modularity_dm <- modularity_dm[testing_ids$V1 %in% info$ecotype_id[info$GrENE.net==T],]

table(testing_grenet == vcf_id[grenet_index])

for(i in 1:19){
  alpha <- read.delim(paste0("data-intermediate/climate_prs/era5_bio",i,"_modularity_alpha.txt"),header=F)
  beta <- read.delim(paste0("data-intermediate/climate_prs/era5_bio",i,"_modularity_beta.txt"),header=F)
  prd <- alpha$V1 + as.matrix(grenet_modularity_dm) %*% beta$V2
  hapfm_summary[i+38,2] <- summary(lm(era5_true[,i]~prd))$adj.r.squared
}

ggplot(hapfm_summary,aes(x=bioclim,y = r2,col=source))+
  geom_boxplot()

combined_summary <- rbind(bioclim_summary[,c(1,6,2)],hapfm_summary)

combined_summary %>%
  filter(source %in% c("era5","KNN","modularity")) %>%
  ggplot(.,aes(x=bioclim,y = r2,col=source))+
  geom_boxplot()
p
