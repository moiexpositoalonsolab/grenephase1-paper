rm(list=ls())

library(dplyr)
library(RColorBrewer)
library(ggplot2)

## EXPERIMENTAL_INTRO

## This script plot the climate PCA of 1TG ecotypes, 231 founders and 43 sites

prefix <- "EXPERIMENTAL_INTRO_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)


ecotype_climate <- read.delim("data-external/bioclimvars_ecotypes_era5.csv",sep=",")
site_climate <- read.delim("data-external/bioclimvars_experimental_sites_era5.csv",sep=",")

ecotype_info <- read.delim("data/1001g_regmap_grenet_ecotype_info_corrected_2024May16.csv",sep=",")
site_info <- read.delim("data/merged_sample_table.csv",sep=",")


grenet_climate <- ecotype_climate %>%
  filter(ecotype %in% ecotype_info$ecotype_id[ecotype_info$GrENE.net==T])


## count the number of years each site has survived
survival <- c()
for(i in 1:nrow(site_climate)){
  site_ <- site_climate$site[i]
  survival <- c(survival,max(site_info$generation[site_info$site==site_]))
}

survival[survival==-Inf] <- 0
site_climate$survival <- as.factor(survival)

p <- ggplot(grenet_climate,aes(x = bio1_new,y=bio12_new))+
  geom_point(pch=3,cex=3,color="grey50")+
  xlab("Annual Mean Temperature (Celsius)") +
  ylab("Annual Mean Percipitation (mm)")+
  geom_point(data = site_climate,mapping = aes(x = bio1,y = bio12,color=survival),pch=16,alpha=0.7,cex=5)+
  scale_color_manual(values=c("salmon","#C7E9C0","#74C476","#006D2C")) + 
  theme_classic()
p
ggsave(p,file=paste0("figs/",prefix,"site_ecotype_climate_overview.pdf"),width = 5,height = 4)
