################################################################################
### Goal
### Visualize some well known gnes

################################################################################
#### If transformed to collab can use the stuff below
# from google.colab import drive
# drive.mount('/content/drive')
# %load_ext rpy2.ipython
# %%R

################################################################################
#### Libraries
library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(patchwork)

################################################################################
#### Location if run locally
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
# myfolder<-"~/grenephase1-analyses"

myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)



load(file = "data-intermediate/alleles-LDpruned-TAIR10_parsedgenes.rda")
load(file = "data-intermediate/allele_allyears_frequencies_long_raw_climate_flowers.rda")

annotated_snps %>% head
a<-allele_allyears_frequencies_long_raw_climate_flowers
  a %>% head

################################################################################
### FLC  AT5G10140
### GIGANTEA  AT1G22770
### Gigantia  AT1G22770
mygenes<-read.csv("data-external/common-genes-and-identifiers.csv")


annotated_snps %>% 
  dplyr::filter(grepl( "AT5G10140",closest_gene)) #flc

annotated_snps %>% 
  dplyr::filter(grepl( "AT1G22770",closest_gene)) # gigantea

annotated_snps %>% 
  dplyr::filter(grepl( "AT1G65480",closest_gene)) # ft

annotated_snps %>% 
  dplyr::filter(grepl( "AT3G58670",closest_gene)) # pco5

annotated_snps %>% 
  dplyr::filter(grepl( "AT5G22090",closest_gene)) # ear1

annotated_snps %>% 
  dplyr::filter(grepl( "AT5G45830",closest_gene)) # dog1


################################################################################
### Have a look at a snp close to FT1

asub<-
  a %>% 
  # dplyr::filter(snp=="1_24334966") # ft
  # dplyr::filter(snp=="1_8062999") # gigante
  dplyr::filter(snp=="5_7319807") # ear1
  

ggplot(asub)+
  geom_point(aes(y=maxfreq, x=latitude))+
  stat_summary(aes(y=maxfreq, x=latitude, group=site))+
  stat_smooth(aes(y=maxfreq, x=latitude),method = "glm")

ggplot(asub)+
  geom_point(aes(y=maxfreq, x=julian_day))+
  stat_summary(aes(y=maxfreq, x=julian_day, group=site))+
  stat_smooth(aes(y=maxfreq, x=julian_day),method = "glm")
