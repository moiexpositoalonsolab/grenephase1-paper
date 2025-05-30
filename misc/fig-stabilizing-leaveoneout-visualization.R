################################################################################
### Goal
### Estimate some basic metrics of ecotypes increasing in frequency across
### sites, means, variances, optimal places

################################################################################
library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
theme_set(theme_minimal())

################################################################################
#### Location
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
# Leave one out 
dat<-read.table("data-intermediate/leave_one_out_comparison_all.txt", header=T)
head(dat)
# dat<-dat %>% dplyr::filter(source=="climate_distance")
# dat<-dat %>% dplyr::filter(source=="stabilizing_climate")
dat<-dat %>% dplyr::filter(source=="climate_PRS")
dat<-
  dat %>% 
  split
  

# Merge with climate information
load("grene/data/worldclim_sitesdata.rda")
load("grene/data/locations_data.rda")
load("data/flowers_survival_long.rda")
# Parse flowers to be able to merge
flowers_survival_long %>% 
  goup

dat<-
dat %>% 
  merge(.,worldclim_sitesdata,by="site") %>% 
  merge(., locations_data,by="site")

# Do visualizations
ggplot(dat)+
  stat_summary(aes(y=sp_r, x= bio1, color=longitude> -15))
ggplot(dat)+
  stat_summary(aes(y=pearson_r, x= bio1, color=longitude> -15))
ggplot(dat)+
  stat_summary(aes(y=r2, x= bio1, color=longitude> -15))
ggplot(dat)+
  stat_summary(aes(y=pearson_r, x= bio12, color=longitude> -15))
ggplot(dat)+
  stat_summary(aes(y=pearson_r, x= bio18, color=longitude> -15))


cor.test(dat$pearson_r, dat$bio1)
cor.test(dat$pearson_r, dat$bio18)


