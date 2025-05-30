################################################################################
### Goal
### Estimate selection coefficient of top SNPs
################################################################################

library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives/MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses" # GOOGLE DRVIE
setwd(myfolder)

################################################################################

# test Ne estimation
sub<-read_merged_allele_subset(
  subsetfile = "data-intermediate/merged_hapFIRE_allele_frequency-AT2G22540-SVP.csv",
  headfile="data-intermediate/merged_hapFIRE_allele_frequency-header.csv",
  snpfile="data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-AT2G22540-SVP.csv",
  addclimate = T, addflowers = T
)
head(sub)
sub$flower %>% unique
sub$year %>% unique
# View(head(sub))

subsub<-sub %>% 
  dplyr::filter(site==4, year %in% c(1,3)) %>% 
  dplyr::filter(rep==1)
