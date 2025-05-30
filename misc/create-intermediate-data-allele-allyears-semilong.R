
################################################################################
### Goal
### Create intermediate data for the LD prunned alleles long and wid
### with climate and with population size

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

################################################################################
# Specify the file path
file_edelta = '/data/delta_ecotype_freq.txt'
file_e0 = '/data/merged_ecotype_frequency.txt'
file_eid = '/data/greneNet_final_v1.0.fam'
file_climate_sites="/grene/data/worldclim_sitesdata.rda"
file_sites_info="/grene/data/sites.clim.rda"

file_pdelta = '/data/merged_hapFIRE_allele_frequency_LDpruned.txt'
file_p0 = '/data/average_seedmix_p0_LDpruned.txt'

################################################################################
# Load frequencies 
# pdelta<-read.table("data/merged_hapFIRE_delta_p_LDpruned.txt", header=T)
p<-read.table("data/merged_hapFIRE_allele_frequency_LDpruned.txt", header=T)
p0<-read.table("data/average_seedmix_p0_LDpruned.txt", header=F)
p0$snp<-paste0(p0$V1,"_",p0$V2)
p$snp<-p0$snp
p0<-p0 %>% dplyr::rename(chr=V1,pos=V2,startfreq=V3)

################################################################################
# pivot partially 
pmidlong<-
  p %>%
  # head %>%
  pivot_longer(
    cols = -snp,
    names_to = c("site", "year", "rep"),
    names_pattern = "X(\\d+)_(\\d+)_(\\d+)",
    values_to = "freq"
  ) %>% 
  dplyr::select(site, year, rep, freq, snp) %>% 
  pivot_wider(
    names_from = year,
    values_from = freq,
    names_prefix = "y"
  ) %>% 
  merge(., p0,by='snp') %>% 
  dplyr::rename(y0=startfreq) %>% 
  dplyr::mutate(y3=y3-y2) %>%  # needs to be reversed so that the calculation is correct
  dplyr::mutate(y2=y2-y1) %>% 
  dplyr::mutate(y1=y1-y0)

pmidlong<-
  pmidlong %>% 
  dplyr::select(site,rep,snp,chr,pos,y0,y1,y2,y3)
  

allele_allyears_frequencies_semilong_raw = a = pmidlong
head(a)

message("saving alleles all years frequency, long format, raw")
write.csv( file= paste0(myfolder,"/data-intermediate/allele_allyears_frequencies_semilong_raw.csv"),
           allele_allyears_frequencies_semilong_raw)
save(file=paste0(myfolder,"/data-intermediate/allele_allyears_frequencies_semilong_raw.rda"),
     allele_allyears_frequencies_semilong_raw)


################################################################################

