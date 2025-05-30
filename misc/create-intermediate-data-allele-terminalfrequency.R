
#### If in google drive
# from google.colab import drive
# drive.mount('/content/drive')
# %load_ext rpy2.ipython

################################################################################
library(sp)
library(raster)
library(rbioclim)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(rgbif)
library(dplyr)
library(moiR)
library(tidyverse)


################################################################################
# Set map limits
xlim = c(-200, +200)
ylim = c(0, +90)
xlim = c(-200, +200)
ylim = c(0, +90)
xlim=c(-10.5,+ 53)
ylim=c(32,65)
Range=extent(c(xlim,ylim))

# Download present data w2
now<-rbioclim::getData(name="worldclim2",var='bio', res=2.5, path = "~/")
now<- now %>% crop(.,Range)

# Download future data Max Planck model
fut<-rbioclim::getData(name="CMIP5",var='bio', res=2.5,model='MP', year=50, rcp=85,path = "~/")
fut<- fut %>% crop(.,Range)
names(fut) <- names(now) # Fix names in future dataset
for(i in 1:11) fut[[i]]<-fut[[i]]/10 # Fix units of temp in future dataset (wordlclim2 temp is not Cx10)


################################################################################
# %%R

# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
myfolder<-"~/grenephase1-analyses"

# Specify the file path

library(moiR)

file_edelta = '/data/delta_ecotype_freq.txt'
file_e0 = '/data/merged_ecotype_frequency.txt'
file_climate_sites="/grene/data/worldclim_sitesdata.rda"
file_sites_info="/grene/data/sites.clim.rda"

file_path_pdelta = paste0(myfolder, '/data/merged_hapFIRE_allele_frequency_LDpruned.txt' )
file_path_p0 = paste0(myfolder,'/data/average_seedmix_p0_LDpruned.txt')

file_ecotypes_wodclim="/grene/data/worldclim_ecotypesdata.rda"
file_ecotypes_data="/grene/data/ecotypes_data.rda"


myallele="1_17199262" # example from allele-climate correlation
################################################################################

# Load frequencies 
pdelta<-read.table(file_path_pdelta, header=T)
p0<-read.table(file_path_p0, header=F)
p0$snp<-paste0(p0$V1,"_",p0$V2)
pdelta$snp<-p0$snp
p0<-p0 %>% rename(p0=V3)

# check dataset
head(pdelta)
dim(pdelta)
head(p0)
dim(p0)

################################################################################
# %%R

# Make long to parse the column names
plong<-
  pdelta %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to=c("site","year","rep"),
    values_to = c("freq"),
    names_pattern = "X?(.*)_(.)_(.*)"
  ) %>%
  dplyr::select(site, year, rep, freq, snp)

df<- plong

# Extract the terminal year of frequency change
df<-
  df %>%
  group_by(site, rep, snp) %>%
  mutate(year=as.numeric(year)) %>%
  reframe(max_year = if(all(is.na(freq))) NA else max(year[!is.na(freq)], na.rm = TRUE),
          maxfreq = freq[year==max_year])

# Add starting frequency
df<- merge(
  df,
  select(p0,p0,snp)
)
# Create change in frequency
df<- df %>% 
  mutate(deltafreq= maxfreq - p0)

head(df)

allele_terminal_frequencies_long = df

################################################################################
write.csv( file= paste0(myfolder,"/intermediate-data/allele_terminal_frequencies_long.csv"), allele_terminal_frequencies_long)
save(file=paste0(myfolder,"/intermediate-data/allele_terminal_frequencies_long.rda"),allele_terminal_frequencies_long)
