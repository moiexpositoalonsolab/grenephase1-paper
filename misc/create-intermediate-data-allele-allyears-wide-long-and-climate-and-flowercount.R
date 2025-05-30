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
# 
# file_path_p0 <- paste0(myfolder, file_p0)
# file_path_pdelta <- paste0(myfolder, file_pdelta)
# file_path_edelta <- paste0(myfolder, file_edelta)
# file_path_e0 <- paste0(myfolder, file_e0)
# file_path_eid <- paste0(myfolder, file_eid)
# file_path_climate_sites=paste0(myfolder, file_climate_sites)
# file_path_sites_info=paste0(myfolder, file_sites_info)

################################################################################

# Load frequencies 
# pdelta<-read.table("data/merged_hapFIRE_delta_p_LDpruned.txt", header=T)
p<-read.table("data/merged_hapFIRE_allele_frequency_LDpruned.txt", header=T)
p0<-read.table("data/average_seedmix_p0_LDpruned.txt", header=F)
p0$snp<-paste0(p0$V1,"_",p0$V2)
p$snp<-p0$snp
p0<-p0 %>% dplyr::rename(startfreq=V3,)


################################################################################
# Make a long dataset with all years
# Make long to parse the column names
plong<-
  p %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to=c("site","year","rep"),
    values_to = c("freq"),
    names_pattern = "X?(.*)_(.)_(.*)"
  ) %>%
  dplyr::select(site, year, rep, freq, snp)

df<- plong
# Add starting frequency
df<- merge(
  df,
  dplyr::select(p0,startfreq,snp)
)
# Save **
allele_allyears_frequencies_long_raw = df 
message("saving alleles all years frequency, long format, raw")
write.csv( file= paste0(myfolder,"/data-intermediate/allele_allyears_frequencies_long_raw.csv"), 
           allele_allyears_frequencies_long_raw)
save(file=paste0(myfolder,"/data-intermediate/allele_allyears_frequencies_long_raw.rda"),
          allele_allyears_frequencies_long_raw)

# Make a long dataset with the terminal year
# Extract the terminal year of frequency change
df<-
  df %>%
  group_by(site, rep, snp) %>%
  mutate(year=as.numeric(year)) %>%
  reframe(max_year = if(all(is.na(freq))) NA else max(year[!is.na(freq)], na.rm = TRUE),
          maxfreq = freq[year==max_year])
# # Create change in frequency
# df<- df %>% 
#   mutate(deltafreq= maxfreq - p0)


# Save **
allele_terminal_frequencies_long_raw = df 
message("saving alleles terminal year frequency, long format, raw")
write.csv( file= paste0(myfolder,"/data-intermediate/allele_terminal_frequencies_long_raw.csv"), 
           allele_terminal_frequencies_long_raw)
save(file=paste0(myfolder,"/data-intermediate/allele_terminal_frequencies_long_raw.rda"),
           allele_terminal_frequencies_long_raw)

################################################################################
# Add climate
# Load climates
file_path_climate_sites=paste0(myfolder, "/grene/data/worldclim_sitesdata.rda")
file_path_sites_info=paste0(myfolder, "/grene/data/locations_data.rda")
file_path_climate_ecotypes=paste0(myfolder, "/grene/data/worldclim_ecotypesdata.rda")

load(file_path_climate_sites)
load(file_path_sites_info)
load(file_path_climate_ecotypes)

colnames(worldclim_sitesdata)[-1]<-paste0("sites_",colnames(worldclim_sitesdata)[-1])
colnames(worldclim_ecotypesdata)[-1] <- paste0("ecotypes_",colnames(worldclim_ecotypesdata)[-1])
locations_data<-locations_data %>% 
  dplyr::select(site,longitude,latitude)

# Population
# Add population
# load("data/flowers_survival_long.rda")
# flowers_survival_long
load("grene/data/samples_data.rda")
samples_tmp<-samples_data %>% 
  dplyr::filter(year<2021) %>% 
  dplyr::mutate(year= year-2017) %>% 
  dplyr::rename(flowers=flowerscollected) %>% 
  dplyr::select(site,plot,year, date,month,day,flowers) %>% 
  dplyr::mutate(julian_day = yday(ymd(date)))


################################################################################
### Climate and popsize for terminal years

# Climate
df<-allele_terminal_frequencies_long_raw
a<-df
a<-
  a %>% 
  merge(
    .,
    locations_data, by='site'
  ) %>% 
  merge(
    .,
    worldclim_sitesdata, by="site"
  )

## Save
message("saving alleles all years frequency, long format, delta, with climate")
allele_terminal_frequencies_long_raw_climate = a
a[1:5,1:5]
# write.csv( file= paste0(myfolder,"/data-intermediate/allele_terminal_frequencies_long_raw_climate.csv"), 
#            allele_terminal_frequencies_long_raw_climate)
save(file=paste0(myfolder,"/data-intermediate/allele_terminal_frequencies_long_raw_climate.rda"),
     allele_terminal_frequencies_long_raw_climate)

# Add population
df<-allele_terminal_frequencies_long_raw
a2<-df

a2<-
  a2 %>% 
  merge(
    .,
    locations_data, by='site'
  ) %>% 
  merge(
    .,
    samples_tmp, by.x=c("site","max_year","rep"), by.y=c("site","year","plot")
  )

# Reduce the size to only bio1-19
# Trick to merge worldclim data
a3<-a2
rownames(worldclim_sitesdata)<-worldclim_sitesdata$site
subsetclim <- worldclim_sitesdata[a3$site, 2:20]
a3<-cbind(a3,subsetclim)


message("saving alleles all years frequency, long format, delta, with climate and flower collection size")
allele_terminal_frequencies_long_raw_climate_flowers = a3
a3[1:5,1:15]
a3 %>% head
# write.csv( file= paste0(myfolder,"/data-intermediate/allele_terminal_frequencies_long_raw_climate_flowers.csv"), 
#            allele_terminal_frequencies_long_raw_climate_flowers)
save(file=paste0(myfolder,"/data-intermediate/allele_terminal_frequencies_long_raw_climate_flowers.rda"),
     allele_terminal_frequencies_long_raw_climate_flowers)


################################################################################
### Climate and popsize for all years and frequencies

### Merge with all years
df<-allele_allyears_frequencies_long_raw
a<-df
a<-
  a %>% 
  merge(
    .,
    locations_data, by='site'
  ) %>% 
  merge(
    .,
    worldclim_sitesdata, by="site"
  )

## Save
message("saving alleles all years frequency, long format, delta, with climate")
allele_allyears_frequencies_long_raw_climate = a
a[1:5,1:5]
# write.csv( file= paste0(myfolder,"/data-intermediate/allele_allyears_frequencies_long_raw_climate.csv"), 
#            allele_allyears_frequencies_long_raw_climate)
save(file=paste0(myfolder,"/data-intermediate/allele_allyears_frequencies_long_raw_climate.rda"),
           allele_allyears_frequencies_long_raw_climate)
################################################################################
# Add population
# load("data/flowers_survival_long.rda")
# flowers_survival_long
load("grene/data/samples_data.rda")
samples_tmp<-samples_data %>% 
  dplyr::filter(year<2021) %>% 
  dplyr::mutate(year= year-2017) %>% 
  dplyr::rename(flowers=flowerscollected) %>% 
  dplyr::select(site,plot,year, date,month,day,flowers) %>% 
  dplyr::mutate(julian_day = yday(ymd(date)))

# Reduce the size to only bio1-19
df<-allele_allyears_frequencies_long_raw
a2<-df

a2<-
  a2 %>% 
  merge(
    .,
    locations_data, by='site'
  ) %>% 
  merge(
    .,
    samples_tmp, by.x=c("site","year","rep"), by.y=c("site","year","plot")
  )

# Trick to merge worldclim data
a3<-a2
rownames(worldclim_sitesdata)<-worldclim_sitesdata$site
subsetclim <- worldclim_sitesdata[a3$site, 2:20]
a3<-cbind(a3,subsetclim)
  

message("saving alleles all years frequency, long format, delta, with climate and flower collection size")
allele_allyears_frequencies_long_raw_climate_flowers = a3
a3[1:5,1:15]
a3 %>% head
# write.csv( file= paste0(myfolder,"/data-intermediate/allele_allyears_frequencies_long_raw_climate_flowers.csv"), 
#            allele_allyears_frequencies_long_raw_climate_flowers)
save(file=paste0(myfolder,"/data-intermediate/allele_allyears_frequencies_long_raw_climate_flowers.rda"),
     allele_allyears_frequencies_long_raw_climate_flowers)
