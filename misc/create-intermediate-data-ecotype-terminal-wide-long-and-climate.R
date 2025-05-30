################################################################################
### Goal
### Create an ecotype frequency in long format for all sites
### with ecotype and site climate metadata
################################################################################
library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)

################################################################################
#### Location
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)


################################################################################
#### Files
file_path_climate_sites=paste0(myfolder, "/grene/data/worldclim_sitesdata.rda")
file_path_sites_info=paste0(myfolder, "/grene/data/locations_data.rda")
file_path_climate_ecotypes=paste0(myfolder, "/grene/data/worldclim_ecotypesdata.rda")

# Intermediate files
load("data-intermediate/ecotype_terminal_frequencies_long.rda")
load("data-intermediate/ecotype_terminal_frequencies_long_raw.rda")

e<-ecotype_terminal_frequencies_long_raw
e %>% head
################################################################################

###############
# Load climates
load(file_path_climate_sites)
load(file_path_sites_info)
load(file_path_climate_ecotypes)

colnames(worldclim_sitesdata)[-1]<-paste0("sites_",colnames(worldclim_sitesdata)[-1])
colnames(worldclim_ecotypesdata)[-1] <- paste0("ecotypes_",colnames(worldclim_ecotypesdata)[-1])
locations_data<-locations_data %>% 
  dplyr::select(site,longitude,latitude)

### Merge

e<-
e %>% 
  merge(
    .,
    locations_data, by='site'
  ) %>% 
  merge(
    .,
    worldclim_sitesdata, by="site"
  ) %>% 
  merge(
    .,
    worldclim_ecotypesdata, by.x="id", by.y="ecotypeid"
  )

ecotype_terminal_frequencies_long_raw_climate = e
write.csv( file= paste0(myfolder,"/data-intermediate/ecotype_terminal_frequencies_long_raw_climate.csv"), 
           ecotype_terminal_frequencies_long_raw_climate)
save(file=paste0(myfolder,"/data-intermediate/ecotype_terminal_frequencies_long_raw_climate.rda"),
     ecotype_terminal_frequencies_long_raw_climate)

