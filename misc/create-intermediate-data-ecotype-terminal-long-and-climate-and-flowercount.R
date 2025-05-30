################################################################################
### Goal
### Create an ecotype frequency in long format for all sites
### with climates but also flowers
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
load("data-intermediate/ecotype_terminal_frequencies_long_raw_climate.rda")
# load("data-intermediate/ecotype_terminal_frequencies_long_raw.rda")

e<-ecotype_terminal_frequencies_long_raw_climate
e %>% head
e[1:5,1:10]
################################################################################
#### Include the number of flowers

# This is the population size created by Tati
popsize=read.csv(paste0(myfolder, "/data-intermediate//pop_size_estimations.csv"))
popsize<-data.frame(popsize)
popsize<-popsize %>% dplyr::select(-X,-longitude,-latitude, -predicted)
popsize[1:5,1:8]

# Merge with ecotype frequency
epop=merge(
  x=e, by.x=c("site","rep","max_year"),
  y=popsize, by.y=c("site","plot","generation")
)
head(epop)
epop[1:5,1:10]
colnames(epop)
dim(epop)


ecotype_terminal_frequencies_long_raw_climate_popsize = epop
write.csv( file= paste0(myfolder,"/data-intermediate/ecotype_terminal_frequencies_long_raw_climate_popsize.csv"), 
           ecotype_terminal_frequencies_long_raw_climate_popsize)
save(file=paste0(myfolder,"/data-intermediate/ecotype_terminal_frequencies_long_raw_climate_popsize.rda"),
     ecotype_terminal_frequencies_long_raw_climate_popsize)

