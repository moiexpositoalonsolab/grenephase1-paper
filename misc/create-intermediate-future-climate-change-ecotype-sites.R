################################################################################
### Goal
### Extract future climate data for sites and ecotypes

################################################################################
################################################################################
library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(cowplot)

################################################################################
#### Location
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
# myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)


################################################################################
# Get ecotype data

load("grene/data/worldclim_ecotypesdata.rda")
load("grene/data/ecotypes_data.rda")

load("grene/data/locations_data.rda")
load("grene/data/worldclim_sitesdata.rda")

tmp= merge(ecotypes_data, by.x="ecotypeid", worldclim_ecotypesdata, by.y="ecotypeid")

tpm2 = merge(dplyr::select(locations_data, site, latitude,longitude),
             worldclim_sitesdata, by="site"
             )

################################################################################
## Worldwide. maps may be nice but we also can do zoom in, here are some good
# ranges for Europe

# Broad European range
xlim = c(-15, +90)
ylim = c(25, +65)
Range=extent(c(xlim,ylim))
# # Small range
# xlim=c(-10.5,+ 53)
# ylim=c(32,65)
# Range=extent(c(xlim,ylim))


################################################################################
# Load climates

# now<-rbioclim::getData(name="worldclim2",var='bio', res=2.5,path = "data-external/")
# The present data is already downloaded from worldclim 2
now<-stack(list.files("data-external/wc2.1_2.5/", pattern = "\\.tif$", full.names = TRUE))
names(now)<-paste0("bio",1:19)
now<- now %>% crop(.,Range)

# Download using rbioclim
# fut<-stack("data-external/cmip6_2.5/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2021-2040.tif")
library(rbioclim)
library(raster)
var="bio"
model="INM-CM5-0"
model="MPI-ESM1-2-HR"
res="2.5"
year="2021-2040"
rcp=245
downloadmethod="wget"
# path='/Users/moisesexpositoalonso'
# fut<-rbioclim::getData(name="CMIP6",var=var,model=model,rcp=rcp,year=year,res=res,path = path,download = TRUE, downloadmethod="auto")
path=paste0(myfolder, 'data-external/')
fut<-rbioclim::getData(name="CMIP6",var=var,model=model,rcp=rcp,year=year,res=res,path = path,download = TRUE, downloadmethod="auto")
names(fut)<-paste0("bio",1:19)
fut<- fut %>% crop(.,Range)

year="2041-2060"
rcp=585
downloadmethod="wget"
# path='/Users/moisesexpositoalonso'
# fut<-rbioclim::getData(name="CMIP6",var=var,model=model,rcp=rcp,year=year,res=res,path = path,download = TRUE, downloadmethod="auto")
path=paste0(myfolder, 'data-external/')
futbad<-rbioclim::getData(name="CMIP6",var=var,model=model,rcp=rcp,year=year,res=res,path = path,download = TRUE, downloadmethod="auto")
names(futbad)<-paste0("bio",1:19)
futbad<- futbad %>% crop(.,Range)

# The most simple is to do the distance square
temperaturechange<- (fut[[1]]-now[[1]])^2
temperaturechangebad<- (futbad[[1]]-now[[1]])^2

# Extract the future trajectory for each ecotype
tmp$temperaturechange<-extract(temperaturechange,cbind(tmp$longitude,tmp$latitude))
# Extract future trajectory bad scenario
tmp$temperaturechangebad<-extract(temperaturechangebad,cbind(tmp$longitude,tmp$latitude))

# Extract the data for ecotypes
rcp245<-extract(fut,cbind(ecotypes_data$longitude,ecotypes_data$latitude))
colnames(rcp245)<-paste0("ecotypes_",colnames(rcp245),"_rcp245")
rcp585<-extract(futbad,cbind(ecotypes_data$longitude,ecotypes_data$latitude))
colnames(rcp585)<-paste0("ecotypes_",colnames(rcp585),"_rcp585")

# Extract the data for sites

sites_rcp245<-extract(fut,cbind(locations_data$longitude,locations_data$latitude))
colnames(sites_rcp245)<-paste0("sites_",colnames(sites_rcp245),"_rcp245")
sites_rcp585<-extract(futbad,cbind(locations_data$longitude,locations_data$latitude))
colnames(sites_rcp585)<-paste0("sites_",colnames(sites_rcp585),"_rcp585")

### SAVE INTERMEDIATE
sitesdata<-cbind(sites_rcp245,sites_rcp585)%>% data.frame
sitesdata$site<-locations_data$site
write.csv( file= paste0(myfolder,"/data-intermediate/future_wordlclim_sites.csv"), 
          sitesdata           
           )
ecotypesdata<-cbind(rcp245,rcp585) %>% data.frame
ecotypesdata$id<-ecotypes_data$ecotypeid
write.csv( file= paste0(myfolder,"/data-intermediate/future_wordlclim_ecotypes.csv"), 
           ecotypesdata
            )


################################################################################
# Create ecotypes by sites combinations
cc<-
    expand.grid(ecotypes_data$ecotypeid, locations_data$site) %>% 
    dplyr::rename(id=Var1, site=Var2) %>% 
    merge(.,sitesdata,by="site",all.x=T) %>% 
    merge(.,ecotypesdata,by="id",all.x=T)

# Create the equivalent of the present
climsite<-
  worldclim_sitesdata %>% 
  dplyr::select(site,starts_with("bio"))
  colnames(climsite)[-1]<-paste0("sites_",colnames(climsite)[-1])
climeco<-
  worldclim_ecotypesdata %>% 
    dplyr::select(ecotypeid,starts_with("bio"))
  colnames(climeco)[-1]<-paste0("ecotypes_",colnames(climeco)[-1])
  
clim<-
  expand.grid(ecotypes_data$ecotypeid, locations_data$site) %>% 
  dplyr::rename(id=Var1, site=Var2) %>% 
  merge(.,climsite,by="site",all.x=T) %>% 
  merge(.,climeco,by.x="id",by.y="ecotypeid",all.x=T)
  
dim(cc)
dim(clim)

ccc<-cbind(
    clim, 
    dplyr::select(cc,-site,-id))
  
future_climate_sites_ecotypes<-ccc
save(file = "data-intermediate/future_climate_sites_ecotypes.rda",future_climate_sites_ecotypes)
write.csv(file = "data-intermediate/future_climate_sites_ecotypes.csv",future_climate_sites_ecotypes, row.names = F)


# ################################################################################
# # Merge with ecotyeps
# load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize_flowers.rda")
# load("data-intermediate/ecotype_terminal_frequencies_long_raw_climate.rda")
# 
# e= ecotype_allyears_frequencies_long_raw_climate_popsize_flowers
# e= ecotype_terminal_frequencies_long_raw_climate %>% 
#     dplyr::select(id,site,rep,max_year, startfreq,maxfreq,
#                   starts_with("sites_bio"),starts_with("ecotypes_bio"))
# 
# futecotypes<-read.csv("data-intermediate/future_wordlclim_ecotypes.csv")
# futsites<-read.csv("data-intermediate/future_wordlclim_sites.csv")
# 
# emerge<-
#   e %>% 
#   merge(.,futsites,by="site") %>% 
#   merge(.,futecotypes,by="id")
# 
# ecotype_allyears_frequencies_long_raw_climate_popsize_flowers_futureclimate<-emerge
# save("data-intermediate/ecotype_allyears_frequencies_long_raw_climate_futureclimate.rda",ecotype_allyears_frequencies_long_raw_climate_popsize_flowers_futureclimate)
# save("data-intermediate/ecotype_allyears_frequencies_long_raw_climate_futureclimate.csv",emerge)
