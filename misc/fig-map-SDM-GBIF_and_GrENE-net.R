################################################################################
### Goal
### Get the species distribution of Arabidopsis
### Make something similar but with GRENE net
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
## Worldwide. maps may be nice but we also can do zoom in, here are some good
# ranges for Europe

# Broad European range
xlim = c(-10.5, +90)
ylim = c(25, +65)
Range=extent(c(xlim,ylim))
# Small range
xlim=c(-10.5,+ 53)
ylim=c(32,65)
Range=extent(c(xlim,ylim))


################################################################################
# Load climates

# now<-rbioclim::getData(name="worldclim2",var='bio', res=2.5,path = "data-external/")
# The present data is already downloaded from worldclim 2
now<-stack(list.files("data-external/wc2.1_2.5/", pattern = "\\.tif$", full.names = TRUE))
names(now)<-paste0("bio",1:19)
now<- now %>% crop(.,Range)


# fut<-getData(name="CMIP5",var='bio', res=5,model='MP', year=50, rcp=85)
# Same, future data is already downloaded from worldclim 2
fut<-stack("data-external/cmip6/2.5m/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2021-2040.tif")
names(fut)<-paste0("bio",1:19)
fut<- fut %>% crop(.,Range)

# Save future climate data for point location of ecotypes and sites
load("grene/data/worldclim_ecotypesdata.rda")
load("grene/data/ecotypes_data.rda")
load("grene/data/sites.clim.rda")
worldclim_ecotypesdata %>% head
ecotypes_data %>% head
sites.clim %>% head

ecotypes_future_climate<-extract(fut,cbind(ecotypes_data$longitude,ecotypes_data$latitude)) %>% data.frame
ecotypes_future_climate$id<-ecotypes_data$ecotypeid
sites_future_climate<-extract(fut,cbind(sites.clim$LONGITUDE,sites.clim$LATITUDE))%>% data.frame
sites_future_climate$site<-sites.clim$SITE_CODE
save(file = "data-intermediate/ecotypes_future_climate.rda",ecotypes_future_climate)
save(file = "data-intermediate/sites_future_climate.rda",sites_future_climate)
  
################################################################################
# Plot SDM
# The SDM was fited as this: https://github.com/moiexpositoalonsolab/arabidopsisrange
predictagain=F
if(predictagain){
# To make a prediction, you can load the fitted model below and provide the present climate
  library(dismo)
  library(raster)
  library(rJava)
  load("data-external/mod.rda")
  nowtmp<-now
  names(nowtmp)<-c("wc2.1_2.5m_bio_1","wc2.1_2.5m_bio_2","wc2.1_2.5m_bio_3","wc2.1_2.5m_bio_4","wc2.1_2.5m_bio_5","wc2.1_2.5m_bio_6","wc2.1_2.5m_bio_7","wc2.1_2.5m_bio_8","wc2.1_2.5m_bio_9","wc2.1_2.5m_bio_10","wc2.1_2.5m_bio_11","wc2.1_2.5m_bio_12","wc2.1_2.5m_bio_13","wc2.1_2.5m_bio_14","wc2.1_2.5m_bio_15","wc2.1_2.5m_bio_16","wc2.1_2.5m_bio_17","wc2.1_2.5m_bio_18","wc2.1_2.5m_bio_19")
  # sdm=predict(mod,now)
  sdm=dismo::predict(object=mod, x=nowtmp)
}else{
  load("data-external/sdm.rda")
  sdm <- sdm %>% crop(.,Range)
}
sdmcolor=brewer.pal(9,"Greys")[-1]

# pdf("figs/fig-map-SDM.pdf",height = 8,width = 9)
# plot(sdm, col=sdmcolor)
# dev.off()
# 
# png("figs/fig-map-SDM.png",height = 8,width = 9)
# plot(sdm, col=sdmcolor)
# dev.off()

# Plot the SDM with a high quality function
source("functions-utilities-map-highquality.R")
basemap<- coolmaptile(sdm)
basemap

################################################################################
# Test predictability with the actual number of flowers collected
# flowers <- 
read.csv("data/merged_sample_table.csv") %>% head

load("grene/data/samples_data.rda")
samples_data %>% head

################################################################################
# Generate GRENE_net SDM
# Test predictability with the estimated number of flowers

## Survival data
sur=read.csv("data/survival.csv")

## Survival predicted
popsize=read.csv("pop_size_estimation/pop_size_estimations.csv")
popsizeraw=popsize

## Cliamte 
load("grene/data/worldclim_sitesdata.rda")
load("grene/data/locations_data.rda")
## Merge survival and climate and locations

sur=merge(sur,worldclim_sitesdata,by='site')
sur=merge(sur,locations_data,by='site')

# Get a summary of survival
sursummary<-
  sur %>% 
  group_by(site,latitude,longitude) %>% 
  filter(X2_survival != "-1") %>% 
  summarise(survival_proportion=mean(X2_survival))

###### 
## Merge 
popsize=merge(popsize,worldclim_sitesdata,by='site')
popsize=merge(popsize,locations_data,by='site')

#Get a summary
popsizesummary<-
  popsize %>% 
  group_by(site,latitude,longitude) %>% 
  summarise(
    avg_flowers=sum(flowerscollected_corrected),
    avg_popsize=sum(totalplantnumber_complete)
  )

################################################################################
# Fit GLMs or Random Forest to census data

##########################################
## Try to predict number of flowers
tofit<-
popsize %>% 
  dplyr::select(flowerscollected_corrected, starts_with("bio"), generation) %>% 
  dplyr::rename(flowers=flowerscollected_corrected) 

# LM
library(stats)
lm(data=tofit,flowers~.)->
  flowermod
summary(flowermod)

tofit$flowers_predicted
lm(tofit$flowers~predict(flowermod)) %>% summary
plot(tofit$flowers~predict(flowermod)) 



# RF
library(randomForest)
randomForest(data=tofit,flowers~.)->
  flowerrf
flowerrf
plot(tofit$flowers~predict(flowerrf)) 
lm(tofit$flowers~predict(flowerrf)) %>% summary

# Project number of flowers in a map
flowerpred<-predict(now,flowermod)
flowerpred<-predict(now,flowerrf)
plot(flowerpred, col=sdmcolor)

##########################################
## Try to predict total number of plants from climate
tofit<-
  popsize %>% 
  dplyr::select(totalplantnumber_complete, generation,starts_with("bio")) %>% 
  dplyr::rename(plants=totalplantnumber_complete) 

# LM
lm(data=tofit,plants~.)->
  plantmod
summary(plantmod)

lm(tofit$plants~predict(plantmod)) %>% summary
plot(tofit$plants~predict(plantmod)) 

# RF
randomForest(data=tofit,plants~.)->
  plantrf
plantrf

# Project number of flowers in a map
plantpred<-predict(now,plantrf)
plot(plantpred, col=sdmcolor)

plantpredglm<-predict(now,plantmod)
plot(plantpredglm, col=sdmcolor)

# What about a GLM that is polynomial?
## Did not improve

# tofit<-
#   popsize %>% 
#   dplyr::select(totalplantnumber_complete, starts_with("bio")) %>% 
#   dplyr::rename(plants=totalplantnumber_complete) 
# # generate square variables
# tofitsquare<-tofit[,-1]^2
# colnames(tofitsquare)<-paste0("bio",1:19,"22")
# tofit<-cbind(tofit,tofitsquare)
# # also for the map
# nowsquare<-now
# values(nowsquare)<-values(nowsquare)^2
# names(nowsquare)<-paste0("bio",1:19,"22")
# nowsquare<-stack(now,nowsquare)
# 
# # # fit polynomial
# tofit %>% 
#   dplyr::select(plants,bio1,bio12,bio122,bio1222) %>% 
# lm(data=.,plants~.)->
#   plantmodsq
# summary(plantmodsq)
# plot(predict(nowsquare,plantmodsq), col=sdmcolor)

##########################################
## Try to predict survival per replicate in the 3rd year

tofit<-
  sur %>% 
  dplyr::select(X3_survival, starts_with("bio")) %>% 
  dplyr::rename(survival=X3_survival) %>% 
  dplyr::filter(survival>= 0) %>% 
  na.omit()
dim(tofit)

# Fit survival
library(randomForest)
randomForest(data=tofit,factor(survival)~.)->
  survivalrf
survivalrf
ggplot(data.frame(observed=tofit$survival, predicted=predict(survivalrf)))+
         geom_jitter(aes(x=predicted,y=observed),width = 0.1, height = 0.1)

# Fit survival per replicate with GLM
glm(data=tofit,factor(survival)~., family = "binomial")->
  survivalglm
survivalglm
survivalglm %>% summary()
survivalglm %>% R2

lm(tofit$survival~predict.glm(survivalglm)) %>% summary
plot(tofit$survival~predict.glm(survivalglm)) 

# plot survival in a map
survivalpredglm<-predict(now,survivalglm)
plot(survivalpredglm, col=sdmcolor)


##########################################
## Try to predict survival per replicate in the 3rd year
