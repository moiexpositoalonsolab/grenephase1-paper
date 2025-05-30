rm(list=ls())

library(dplyr)
library(geodata)
require(raster)
require(sp)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(ggthemes)

source("functions-utilities-map-highquality.R")


## This script use stabilizing selection model to predict the fitness loss from 2000 to 2030 

prefix <- "FUTURE_PREDICTION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)


future_current_fitness_ratio <- function(z_future,z_current,zo,vs){
  ratio <- exp(((z_current-zo)^2 - (z_future-zo)^2)/vs)
  return(ratio)
}

## first predict the 231 grenenet founder

grenet_prs <- read.delim("data-intermediate/CLIMATE_GWAS_bslmm_prs_ecotype_climate.txt")

ecotype_info <- read.delim("data/1001g_regmap_grenet_ecotype_info_corrected_2024May16.csv",sep=",") %>%
  filter(ecotype_id %in% grenet_prs$ecotype) %>%
  arrange(match(ecotype_id, grenet_prs$ecotype))

## future climate
long <- ecotype_info$Longitude_corrected
lat <- ecotype_info$Latitude_corrected

coords <- data.frame(x=long, y = lat)
mycrs <- "+proj=longlat +datum=WGS84"

wc <- raster::stack(x = "data-external/cmip6/2.5m/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp585_2041-2060.tif")

pt <- sp::SpatialPoints(coords = coords, proj4string = CRS(mycrs))
# extract from the given worldclim raster
wc_value <- raster::extract(x = wc, y = pt, method = 'simple', buffer=NULL)
wc_value[,1]

## collection climate
collection_climate <- read.delim("data-external/bioclimvars_ecotypes_era5.csv",sep=",") %>%
  filter(ecotype %in% grenet_prs$ecotype) %>%
  arrange(match(ecotype,grenet_prs$ecotype))

## ecotype stabilizing selection models
stabilizing_selection <- read.delim("data-intermediate/LOCAL_ADAPTATION_ecotype_Wmax_Vs.txt")

ratio <- c()

for(i in 1:nrow(stabilizing_selection)){
  #ratio <- c(ratio,future_current_fitness_ratio(z_future =wc_value[i,1],z_current = collection_climate$bio1_new[i],zo = grenet_prs$PRS_bio1[i],vs = stabilizing_selection$Vs[i]))
  ratio <- c(ratio,future_current_fitness_ratio(z_future =wc_value[i,1],z_current = collection_climate$bio1_new[i],zo = collection_climate$bio1_new[i],vs = stabilizing_selection$Vs[i]))
}

hist(ratio)

df <- tibble(ecotype = grenet_prs$ecotype,
                  lat = ecotype_info$Latitude_corrected,
                  long = ecotype_info$Longitude_corrected,
                  temperaturechange = (wc_value[,1] - collection_climate$bio1_new)^2,
                  fitness_loss = ratio)


########### draw the map #########

## Broad European range
xlim = c(-15, +90)
ylim = c(25, +65)
Range=extent(c(xlim,ylim))

## future climate map 
fut<-stack("data-external/cmip6/2.5m/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp585_2041-2060.tif")
names(fut) <- paste0("bio",1:19)
fut<- fut %>% crop(.,Range)

## current climate map
now<-stack(list.files("data-external/wc2.1_2.5/", pattern = "\\.tif$", full.names = TRUE))
names(now)<-paste0("bio",1:19)
now<- now %>% crop(.,Range)

temperaturechange<- (fut[[1]]-now[[1]])


source("functions-utilities-map-highquality.R")
basemap<- coolmaptile(temperaturechange, rev(c('black',rev(brewer.pal(name='Greys',4))) ),fillname = "Temp change ssp585 2050")
basemap

basemap+
  geom_point(data=df,aes(y=lat, x=long,color=fitness_loss))+
  scale_color_gradientn(colors = c(rev(RColorBrewer::brewer.pal(5,"Reds"))),
                        values = scales::rescale(c(0.8,0.95,0.98,0.99,1)))



################# draw the map for 1001g ##########

stabilizing_selection_1001g <- read.delim("data-intermediate/FUTURE_PREDICTION_1001g_grenet_predicted_Wmax_Vs.txt",sep=" ")

ecotype_info_1001g <- read.delim("data/1001g_regmap_grenet_ecotype_info_corrected_2024May16.csv",sep=",")  %>%
  filter(ecotype_id %in% stabilizing_selection_1001g$ecotype) %>%
  arrange(match(ecotype_id,stabilizing_selection_1001g$ecotype))

## future climate
long <- ecotype_info_1001g$Longitude_corrected
lat <- ecotype_info_1001g$Latitude_corrected

coords <- data.frame(x=long, y = lat)
mycrs <- "+proj=longlat +datum=WGS84"

wc <- raster::stack(x = "data-external/cmip6/2.5m/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp585_2041-2060.tif")

pt <- sp::SpatialPoints(coords = coords, proj4string = CRS(mycrs))
# extract from the given worldclim raster
wc_value <- raster::extract(x = wc[[1]], y = pt, method = 'simple', buffer=NULL,)
wc_value

## collection climate
collection_climate <- read.delim("data-external/bioclimvars_ecotypes_era5.csv",sep=",") %>%
  filter(ecotype %in% stabilizing_selection_1001g$ecotype) %>%
  arrange(match(ecotype,stabilizing_selection_1001g$ecotype))

ratio_1001g <- c()

for(i in 1:nrow(stabilizing_selection_1001g)){
  #ratio <- c(ratio,future_current_fitness_ratio(z_future =wc_value[i,1],z_current = collection_climate$bio1_new[i],zo = grenet_prs$PRS_bio1[i],vs = stabilizing_selection$Vs[i]))
  ratio_1001g <- c(ratio_1001g,future_current_fitness_ratio(z_future =wc_value[i],z_current = collection_climate$bio1_new[i],zo = collection_climate$bio1_new[i],vs = stabilizing_selection_1001g$Vs[i]))
}

hist(ratio_1001g)

df_1001g <- tibble(ecotype = stabilizing_selection_1001g$ecotype,
             lat = ecotype_info_1001g$Latitude_corrected,
             long = ecotype_info_1001g$Longitude_corrected,
             temperaturechange = (wc_value - collection_climate$bio1_new)^2,
             fitness_ratio = ratio_1001g)

########### draw the map #########

## future climate map 
fut<-stack("data-external/cmip6/2.5m/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp585_2041-2060.tif")
names(fut) <- paste0("bio",1:19)
#fut<- fut %>% crop(.,Range)

## current climate map
now<-stack(list.files("data-external/wc2.1_2.5/", pattern = "\\.tif$", full.names = TRUE))
names(now)<-paste0("bio",1:19)
#now<- now %>% crop(.,Range)

temperaturechange<- (fut[[1]]-now[[1]])


basemap<- coolmaptile(temperaturechange, rev(c('black',rev(brewer.pal(name='Greys',4))) ),fillname = "Temp change ssp585 2050")

future_1001g_ssp585_2050 <- basemap+
  geom_point(data=df_1001g,aes(y=lat, x=long,color=fitness_loss))+
  scale_color_gradientn(colors = c(rev(RColorBrewer::brewer.pal(5,"Reds"))),
                        values = scales::rescale(c(0.8,0.95,0.98,0.99,1.02)))
ggsave("figs/FUTURE_PREDICTION_1001g_2050_ssp585.pdf",plot = future_1001g_ssp585_2050,height = 8,width = 10,device = "pdf",units = "in")
ggsave("figs/FUTURE_PREDICTION_1001g_2050_ssp585.png",plot = future_1001g_ssp585_2050,height = 8,width = 10,device = "png",units = "in",dpi = 800)
