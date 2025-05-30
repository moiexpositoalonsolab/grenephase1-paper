################################################################################
### Goal
### Estimate some basic metrics of ecotypes increasing in frequency across
### sites, means, variances, optimal places. Including poulation size

################################################################################
library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(patchwork)

################################################################################
#### Location
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
# myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
# Load the long format with ecotype frequencies and climates

load("grene/data/sites.clim.rda")
load("grene/data/ecotypes_data.rda")

# load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate.rda")
# load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize.rda")

load(file=paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize.rda"))
load(file=paste0(myfolder,"/data-intermediate/ecotype_terminal_frequencies_long_raw_climate_popsize.rda"))
e<-enew <- ecotype_allyears_frequencies_long_raw_climate_popsize 
# e<-enew <- ecotype_terminal_frequencies_long_raw_climate_popsize %>% 
#   rename(year=max_year, freq=maxfreq) # Use the terminal
# enew[1:5,1:10]

#####*********************************************************************######
# Merge with location of origin
em<-
  e %>% 
  dplyr::select(-latitude,-longitude) %>% 
  merge(.,by.x="site",
        rename(dplyr::select(sites.clim,SITE_CODE,LONGITUDE,LATITUDE),
               sites_latitude=LATITUDE,
               sites_longitude=LONGITUDE),
        by.y="SITE_CODE"
        ) %>% 
  merge(.,by.x="id",
        rename(dplyr::select(ecotypes_data,ecotypeid,longitude,latitude), 
                ecotypes_latitude=latitude,
                ecotypes_longitude=longitude),
        by.y="ecotypeid"
        )

# Plot distance  
ggplot(em)+
  geom_point(aes(y=freq, x=sites_latitude-ecotypes_latitude))+
  geom_vline(xintercept = 0, lty='dotted')+
  geom_vline(xintercept = mean(em$sites_latitude-em$ecotypes_latitude), lty='dashed')+
  geom_vline(xintercept = weighted.mean(em$sites_latitude-em$ecotypes_latitude, w = em$freq), lty='solid')+
  xlab("Latitude transplant (site - ecotype)")+
  ylab("Change in ecotype freq")+
  theme_minimal()


ggplot(em)+
  geom_point(aes(y=freq, x=sites_longitude-ecotypes_longitude))+
  geom_vline(xintercept = 0, lty='dotted')+
  # geom_vline(xintercept = mean(em$sites_latitude-em$ecotypes_latitude), lty='dashed')+
  # geom_vline(xintercept = weighted.mean(em$sites_latitude-em$ecotypes_latitude, w = em$freq), lty='solid')+
  xlab("Latitude transplant (site - ecotype)")+
  ylab("Change in ecotype freq")+
  theme_minimal()


latlon_stabilizing<-
  em %>% 
  arrange((freq-startfreq)) %>% 
  ggplot(. ) +
  geom_point(aes(x=sites_longitude-ecotypes_longitude,
                 y=sites_latitude-ecotypes_latitude,
                 color=(freq-startfreq), size=(freq-startfreq)),
             alpha=0.9, shape=16)+
  scale_color_gradientn("freq. change (p1-p0)",colors=brewer.pal(9,"Greens")[-1])+
  scale_size_continuous("freq. change (p1-p0)", range = c(1,4))+
  geom_hline(yintercept = 0,lty='dotted')+
  geom_vline(xintercept = 0,lty='dotted')+
  theme_minimal()+
  labs(x = "Longitude distance (°)", y = "Latitude distance (°)")
latlon_stabilizing
save_plot("figs/fig_stabilizing_D2_latitude-longitude.pdf",latlon_stabilizing, base_width = 6,base_height =5)
save_plot("figs/fig_stabilizing_D2_latitude-longitude.png",latlon_stabilizing, base_width = 6,base_height =5)

# quantify how mnay samples
em %>% 
  dplyr::filter(sites_longitude> -15) %>% 
  dplyr::select(site) %>% unlist %>% unique %>% length
em %>% 
  dplyr::filter(sites_longitude> -15) %>% 
  nrow
em %>% 
  # dplyr::filter(sites_longitude> -15) %>% 
  dplyr::select(site) %>% unlist %>% unique %>% length
em %>% 
  # dplyr::filter(sites_longitude> -15) %>% 
  nrow



latlon_stabilizing_eu<-
  em %>% 
  dplyr::filter(sites_longitude> -15) %>% 
  arrange((freq-startfreq)) %>% 
  ggplot(. ) +
  geom_point(aes(x=sites_longitude-ecotypes_longitude,
                 y=sites_latitude-ecotypes_latitude,
                 color=(freq-startfreq), size=(freq-startfreq)),
             alpha=0.9, shape=16)+
  scale_color_gradientn("freq. change (p1-p0)",colors=brewer.pal(9,"Greens")[-1])+
  scale_size_continuous("freq. change (p1-p0)", range = c(1,4))+
  geom_hline(yintercept = 0,lty='dotted')+
  geom_vline(xintercept = 0,lty='dotted')+
  theme_minimal()+
  labs(x = "Longitude distance (°)", y = "Latitude distance (°)")  
latlon_stabilizing_eu
save_plot("figs/fig_stabilizing_D2_latitude-longitude-onlyEU.pdf",latlon_stabilizing_eu, base_width = 6,base_height =5)
save_plot("figs/fig_stabilizing_D2_latitude-longitude-onlyEU.pdf",latlon_stabilizing_eu, base_width = 6,base_height =5)

