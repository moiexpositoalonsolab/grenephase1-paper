rm(list=ls())

library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)


myfolder<-"~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)


################################################################################
# Figure out map first
# Broad European range
xlim = c(-15, +90)
ylim = c(25, +65)
Range=extent(c(xlim,ylim))
# Small range
xlim=c(-10.5,+ 53)
ylim=c(32,65)
# Range=extent(c(xlim,ylim))
# Load the SDM for 
load("data-external/sdm.rda")
sdm <- sdm %>% crop(.,Range)

sdmcolor=brewer.pal(9,"Greys")[-1]

source("functions-utilities-map-highquality.R")
basemap<- coolmaptile(sdm)
basemap


################################################################################

## load Vs Wmax 
Vsresults<-read.table("data-intermediate/LOCAL_ADAPTATION_ecotype_Wmax_Vs.txt",header = T)

## load ecotype_info
ecotype_climate <- read.delim("data-external/bioclimvars_ecotypes_era5.csv",sep=",")
colnames(ecotype_climate) <- c("ecotype",paste0("bio",1:19,"_ecotype"))

ecotype_info <- read.delim("data/1001g_regmap_grenet_ecotype_info_corrected_2024May16.csv",sep = ",")
colnames(ecotype_info)[c(1,5,6)] <- c("ecotype","latitude_ecotype","longitude_ecotype")

ecotype_info <- ecotype_info %>%
  filter(GrENE.net == TRUE) %>%
  arrange(match(ecotype,Vsresults$ecotype)) %>%
  left_join(.,Vsresults,by="ecotype") %>%
  left_join(.,ecotype_climate,by="ecotype") %>%
  mutate(inv_Vs = 1/Vs)


quantiles <- quantile(ecotype_info$inv_Vs, probs = seq(0, 1, 0.1))
ecotype_info$quantile_group <- cut(ecotype_info$inv_Vs, breaks = quantiles, include.lowest = TRUE, labels = 1:10)

# PLOT VS on a map

# Make the geographic plot
map_and_vs<-
  basemap+
  geom_point(data=ecotype_info, aes(x=longitude_ecotype, y=latitude_ecotype, color=inv_Vs),cex=3,shape=16)+
  scale_color_gradientn(colors=c(rev(RColorBrewer::brewer.pal(5,"Reds")),RColorBrewer::brewer.pal(5,"Greens")),
                        values = rescale(c(0,0.001,0.002,0.004,0.006,0.008,0.011,0.015,0.02)))
map_and_vs

ggsave(plot = map_and_vs,filename = "figs/LOCAL_ADAPTATION_inv_Vs_map.pdf",device = "png", height = 7,width = 9,dpi = 800,units = "in")
#ggsave(plot = map_and_vs,filename = "figs/LOCAL_ADAPTATION_Vs_map.pdf",device = "pdf", height = 5,width = 7)


# Make the plot against latitude
fig_vs_and_latitude<-
  ggplot(data=ecotype_info)+
  # geom_point(aes(x=latitude, y=1/Vs),shape=16,color="black")+
  geom_point(aes(x=latitude_ecotype, y=inv_Vs, color=inv_Vs),shape=16,cex=3)+
  scale_color_gradientn(colors=c(rev(RColorBrewer::brewer.pal(5,"Reds")),RColorBrewer::brewer.pal(5,"Greens")),
                        values = rescale(c(0,0.001,0.002,0.004,0.006,0.008,0.011,0.015,0.02)))+
  stat_smooth(aes(x=latitude_ecotype, y=inv_Vs, color=inv_Vs),method = "glm", color='grey')+
  annotate(geom = "text",x=25,y=0.008,label="y=0.0167 - 0.0002249x \n R-squared: 0.207",size=5)+
  xlab("Latitude (°N)")+
  theme_minimal() 
fig_vs_and_latitude



fig_vs_and_latitude_flip<-fig_vs_and_latitude+coord_flip() # in order to put together in illustrator
ggsave("figs/LOCAL_ADAPTATION_Vs_latitude_narrow.pdf",plot = fig_vs_and_latitude_flip,height = 6,width = 5,device = "pdf")
ggsave("figs/LOCAL_ADAPTATION_Vs_latitude.pdf",plot=fig_vs_and_latitude, height = 5,width = 6,device = "pdf")

# Make the plot against longitude
fig_vs_and_longitude<-
  ggplot(data=ecotype_info)+
  # geom_point(aes(x=longitude, y=1/Vs),shape=16,color="black")+
  geom_point(aes(x=longitude_ecotype, y=inv_Vs, color=inv_Vs),shape=16,cex=3)+
  # scale_color_gradientn(colors=rev(brewer.pal(9,"Spectral")))+
  scale_color_gradientn(colors=c(rev(RColorBrewer::brewer.pal(5,"Reds")),RColorBrewer::brewer.pal(5,"Greens")),
                        values = rescale(c(0,0.001,0.002,0.004,0.006,0.008,0.011,0.015,0.02)))+
  stat_smooth(aes(x=longitude_ecotype, y=inv_Vs, color=inv_Vs),method = "glm", color='grey')+
  # xlim(30,65)+
  xlab("Longitude (°E)")+
  theme_minimal() 
fig_vs_and_longitude
# fig_vs_and_longitude<-fig_vs_and_longitude+coord_flip() # in order to put together in illustrator
ggsave("figs/LOCAL_ADAPTATION_Vs_longitude.pdf",plot=fig_vs_and_longitude, height = 5,width = 6,device = "pdf")


################################################################################
# Compare to bio1 and bio12

fig_vs_bio1<-
  ggplot(ecotype_info,
         aes(y=inv_Vs, x=bio1_ecotype,color=inv_Vs))+
  geom_point(cex=3)+
  scale_color_gradientn(colors=c(rev(RColorBrewer::brewer.pal(5,"Reds")),RColorBrewer::brewer.pal(5,"Greens")),
                        values = rescale(c(0,0.001,0.002,0.004,0.006,0.008,0.011,0.015,0.02)))+
  stat_smooth(formula = y~x, method="glm", col="grey")+
  
  xlab("Temperature (°C)")+
  annotate(geom = "text",x=0,y=0.015,label="y=0.00164+0.000489x \n R-squared: 0.3993",cex=5)+
  theme_minimal(base_size = 18)
fig_vs_bio1 

ggsave(filename = "figs/LOCAL_ADAPTATION_inv_Vs_ecotype_bio1.pdf",plot = fig_vs_bio1, height = 5,width = 6,device = "pdf")



################################################################################

# Compare Vs with SDM
suitability<-raster::extract(sdm,y = cbind(ecotype_info$longitude_ecotype, ecotype_info$latitude_ecotype))
ecotype_info$suitability=suitability

fig_vs_suitability<-
  ggplot(ecotype_info,
         aes(y=inv_Vs, x=suitability,color=inv_Vs))+
  geom_point(cex=3)+
  scale_color_gradientn(colors=c(rev(RColorBrewer::brewer.pal(5,"Reds")),RColorBrewer::brewer.pal(5,"Greens")),
                        values = rescale(c(0,0.001,0.002,0.004,0.006,0.008,0.011,0.015,0.02)))+
  #stat_smooth(formula = y~x, method="lm", col="grey")+
  #stat_smooth(method="glm", col="grey",formula=y~poly(x,2))+
  xlab("Suitability")+
  geom_rug(aes(x=suitability),col=  "black" )+
  theme_minimal(base_size = 18)
fig_vs_suitability 
ggsave(filename = "figs/LOCAL_ADAPTATION_inv_Vs_ecotype_suitability.pdf",plot = fig_vs_suitability, height = 5,width = 6,device = "pdf")

