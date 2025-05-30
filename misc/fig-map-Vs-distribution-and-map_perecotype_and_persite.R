library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)

myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives/MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
## Raw data

load("grene/data/worldclim_ecotypesdata.rda")
load("grene/data/worldclim_sitesdata.rda")
load("grene/data/ecotypes_data.rda")

## Intermediate data
Vsresults<-read.table("data-intermediate/stabilizing_selection_wmax_vs.txt",header = T) %>% 
          dplyr::rename(Vs=vs_mean, wmax=wmax_mean)

# Vssite<-read.csv("data-intermediate/stabilizing_selection_simple_vspersite.csv"
#                    ,header = T)
load("data-intermediate/vssitesterminal.rda")
load("data-intermediate/vssitesyears.rda")
Vssite=vssitesyears
Vssite=vssitesterminal

# Generate joint matrix of Vs per ecotype
tmp= merge(ecotypes_data, by.x="ecotypeid", worldclim_ecotypesdata, by.y="ecotypeid")
dfeco= merge(Vsresults, by.x="ecotype", tmp, by.y="ecotypeid")


################################################################################
# Simple GGplot to start checking
customredgreen=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) )
customredgreen=rev(customredgreen)

Vssite<-
  Vssite %>% dplyr::mutate(Vssitecorrected=ifelse(Vssite<0, Vssite,NA))

library(ggnewscale)
ggplot()+
  geom_point(data=dfeco, aes(x=longitude, y=latitude, color=-1/Vs), 
             shape=3)+
  scale_color_gradientn(colors= customredgreen)+
  new_scale_color()+
  geom_point(data=Vssite, aes(x=longitude, y=latitude, color=Vssitecorrected))+
  # scale_color_gradient2(midpoint = 0,high = "black")+
  scale_color_gradientn(colors=rev(brewer.pal(9,"Greys")))+
  theme_minimal()


ggplot()+
  geom_point(data=dfeco, aes(x=longitude, y=latitude, color=-1/Vs), 
             shape=3)+
  scale_color_gradientn(colors= customredgreen)+
  geom_point(data=Vssite, aes(x=longitude, y=latitude, color=Vssitecorrected),
             size=3
             )+
  # new_scale_color()+
  # scale_color_gradient2(midpoint = 0,high = "black")+
  # scale_color_gradientn(colors=rev(brewer.pal(9,"Greys")))+
  theme_minimal()

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

# # Generate the colors
# # Plot the points with the colors
# colors <- colorRampPalette(brewer.pal(9, "Spectral"))(length(unique(toplot$Vs)))
# 
# plot(sdm, col=sdmcolor)
# points(x = toplot$longitude, y=toplot$latitude , 
#        col=rev(colors)[rank(1/toplot$Vs)] , 
#        pch=16)
# pdf("figs/fig-map-Vs-ecotype-mcmcglmm.pdf",height = 8,width = 9)
# plot(sdm, col=sdmcolor)
# points(x = toplot$longitude, y=toplot$latitude , 
#        col=rev(colors)[rank(1/toplot$Vs)] , 
#        pch=16)
# dev.off()

# Plot the SDM with a high quality function
source("functions-utilities-map-highquality.R")
basemap<- coolmaptile(sdm)


################################################################################
# PLOT VS on a map

# Make the geographic plot
map_and_vs<-
basemap+
  geom_point(data=dfeco, aes(x=longitude, y=latitude, color=1/Vs))+
  # scale_color_gradientn(colors=rev(brewer.pal(9,"Spectral")))
  scale_color_gradientn(colors=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ))
  
map_and_vs
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-map.pdf",map_and_vs, base_height = 4,base_width = 6)
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-map.png",map_and_vs, base_height = 5,base_width = 7)


# Make the plot against latitude
fig_vs_and_latitude<-
ggplot(data=dfeco)+
  # geom_point(aes(x=latitude, y=1/Vs),shape=16,color="black")+
  geom_point(aes(x=latitude, y=1/Vs, color=1/Vs),shape=16)+
  # scale_color_gradientn(colors=rev(brewer.pal(9,"Spectral")))+
  scale_color_gradientn(colors=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ))+
  stat_smooth(aes(x=latitude, y=1/Vs, color=1/Vs),method = "glm", color='grey')+
  stat_smooth(aes(x=latitude, y=1/Vs, color=1/Vs),method = "glm", formula = y~poly(x,2), color='grey')+
  # xlim(30,65)+
  xlab("Latitude (°N)")+
  theme_minimal() 
fig_vs_and_latitude
fig_vs_and_latitude_flip<-fig_vs_and_latitude+coord_flip() # in order to put together in illustrator
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-latitude-narrow.pdf",fig_vs_and_latitude_flip, base_height = 4,base_width = 3)
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-latitude-narrow.png",fig_vs_and_latitude_flip, base_height = 4,base_width = 3)
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-latitude.pdf",fig_vs_and_latitude, base_height = 4,base_width = 5)
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-latitude.png",fig_vs_and_latitude, base_height = 4,base_width = 5)

# Make the plot against longitude
fig_vs_and_longitude<-
  ggplot(data=dfeco)+
  # geom_point(aes(x=longitude, y=1/Vs),shape=16,color="black")+
  geom_point(aes(x=longitude, y=1/Vs, color=1/Vs),shape=16)+
  # scale_color_gradientn(colors=rev(brewer.pal(9,"Spectral")))+
  scale_color_gradientn(colors=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ))+
  stat_smooth(aes(x=longitude, y=1/Vs, color=1/Vs),method = "glm", color='grey')+
  stat_smooth(aes(x=longitude, y=1/Vs, color=1/Vs),method = "glm", formula = y~poly(x,2), color='grey')+
  # xlim(30,65)+
  xlab("Longitude (°E)")+
  theme_minimal() 
fig_vs_and_longitude
# fig_vs_and_longitude<-fig_vs_and_longitude+coord_flip() # in order to put together in illustrator
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-longitude.pdf",fig_vs_and_longitude, base_height = 4,base_width = 5)
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-longitude.png",fig_vs_and_longitude, base_height = 4,base_width = 5)


################################################################################
# Compare to bio1 and bio12

fig_vs_bio1<-
  ggplot(dfeco,
         aes(y=1/Vs, x=bio1, color=1/Vs))+
  geom_point()+
  scale_color_gradientn("1/Vs",colors=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ))+
  stat_smooth(formula = y~poly(x,2), method="glm", col="grey")+
  stat_smooth(formula = y~x, method="glm", col="grey")+
  xlab("Temperature (°C)")+
  theme_minimal()
fig_vs_bio1
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-bio1.pdf",fig_vs_bio1, base_height = 5,base_width = 6)
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-bio1.png",fig_vs_bio1, base_height = 5,base_width = 6)


fig_vs_bio1_wmaxsize<-
  ggplot(dfeco,
         aes(y=1/Vs, x=bio1, color=1/Vs, size=wmax))+
  geom_point()+
  scale_color_gradientn("1/Vs",colors=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ))+
  stat_smooth(formula = y~poly(x,2), method="glm", col="grey")+
  stat_smooth(formula = y~x, method="glm", col="grey")+
  xlab("Temperature (°C)")+
  theme_minimal()
fig_vs_bio1_wmaxsize
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-bio1-wmaxsize.pdf",fig_vs_bio1_wmaxsize, base_height = 5,base_width = 6)
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-bio1-wmaxsize.png",fig_vs_bio1_wmaxsize, base_height = 5,base_width = 6)




fig_vs_bio12<-
  ggplot(dfeco,
         aes(y=1/Vs, x=bio12, color=1/Vs))+
  geom_point()+
  scale_color_gradientn("1/Vs",colors=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ))+
  stat_smooth(formula = y~poly(x,2), method="glm", col="grey")+
  stat_smooth(formula = y~x, method="glm", col="grey")+
  xlab("Precipitation (mm)")+
  theme_minimal()
fig_vs_bio12
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-bio12.pdf",fig_vs_bio12, base_height = 5,base_width = 6)
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-bio12.png",fig_vs_bio12, base_height = 5,base_width = 6)


fig_vs_bio18<-
  ggplot(dfeco,
         aes(y=1/Vs, x=bio18, color=1/Vs))+
  geom_point()+
  scale_color_gradientn("1/Vs",colors=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ))+
  stat_smooth(formula = y~poly(x,2), method="glm", col="grey")+
  stat_smooth(formula = y~x, method="glm", col="grey")+
  xlab("Precipitation (mm)")+
  theme_minimal()
fig_vs_bio18
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-bio18.pdf",fig_vs_bio18, base_height = 5,base_width = 6)
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-bio18.png",fig_vs_bio18, base_height = 5,base_width = 6)


################################################################################

# Compare Vs with SDM
suitability<-raster::extract(sdm,y = cbind(dfeco$longitude, dfeco$latitude))
dfeco$suitability=suitability

fig_vs_suitability<-
ggplot(dfeco)+
  geom_point(aes(y=1/Vs, x=suitability, color=1/Vs))+
  scale_color_gradientn("1/Vs",colors=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ))+
  stat_smooth(aes(y=1/Vs, x=suitability), formula = y~poly(x,2), method="glm", col="grey")+
  stat_smooth(aes(y=1/Vs, x=suitability), formula = y~x, method="glm", col="grey")+
  geom_rug(aes(x=suitability),col=  "black" )+
  # geom_rug(aes(x=bio1),col= color_palette() )+
  theme_minimal()
fig_vs_suitability
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-sdm-habitat-suitability.pdf",fig_vs_suitability, base_height = 5,base_width = 6)
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-sdm-habitat-suitability.png",fig_vs_suitability, base_height = 5,base_width = 6)


fig_vs_suitability_wmaxsize<-
  ggplot(dfeco)+
  geom_point(aes(y=1/Vs, x=suitability, color=1/Vs, size=wmax))+
  scale_color_gradientn("1/Vs",colors=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ))+
  stat_smooth(aes(y=1/Vs, x=suitability), formula = y~poly(x,2), method="glm", col="grey")+
  stat_smooth(aes(y=1/Vs, x=suitability), formula = y~x, method="glm", col="grey")+
  geom_rug(aes(x=suitability),col=  "black" )+
  # geom_rug(aes(x=bio1),col= color_palette() )+
  theme_minimal()
fig_vs_suitability
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-sdm-habitat-suitability-wmaxsize.pdf",fig_vs_suitability_wmaxsize, base_height = 5,base_width = 6)
save_plot("figs/fig-Vs-ecotype-mcmcglmm-and-sdm-habitat-suitability-wmaxsize.png",fig_vs_suitability_wmaxsize, base_height = 5,base_width = 6)

################################################################################
################################################################################
## Now Wmax and relationship with Vs
################################################################################


fig_wmax_latitude<-
ggplot(data=dfeco, aes(y=wmax, x=latitude, color=wmax))+
  geom_point()+
  scale_color_gradientn(colors=(brewer.pal(9,"BrBG")))+
  stat_smooth(method="glm", col="grey", formula=y~poly(x,2))+
  stat_smooth(method="glm", col="grey")+
  ylab("Wmax")+
  xlab("Latitude (°N)")+
  theme_minimal()
fig_wmax_latitude


fig_vs_latitude<-
ggplot(data=dfeco, aes(y=1/Vs, x=latitude, color=1/Vs))+
  geom_point()+
  scale_color_gradientn(colors=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ))+
  stat_smooth(method="glm", col="grey", formula=y~poly(x,2))+
  stat_smooth(method="glm", col="grey")+
  xlab("Latitude (°N)")+
  theme_minimal()
fig_vs_latitude


fig_wmax_longitude<-
  ggplot(data=dfeco, aes(y=wmax, x=longitude, color=wmax))+
  geom_point()+
  scale_color_gradientn(colors=(brewer.pal(9,"BrBG")))+
  stat_smooth(method="glm", col="grey", formula=y~poly(x,2))+
  stat_smooth(method="glm", col="grey")+
  ylab("Wmax")+
  xlab("Longitude (°E)")+
  theme_minimal()
fig_wmax_longitude


fig_vs_longitude<-
  ggplot(data=dfeco, aes(y=1/Vs, x=longitude, color=1/Vs))+
  geom_point()+
  scale_color_gradientn(colors=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ))+
  stat_smooth(method="glm", col="grey", formula=y~poly(x,2))+
  stat_smooth(method="glm", col="grey")+
  xlab("Longitude (°E)")+
  theme_minimal()
fig_vs_longitude

fig_wmax_vs<-
ggplot(data=dfeco, aes(y=1/Vs, x=wmax, color=1/Vs))+
  geom_point()+
  scale_color_gradientn(colors=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ))+
  xlab("Wmax")+
  ylab("1/Vs")+
  stat_smooth(method="glm", col="grey")+
  theme_minimal()
fig_wmax_vs

fig_vs_wmax_distribution<- fig_wmax_
fig_wmax_latitude<- + fig_vs_latitude + fig_wmax_vs
save_plot("figs/fig-Vs-wmax-ecotype-mcmcglmm-distribution.pdf",fig_vs_wmax_distribution, base_height = 4,base_width = 10)
save_plot("figs/fig-Vs-wmax-ecotype-mcmcglmm-distribution.png",fig_vs_wmax_distribution, base_height = 4,base_width = 10)

save_plot("figs/fig-Vs-wmax-ecotype-mcmcglmm-distribution-a.pdf",fig_wmax_latitude, base_height = 4,base_width = 5)
save_plot("figs/fig-Vs-wmax-ecotype-mcmcglmm-distribution-b.pdf",fig_vs_latitude, base_height = 4,base_width = 5)
save_plot("figs/fig-Vs-wmax-ecotype-mcmcglmm-distribution-c.pdf",fig_wmax_vs, base_height = 4,base_width = 5)
save_plot("figs/fig-Vs-wmax-ecotype-mcmcglmm-distribution-a.png",fig_wmax_latitude, base_height = 4,base_width = 5)
save_plot("figs/fig-Vs-wmax-ecotype-mcmcglmm-distribution-b.png",fig_vs_latitude, base_height = 4,base_width = 5)
save_plot("figs/fig-Vs-wmax-ecotype-mcmcglmm-distribution-c.png",fig_wmax_vs, base_height = 4,base_width = 5)


save_plot("figs/fig-Vs-wmax-ecotype-mcmcglmm-distribution-e.pdf",fig_vs_longitude, base_height = 4,base_width = 5)
save_plot("figs/fig-Vs-wmax-ecotype-mcmcglmm-distribution-f.pdf",fig_wmax_longitude, base_height = 4,base_width = 5)
save_plot("figs/fig-Vs-wmax-ecotype-mcmcglmm-distribution-e.png",fig_vs_longitude, base_height = 4,base_width = 5)
save_plot("figs/fig-Vs-wmax-ecotype-mcmcglmm-distribution-f.png",fig_wmax_longitude, base_height = 4,base_width = 5)
