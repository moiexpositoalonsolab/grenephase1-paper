################################################################################
### Goal
### Example file to generate a future projection of vulnerability
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

## Raw data
file_edelta = 'data/delta_ecotype_freq.txt' 
file_e0 = 'data/merged_ecotype_frequency.txt'
file_eid = 'data/greneNet_final_v1.0.fam'
file_climate_sites= "grene/data/worldclim_sitesdata.rda"
file_info_sites=  "grene/data/sites.clim.rda"
files_climate_ecotypes="grene/data/worldclim_ecotypesdata.rda"
files_info_ecotypes="grene/data/ecotypes_data.rda"

load(files_climate_ecotypes)
load(file_climate_sites)
load(files_info_ecotypes)

## Intermediate data
Vsresults<-read.table("data-intermediate/stabilizing_selection_wmax_vs.txt",header = T)


# Merge
tmp= merge(ecotypes_data, by.x="ecotypeid", worldclim_ecotypesdata, by.y="ecotypeid")
dfeco= merge(Vsresults, by.x="ecotype", tmp, by.y="ecotypeid") %>% 
        rename(Vs=vs_mean, wmax=wmax_mean)


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


################################################################################
# Prediction of vulnerability

# The most simple is to do the distance square
temperaturechange<- (fut[[1]]-now[[1]])^2
temperaturechangebad<- (futbad[[1]]-now[[1]])^2

# Extract the future trajectory for each ecotype
dfeco$temperaturechange<-extract(temperaturechange,cbind(dfeco$longitude,dfeco$latitude))
# Extract future trajectory bad scenario
dfeco$temperaturechangebad<-extract(temperaturechangebad,cbind(dfeco$longitude,dfeco$latitude))

### SAVE INTERMEDIATE



# Calculate a vulnerability metric
# dfeco$vulnerability <- dfeco$wmax * exp(-(dfeco$temperaturechange) / dfeco$Vs)
dfeco$vulnerability_perecotype <- exp(-(dfeco$temperaturechange) / dfeco$Vs) # W= Wmax * exp( - (z1-z0)^2 / Vs )
dfeco$vulnerability_general <- exp(-(dfeco$temperaturechange) / mean(dfeco$Vs)) # Using mean decay of fitness
# Add bad scenario
dfeco$vulnerability_perecotype_bad <- exp(-(dfeco$temperaturechangebad) / dfeco$Vs) # W= Wmax * exp( - (z1-z0)^2 / Vs )
dfeco$vulnerability_general_bad <- exp(-(dfeco$temperaturechangebad) / mean(dfeco$Vs)) # Using mean decay of fitness

dfeco<-na.omit(dfeco)

# Create a base
source("functions-utilities-map-highquality.R")
basemap<- coolmaptile(temperaturechange, rev(c('black',rev(brewer.pal(name='Greys',4))) ))
basemap

save_plot("figs/fig_map_change_bio1_2030_spp245_averageVs.pdf",basemap, base_width = 7,base_height =5)
save_plot("figs/fig_map_change_bio1_2030_spp245_averageVs.png",basemap, base_width = 7,base_height =5)

# Plot the map with the predicted

# Test
ggplot()+
  geom_point(
    data=dfeco,
    aes(
      y=latitude, x=longitude,
      color=vulnerability_general
    ))+
  # scale_color_gradientn(colours = rev(c("black",rev(brewer.pal(9,"Spectral")))))
  scale_color_gradientn(colours = c("black",rev(brewer.pal(5,"Reds")) ))
ggplot()+
  geom_point(
    data=dfeco,
    aes(
      y=latitude, x=longitude,
      color=vulnerability_general_bad
    ))+
  # scale_color_gradientn(colours = rev(c("black",rev(brewer.pal(9,"Spectral")))))
  scale_color_gradientn(colours = c("black",rev(brewer.pal(5,"Reds")) ))

# Vulnerability using average Vs across ecotypes
fig_vulnerability_averageecotype<-
basemap+
  geom_point(
    data=dfeco,
    aes(
      y=latitude, x=longitude,
      color=vulnerability_general
  ))+
  # scale_color_gradientn(colours = rev(c("black",rev(brewer.pal(9,"Spectral")))))
  scale_color_gradientn(colours = c("black",rev(brewer.pal(5,"Reds")) ))

fig_vulnerability_averageecotype<-
  fig_vulnerability_averageecotype+
  guides(color = guide_colorbar(
    # barwidth = 1,  # Adjust the width of the color bar
    # barheight = 10,  # Adjust the height of the color bar
    label.theme = element_text(size = 8)  # Adjust the size of the text labels
  )) +
  theme(
    # legend.key.size = unit(0.5, 'lines'),  # Adjust size of legend keys
    legend.text = element_text(size = 8)  # Adjust size of legend text
  )


save_plot("figs/fig_map_vulnerability_bio1_2030_spp245_averageVs.pdf",fig_vulnerability_averageecotype, base_width = 9,base_height =5)
save_plot("figs/fig_map_vulnerability_bio1_2030_spp245_averageVs.png",fig_vulnerability_averageecotype, base_width = 9,base_height =5)

# Vulnerability using ecotype specific Vs
fig_vulnerability_perecotype<-
  basemap+
  geom_point(
    data=dfeco,
    aes(
      y=latitude, x=longitude,
      color=vulnerability_perecotype
    ))+
  # scale_color_gradientn(colours = rev(c("black",rev(brewer.pal(9,"Spectral")))))
  scale_color_gradientn(colours = c("black",rev(brewer.pal(5,"Reds")) ))

fig_vulnerability_perecotype<-
fig_vulnerability_perecotype+
  guides(color = guide_colorbar(
    # barwidth = 1,  # Adjust the width of the color bar
    # barheight = 10,  # Adjust the height of the color bar
    label.theme = element_text(size = 8)  # Adjust the size of the text labels
  )) +
  theme(
    # legend.key.size = unit(0.5, 'lines'),  # Adjust size of legend keys
    legend.text = element_text(size = 8)  # Adjust size of legend text
  )

save_plot("figs/fig_map_vulnerability_bio1_2030_spp245_perecotypeVs.pdf",fig_vulnerability_perecotype, base_width = 9,base_height =5)
save_plot("figs/fig_map_vulnerability_bio1_2030_spp245_perecotypeVs.png",fig_vulnerability_perecotype, base_width = 9,base_height =5)


