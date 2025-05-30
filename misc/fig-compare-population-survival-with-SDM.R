################################################################################
### Goal
### Compare species distribution model with population survival
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
## Survival data
sur=read.csv("data/survival.csv")

## Survival predicted
popsize=read.csv("pop_size_estimation/pop_size_estimations.csv")
popsizeraw=popsize

## Cliamte 
load("grene/data/worldclim_sitesdata.rda")
load("grene/data/locations_data.rda")
  # locations_data %>% View
################################################################################
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

###### 
## Merge both survival proportion and estimated popualtion size
# popsizesur<-merge(sur,popsizeraw)


################################################################################
# PLOT map against survivals

# Load the SDM for 
load("data-external/sdm.rda")
sdmcolor=brewer.pal(9,"Greys")[-1]


# Generate the colors
colors <- colorRampPalette(brewer.pal(9, "Spectral"))(length(unique(toplot$Vs)))

# Plot the points with the colors
plot(sdm, col=sdmcolor)
# points(x = sur$longitude, y=sur$latitude , 
#        col=rev(colors)[rank(1/sursummary$survival_proportion)] , 
#        pch=16)


# Get the SDM value at the geographic locations
sdmvals=extract(sdm,cbind(sursummary$longitude,sursummary$latitude))
sursummary$SDM<-sdmvals

fig_smd_survival_replicates<-
ggplot(sursummary)+
  geom_point(aes(x=SDM, y=survival_proportion))+
  stat_smooth(aes(x=SDM, y=survival_proportion), method = 'glm', color="grey")+
  xlab("SDM habitat suitability")+
  ylab("GrENE-net proportion of replicate survival")+
  theme_minimal()
fig_smd_survival_replicates

save_plot("figs/fig_smd_survival_replicates.pdf",fig_smd_survival_replicates, base_width = 7,base_height =5)
save_plot("figs/fig_smd_survival_replicates.png",fig_smd_survival_replicates, base_width = 7,base_height =5)


# Get the SDM value at the geographic locations
sdmvals=extract(sdm,cbind(popsizesummary$longitude,popsizesummary$latitude))
popsizesummary$SDM<-sdmvals

ggplot(popsizesummary)+
  geom_point(aes(x=SDM, y=avg_flowers))+
  stat_smooth(aes(x=SDM, y=avg_flowers), method = 'glm', color="grey")+
  # scale_x_log10()+
  # scale_x_log10()+
  xlab("SDM habitat suitability")+
  ylab("GrENE-net sum num. flowers")+
  theme_minimal()->
  fig_smd_sum_flowers
fig_smd_sum_flowers
save_plot("figs/fig_smd_and_sum_flowers_collected.pdf",fig_smd_sum_flowers, base_width = 7,base_height =5)
save_plot("figs/fig_smd_and_sum_flowers_collected.png",fig_smd_sum_flowers, base_width = 7,base_height =5)


ggplot(popsizesummary)+
  geom_point(aes(x=SDM, y=avg_popsize))+
  stat_smooth(aes(x=SDM, y=avg_popsize), method = 'glm', color="grey")+
  xlab("SDM habitat suitability")+
  ylab("GrENE-net sum num. plants")+
  theme_minimal()->
  fig_smd_sum_plants
fig_smd_sum_plants
save_plot("figs/fig_smd_and_sum_plants_estimated.pdf",fig_smd_sum_plants, base_width = 7,base_height =5)
save_plot("figs/fig_smd_and_sum_plants_estimated.png",fig_smd_sum_plants, base_width = 7,base_height =5)
