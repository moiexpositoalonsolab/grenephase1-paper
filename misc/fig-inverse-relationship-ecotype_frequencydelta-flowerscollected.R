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
library(latex2exp)
################################################################################
#### Location
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
# myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
# Load the long format with ecotype frequencies and climates

load(file="data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize_flowers.rda")
ls()

e<-enew <- ecotype_allyears_frequencies_long_raw_climate_popsize_flowers
enew[1:5,1:10]
colnames(e) %>% tail

## How does the change in frequency
# Do the plot
frequency_and_flowers<-
ggplot(e)+
  geom_point(aes(y=(freq-startfreq)^2 , x=flowerscollected_corrected, color=as.numeric(year)), alpha=0.5)+
  # scale_color_continuous("year",type = 'viridis')+
  scale_color_gradientn("year",colors = brewer.pal(8,"Greens")[-1])+
  ylab("Frequency change (a1-a0)^2")+
  xlab("Number of flowers in pool")+
  theme_minimal()

frequency_and_flowers_log<-
e %>% 
  sample_n(1000) %>% 
ggplot(.)+
  geom_point(aes(y=(freq-startfreq)^2 , x=flowerscollected_corrected, color=as.numeric(year)), alpha=0.5)+
  scale_y_log10()+
  scale_x_log10()+
  # stat_smooth(aes(y=freq-startfreq , x=flowerscollected_corrected, color=as.numeric(year)), method='glm', color="grey")+
  # scale_color_continuous("year",type = 'viridis')+
  scale_color_gradientn("year",colors = brewer.pal(8,"Greens")[-1])+
  ylab("Frequency change (a1-a0)^2")+
  xlab("Number of flowers in pool")+
  theme_minimal()

freq_plot_grid<-
  plot_grid(frequency_and_flowers,frequency_and_flowers_log,ncol=2)
freq_plot_grid

save_plot("figs/fig-allele-frequency-inverse-relationship-flower-count.pdf",freq_plot_grid ,base_height = 6,base_width = 12)
save_plot("figs/fig-allele-frequency-inverse-relationship-flower-count.png",freq_plot_grid ,base_height = 6,base_width = 12)


# ################################################################################
# ## Variance in frequency
# evariance<-
#   e %>% 
#   dplyr::filter(year==1) %>% 
#   group_by(site, rep, longitude,latitude, sites_bio1, flowerscollected_corrected, totalplantnumber_complete) %>% 
#   dplyr::rename(p1=freq, p0=startfreq) %>% 
#   summarise(varp= var( (p1-p0)/(p0*(1-p0) ) ) ) 
# # %>% 
# #   left_join(
# #     .,
# #     dplyr::select(e,site,longitude,latitude, starts_with("sites_"),starts_with("flower"), starts_with("plant")),
# #     by='site'
# #   )
#               
# 
# ggplot(evariance,
#        aes(y= varp , x= latitude, alpha=0.5) )+
#   geom_point()+
#   stat_smooth(method='glm', formula = y~poly(x,2))+
#   # scale_color_continuous("year",type = 'viridis')+
#   ylab("Var (p1-p0)")+
#   xlab("Number of flowers in pool")+
#   theme_minimal()
# 
# ggplot(evariance,
#        aes(y= varp , x= flowerscollected_corrected, alpha=0.5) )+
#   geom_point()+
#   stat_smooth(method='glm', formula = y~poly(x,2))+
#   # scale_color_continuous("year",type = 'viridis')+
#   ylab("Var (p1-p0)")+
#   xlab("Number of flowers in pool")+
#   theme_minimal()

