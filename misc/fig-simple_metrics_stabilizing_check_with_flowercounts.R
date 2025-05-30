################################################################################
### Goal
### Estimate some basic metrics of ecotypes increasing in frequency across
### sites, means, variances, optimal places. Including population size

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
# Load the long format with ecotype frequencies and climates

load(file=paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize.rda"))
ls()

e<-enew <- ecotype_allyears_frequencies_long_raw_climate_popsize
enew[1:5,1:10]


# Get some general info, like what is the per ecotype average
simplemetrics <- 
  e %>% 
  dplyr::mutate(count=freq*flowerscollected_corrected) %>% 
  dplyr::group_by(id,year) %>% 
  summarise(avg_freq = mean(count),
            var_freq = var(count)/mean(count),
            opt_clim = weighted.mean(sites_bio1,count/sum(count))
  )
# add back all the other variables
simplemetrics<-
  merge(simplemetrics,e,by=c("id","year")) #

# Plot 
ggplot(simplemetrics ) +
  geom_point(aes(x=avg_freq,y=var_freq), shape=16)+
  stat_smooth(aes(x=avg_freq,y=var_freq), color = "darkgrey")+
  stat_smooth(aes(x=avg_freq,y=var_freq),method = "glm",color = "darkgrey")+
  xlab("Mean in frequency GrENE-net (mean)")+  ylab("CV in frequency GrENE-net")+
  theme_minimal()

ggplot(simplemetrics ) +
  geom_point(aes(x=avg_freq,y=var_freq), shape=16)+
  stat_smooth(aes(x=avg_freq,y=var_freq), color = "darkgrey")+
  stat_smooth(aes(x=avg_freq,y=var_freq),method = "glm",color = "darkgrey")+
  scale_x_log10()+
  scale_y_log10()+
  xlab("Mean in frequency GrENE-net (mean)")+
  ylab("CV in frequency GrENE-net")+
  theme_minimal()

ggplot(simplemetrics ,
       aes(y=opt_clim,x=ecotypes_bio1,color=avg_freq,size=avg_freq)) +
  # geom_vline(xintercept= c(unique(e$sites_bio1)),lty='dotted' )+
  # geom_hline(yintercept= c(unique(e$sites_bio1)),lty='dotted' )+
  geom_point(color='black')+
  geom_point()+
  scale_color_gradientn("Mean frequency \n across sites",colors=brewer.pal(9,"Greens"))+
  stat_smooth(method='glm')+
  stat_smooth(method='glm', formula = y~poly(x,2), col='darkgrey')+
  geom_abline(intercept = 0,slope = 1, col='darkgrey')+
  xlim(c(-5,+24))+
  ylim(c(-5,+24))+
  ylab("Temperature optimum from GrENE-net (C)")+
  xlab("Temperature of origin (C)")+
  guides(size = "none") +
  theme_minimal()
