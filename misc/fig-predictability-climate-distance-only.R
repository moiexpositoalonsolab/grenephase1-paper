################################################################################
### Goal
### Rank sites based on spearman r prediction of cliamte distance

################################################################################
library(tidyverse)
library(dplyr)
library(cowplot)
library(RColorBrewer)

################################################################################

# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives/MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses" # GOOGLE DRVIE

setwd(myfolder)

################################################################################
#### Read the climate and ecotype frequency

# load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate.csv")
load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate.rda")

e=ecotype_terminal_frequencies_long_raw_climate_popsize


#### 

# Correlations using log10 p1 / p0
rhos<-
  e %>% 
  dplyr::filter(longitude> -15) %>% 
  # dplyr::filter(latitude > 32) %>% 
  dplyr::filter(max_year >=2) %>%
  group_by(site,rep,latitude,longitude,altitude) %>% 
  summarise(rho=- cor.test(log10(maxfreq/startfreq),(sites_bio1-ecotypes_bio1)^2 , method='s')$estimate)
ggplot(rhos)+
  geom_point(aes(y=rho,x=latitude, group=site))+
  stat_summary(aes(y=rho,x=latitude, group=site), col='grey')+
  stat_smooth(aes(y=rho,x=latitude),method='glm',color='grey')+
  stat_smooth(aes(y=rho,x=latitude),method='glm', formula = y~poly(x,2), color='grey')+
  ylab("Spearman's r (climate distance)")+
  xlab("Latitude (°N)")+
  geom_hline(yintercept=0,lty='dotted')+
  theme_minimal()

rhos<-
  e %>% 
  dplyr::filter(longitude> -15) %>% 
  dplyr::filter(latitude > 32) %>% 
  group_by(site,rep,latitude,longitude,altitude) %>% 
    summarise(rho=- cor.test(log10(maxfreq/startfreq),(sites_bio1-ecotypes_bio1)^2 , method='s')$estimate)
ggplot(rhos)+
  geom_point(aes(y=rho,x=latitude, group=site))+
  stat_summary(aes(y=rho,x=latitude, group=site), col='grey')+
  stat_smooth(aes(y=rho,x=latitude),method='glm',color='grey')+
  stat_smooth(aes(y=rho,x=latitude),method='glm', formula = y~poly(x,2), color='grey')+
  ylab("Spearman's r (climate distance)")+
  xlab("Latitude (°N)")+
  geom_hline(yintercept=0,lty='dotted')+
  theme_minimal()


ggplot(rhos)+
  geom_point(aes(y=rho,x=altitude, group=site))+
  stat_summary(aes(y=rho,x=altitude, group=site), col='grey')+
  stat_smooth(aes(y=rho,x=altitude),method='glm',color='grey')+
  stat_smooth(aes(y=rho,x=altitude),method='glm', formula = y~poly(x,2), color='grey')+
  ylab("Spearman's r (climate distance)")+
  xlab("Latitude (°N)")+
  geom_hline(yintercept=0,lty='dotted')+
  theme_minimal()

ggplot(rhos)+
  geom_point(aes(y=rho,x=longitude, group=site))+
  stat_summary(aes(y=rho,x=longitude, group=site), col='grey')+
  stat_smooth(aes(y=rho,x=longitude),method='glm',color='grey')+
  stat_smooth(aes(y=rho,x=longitude),method='glm', formula = y~poly(x,2), color='grey')+
  ylab("Spearman's r (climate distance)")+
  xlab("Latitude (°N)")+
  geom_hline(yintercept=0,lty='dotted')+
  theme_minimal()


# Correlations using p1/p0
rhos<-
  e %>% 
  dplyr::filter(longitude> -15) %>% 
  group_by(site,rep,latitude,longitude,altitude) %>% 
  summarise(rho=- cor.test((maxfreq/startfreq),(sites_bio1-ecotypes_bio1)^2 , method='s')$estimate)
ggplot(rhos)+
  geom_point(aes(y=rho,x=latitude, group=site))+
  stat_summary(aes(y=rho,x=latitude, group=site), col='grey')+
  stat_smooth(aes(y=rho,x=latitude),method='glm',color='grey')+
  stat_smooth(aes(y=rho,x=latitude),method='glm', formula = y~poly(x,2), color='grey')+
  ylab("Spearman's r (climate distance)")+
  xlab("Latitude (°N)")+
  geom_hline(yintercept=0,lty='dotted')+
  theme_minimal()

# Correlation using delta p1-p0
rhos2<-
  e %>% 
  dplyr::filter(longitude> -15) %>% 
  group_by(site,rep,latitude,longitude,altitude) %>% 
  summarise(rho=- cor.test(maxfreq-startfreq, (sites_bio1-ecotypes_bio1)^2 , method='s')$estimate)
ggplot(rhos2)+
  geom_point(aes(y=rho,x=latitude, group=site))+
  stat_summary(aes(y=rho,x=latitude, group=site), col='grey')+
  stat_smooth(aes(y=rho,x=latitude),method='glm',color='grey')+
  stat_smooth(aes(y=rho,x=latitude),method='glm', formula = y~poly(x,2), color='grey')+
  ylab("Spearman's r (climate distance)")+
  xlab("Latitude (°N)")+
  geom_hline(yintercept=0,lty='dotted')+
  theme_minimal()


################################################################################
#### Read the climate and ecotype frequency for all years 
# Decided it is better to use all years for this visualization
# load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate.csv")
# load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate.rda")

e <- ecotype_allyears_frequencies_long_raw_climate_popsize

#### 

# Correlations using log10 p1 / p0
rhos<-
  e %>% 
  dplyr::filter(longitude> -15) %>% 
  # dplyr::filter(latitude > 32) %>% 
  # dplyr::filter(maxyear >=2) %>%
  group_by(site,rep, year,latitude,longitude, sites_bio18) %>% 
  summarise(rho=- cor.test(log10(freq/startfreq),(sites_bio1-ecotypes_bio1)^2 , method='s')$estimate)
ggplot(rhos)+
  geom_point(aes(y=rho,x=latitude, group=site, color=as.numeric(year)))+
  stat_summary(aes(y=rho,x=latitude, group=site), col='grey')+
  # stat_smooth(aes(y=rho,x=latitude, group=year,color=year ),method='glm')+
  stat_smooth(aes(y=rho,x=latitude, group=year,color=as.numeric(year)),method='glm', formula = y~poly(x,2))+
  scale_color_gradientn("year", colors=brewer.pal(4,"Greens")[-1])+
  ylab("Spearman's r (climate distance)")+
  xlab("Latitude (°N)")+
  geom_hline(yintercept=0,lty='dotted')+
  theme_minimal()->
  fig_predictability_ecotype_delta_by_climate_naive_bio1
fig_predictability_ecotype_delta_by_climate_naive_bio1

save_plot("figs/fig_predictability_ecotype_delta_by_climate_naive_bio1_across3years.pdf",fig_predictability_ecotype_delta_by_climate_naive_bio1, base_width = 9,base_height =7)
save_plot("figs/fig_predictability_ecotype_delta_by_climate_naive_bio1_across3years.png",fig_predictability_ecotype_delta_by_climate_naive_bio1, base_width = 9,base_height =7)

flowerscollected<-
ggplot(e)+
  geom_point(aes(y=flowerscollected_corrected,x=latitude, group=site, color=as.numeric(year)))+
  stat_summary(aes(y=flowerscollected_corrected,x=latitude, group=site), col='grey')+
  # stat_smooth(aes(y=flowerscollected_corrected,x=latitude, group=year,color=as.numeric(year)),method='glm', formula = y~poly(x,2))+
  scale_color_gradientn("year", colors=brewer.pal(4,"Greens")[-1])+
  ylab("Number of flowers collected")+
  xlab("Latitude (°N)")+
  geom_hline(yintercept=0,lty='dotted')+
  theme_minimal()

save_plot("figs/fig_predictability_ecotype_delta_by_climate_naive_bio1_across3years-flowersamplesize.pdf",
          flowerscollected, base_width = 9,base_height =7)
save_plot("figs/fig_predictability_ecotype_delta_by_climate_naive_bio1_across3years-flowersamplesize.png",
          flowerscollected, base_width = 9,base_height =7)


################################################################################
#### Compare with population size


rhospop<-
  e %>% 
  dplyr::filter(longitude> -15) %>% 
  # dplyr::filter(latitude > 32) %>% 
  # dplyr::filter(maxyear >=2) %>%
  group_by(site,rep, year,latitude,longitude, sites_bio18,flowerscollected_corrected) %>% 
  summarise(rho=- cor.test(log10(freq/startfreq),(sites_bio1-ecotypes_bio1)^2 , method='s')$estimate)

ggplot(rhospop)+
  geom_point(aes(y=rho,x=latitude, group=site, color=as.numeric(year), size=flowerscollected_corrected))+
  stat_summary(aes(y=rho,x=latitude, group=site), col='grey')+
  # stat_smooth(aes(y=rho,x=latitude, group=year,color=year ),method='glm')+
  stat_smooth(aes(y=rho,x=latitude, group=year,color=as.numeric(year)),method='glm', formula = y~poly(x,2))+
  scale_color_gradientn("year", colors=brewer.pal(4,"Greens")[-1])+
  ylab("Spearman's r (climate distance)")+
  xlab("Latitude (°N)")+
  geom_hline(yintercept=0,lty='dotted')+
  theme_minimal()

ggplot(rhospop)+
  geom_point(aes(size=rho,x=latitude, group=site, color=as.numeric(year), y=flowerscollected_corrected))+
  stat_summary(aes(y=flowerscollected_corrected,x=latitude, group=site), col='grey')+
  # stat_smooth(aes(y=rho,x=latitude, group=year,color=year ),method='glm')+
  stat_smooth(aes(y=flowerscollected_corrected,x=latitude, group=year,color=as.numeric(year)),method='glm', formula = y~poly(x,2))+
  scale_color_gradientn("year", colors=brewer.pal(4,"Greens")[-1])+
  ylab("Spearman's r (climate distance)")+
  xlab("Latitude (°N)")+
  geom_hline(yintercept=0,lty='dotted')+
  theme_minimal()->
fig_flowerscollected_predictability_ecotype_delta_by_climate_naive_bio1

save_plot("figs/fig_flowerscollected_predictability_ecotype_delta_by_climate_naive_bio1_across3years.pdf",fig_flowerscollected_predictability_ecotype_delta_by_climate_naive_bio1, base_width = 9,base_height =7)
save_plot("figs/fig_flowerscollected_predictability_ecotype_delta_by_climate_naive_bio1_across3years.png",fig_flowerscollected_predictability_ecotype_delta_by_climate_naive_bio1, base_width = 9,base_height =7)

