################################################################################
### Goal
### Analyse flowering dates from grene-net

################################################################################
#### If transformed to collab can use the stuff below
# from google.colab import drive
# drive.mount('/content/drive')
# %load_ext rpy2.ipython
# %%R

################################################################################
#### Libraries
library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(patchwork)

################################################################################
#### Location if run locally
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
# myfolder<-"~/grenephase1-analyses"

myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
# Just simple sample data

load("grene/data/samples_data.rda")
load("grene/data/locations_data.rda")
load("grene/data/worldclim_sitesdata.rda")

samples_tmp<-samples_data %>%
  dplyr::filter(year<2021) %>%
  dplyr::mutate(year= year-2017) %>%
  dplyr::rename(flowers=flowerscollected) %>%
  dplyr::select(site,plot,year, date,month,day,flowers) %>%
  dplyr::mutate(julian_day = yday(ymd(date)))

flowering<-
  worldclim_sitesdata %>%
  merge(
    .,
    locations_data, by='site'
  ) %>%
  merge(
    .,
    samples_tmp, by="site"
  )

flowering<-
  flowering %>%
  mutate(date=ymd(date))

figallflowers<-
  ggplot(flowering)+
    geom_point(aes(x=latitude, y=date, color=flowers))+
    scale_color_gradientn("# flowers",colors=brewer.pal(9,"Greens")[-1])+
    ylab("Date")+
    xlab("Latitude (N)")+
    theme_minimal()
figallflowers

save_plot("figs/fig-date-flowering-and-latitude.pdf",figallflowers,base_height = 5, base_width = 7)
save_plot("figs/fig-date-flowering-and-latitude.png",figallflowers,base_height = 5, base_width = 7)


sink("tables/correlation-flowering-latitude.txt")
# print("warm neutral")

print("warm")
flowering %>%
  group_by(year) %>%
  summarise(cor=cor.test(julian_day,latitude,method = 'p')$estimate,
            p=cor.test(julian_day,latitude,method = 'p')$p.value,
            n=length(julian_day)
            )

sink()


################################################################################
# Load the ecotypes long format dataset with the most info about locations, climate, population size, flowers
load(file = "data-intermediate//ecotype_allyears_frequencies_long_raw_climate_popsize_flowers.rda")
load(file = "data-intermediate//ecotype_allyears_frequencies_long_raw_climate_popsize_flowers.rda")

e=ecotype_allyears_frequencies_long_raw_climate_popsize_flowers
################################################################################

# quantify average
estimated_flowering<-
  e %>%
  group_by(site, year, latitude, longitude, altitude) %>%
  # filter( julian_day< 244) %>% # separate things that flower before and after september
  # group_by(site, year, latitude, longitude, altitude, julian_day > 244) %>% # separate things that flower before and after september
  # dplyr::mutate(seasonday=ifelse(julian_day<264,julian_day+(365-264), julian_day-264 )) %>% # to exclude the Netherland?
  # summarise(avg_flowering=weighted.mean(w = flowers,seasonday))
summarise(avg_flowering=weighted.mean(w = flowers,julian_day))

# plot
flowering_correlations_average<-
ggplot(estimated_flowering,
       aes(y=avg_flowering, x=latitude, group=year, color=as.numeric(year))
       )+
  geom_point()+
  # stat_smooth(method='glm', col="grey")+
  stat_smooth(method='glm')+
  # scale_color_continuous("year",type = 'viridis')+
  scale_color_gradientn("year",colors=brewer.pal(9,"Greens")[-1])+
  # geom_hline(lty="dotted",yintercept = 1)+ # equinox
  # geom_hline(lty="dotted",yintercept = 355-264)+ # solstice
  # geom_hline(lty="dotted",yintercept = 80+(365-264))+ # march equinox
  # facet_wrap(~year)+
  labs(y="Flowering day")+
  labs(x="Latitude (N)")+
  theme_minimal()
flowering_correlations_average

# correlation test
source("functions-regressions.R")
sink("tables/test-correlation-latitude-grenenet-flowercollections,txt")
cat("Taking the average flowering [julian day] per grene-net site and every year and correlating with latitude")
estimated_flowering %>%
  ungroup() %>%
  group_by(year) %>%
  summarise(r=lm_eq(avg_flowering,latitude))
sink()


save_plot("figs/fig-date-flowering-average-and-latitude.pdf",flowering_correlations_average,base_height = 5, base_width = 7)
save_plot("figs/fig-date-flowering-average-and-latitude.png",flowering_correlations_average,base_height = 5, base_width = 7)


################################################################################
# Flowering and annual tempereature
estimated_flowering<-
  e %>%
  group_by(site, year, latitude, longitude, altitude, sites_bio1) %>%
  # filter( julian_day< 244) %>% # separate things that flower before and after september
  # group_by(site, year, latitude, longitude, altitude, julian_day > 244) %>% # separate things that flower before and after september
  # dplyr::mutate(seasonday=ifelse(julian_day<264,julian_day+(365-264), julian_day-264 )) %>% # to exclude the Netherland?
  # summarise(avg_flowering=weighted.mean(w = flowers,seasonday))
  summarise(avg_flowering=weighted.mean(w = flowers,julian_day))

flowering_correlations_average_bio1<-
ggplot(estimated_flowering,
       aes(y=avg_flowering, x=sites_bio1, group=year, color=as.numeric(year))
)+
  geom_point()+
  # stat_smooth(method='glm', col="grey")+
  stat_smooth(method='glm')+
  scale_color_continuous("year",type = 'viridis')+
  # geom_hline(lty="dotted",yintercept = 1)+ # equinox
  # geom_hline(lty="dotted",yintercept = 355-264)+ # solstice
  # geom_hline(lty="dotted",yintercept = 80+(365-264))+ # march equinox
  # facet_wrap(~year)+
  labs(y="Flowering day")+
  labs(x="Flowering day")+
  theme_minimal()
flowering_correlations_average_bio1

save_plot("figs/fig-date-flowering-average-and-bio1.pdf",flowering_correlations_average_bio1,base_height = 5, base_width = 7)
save_plot("figs/fig-date-flowering-average-and-bio1.png",flowering_correlations_average_bio1,base_height = 5, base_width = 7)

#### Flowering blue green colors
estimated_flowering<-
  e %>%
  group_by(site, year, latitude, longitude, altitude, sites_bio1) %>%
  # filter( julian_day< 244) %>% # separate things that flower before and after september

  # group_by(site, year, latitude, longitude, altitude, julian_day > 244) %>% # separate things that flower before and after september
  # dplyr::mutate(seasonday=ifelse(julian_day<264,julian_day+(365-264), julian_day-264 )) %>% # to exclude the Netherland?
  # summarise(avg_flowering=weighted.mean(w = flowers,seasonday))
  summarise(avg_flowering=weighted.mean(w = flowers,julian_day))


flowering_correlations_average_bio1_exp<-
  ggplot(estimated_flowering)+
  # stat_smooth(
  #   aes(y=avg_flowering, x=sites_bio1, group=year, color=avg_flowering),
  #   method='glm', col="grey", se = "lightgrey")+
  # geom_point(
  #   aes(y=avg_flowering, x=sites_bio1, group=year),
  #   size=5.2, color="black"
  # )+
  geom_point(
    aes(y=as.Date(avg_flowering,), x=sites_bio1, group=year, color=avg_flowering),
    size=5
  )+
  scale_y_date(date_breaks = "1 month",date_labels = "%b")+
  scale_colour_gradient2("FT",low = "#76AB33",mid = "lightgrey",high = "#798AE9",midpoint = 125)+
  # geom_hline(lty="dotted",yintercept = 1)+ # equinox
  # geom_hline(lty="dotted",yintercept = 355-264)+ # solstice
  # geom_hline(lty="dotted",yintercept = 80+(365-264))+ # march equinox
  # facet_wrap(~year)+
  labs(y="Flowering day")+
  xlab("Annual temperature (C)")+
  theme_minimal()
flowering_correlations_average_bio1_exp

save_plot("figs/fig-date-flowering-average-and-bio1-bluegreen.pdf",
          flowering_correlations_average_bio1_exp,base_height = 4, base_width = 5)
save_plot("figs/fig-date-flowering-average-and-bio1-bluegreen.png",
          flowering_correlations_average_bio1_exp,base_height = 4, base_width = 5)
