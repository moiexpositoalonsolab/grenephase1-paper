################################################################################
### Goal
### Estimate some basic metrics of per sites summary statistics
### sites, means, variances, optimal places

################################################################################
library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(cowplot)

################################################################################
#### Location
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
# Load the long format with ecotype frequencies and climates
load("data-intermediate/ecotype_terminal_frequencies_long_raw_climate.rda")
load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize_flowers.rda")
load("grene/data/worldclim_sitesdata.rda")

e = ecotype_allyears_frequencies_long_raw_climate_popsize_flowers
et = ecotype_terminal_frequencies_long_raw_climate
head(e)

################################################################################
# There are 3 things in these ecotypes, first the distance

# Create a subset to get things easier
es= e %>% 
  dplyr::select(id,rep,year,site,longitude,latitude, ecotypes_bio1, sites_bio1, freq, startfreq) 
head(es)

# Get some general info, like what is the per ecotype average

simplemetricsallyears <-
  es %>% 
  dplyr::group_by(site) %>% 
  dplyr::mutate(relative= freq/sum(freq)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(site, year, sites_bio1) %>% 
  dplyr::summarise(avg_freq = mean(freq),
            cv_freq = var(freq)/mean(freq),
            var_freq = var(freq),
            opt_clim = weighted.mean(ecotypes_bio1,relative,na.rm=T)
  ) 
simplemetricsallyears$opt_clim
simplemetricsallyears %>% head
save(file = "data-intermediate/simplesummaries-persite-allyears.rda",
     simplemetricsallyears)

simplemetricsallyearsperreplicate <- 
  es %>% 
  dplyr::group_by(site) %>% 
  dplyr::mutate(relative= freq/sum(freq)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(site, year, rep, sites_bio1) %>% 
  dplyr::summarise(avg_freq = mean(freq),
                   cv_freq = var(freq)/mean(freq),
                   var_freq = var(freq),
                   opt_clim = weighted.mean(ecotypes_bio1,relative,na.rm=T)
  ) 
simplemetricsallyearsperreplicate %>% head(2)
simplemetricsallyearsperreplicate$opt_clim
save(file = "data-intermediate/simplesummaries-persite-perreplicate-allyears.rda",
     simplemetricsallyearsperreplicate)

################################################################################

ggplot(simplemetricsallyearsperreplicate)+
  geom_point(aes(y=avg_freq, x=sites_bio1))

ggplot(simplemetricsallyears)+
  geom_point(aes(y=var_freq, x=sites_bio1))

ggplot(simplemetricsallyearsperreplicate)+
  geom_point(aes(y=opt_clim, x=sites_bio1, color=year))+
  stat_smooth(aes(y=opt_clim, x=sites_bio1, color=year),method='glm')+
  geom_abline(slope = 1,intercept = 0)+
  ylab("Temperature of optimal accession (C)")+xlab("Temperature of garden (C)")+
  xlim(c(4,21))+ylim(c(4,21))+
  geom_hline( yintercept = min(es$ecotypes_bio1,na.rm = T),ltd='dotted')+
  geom_hline( yintercept = max(es$ecotypes_bio1,na.rm = T),ltd='dotted')+
  theme_minimal()

ggplot(simplemetricsallyearsperreplicate)+
  geom_histogram(aes(x=opt_clim),color="lightgrey")+
  geom_histogram(aes(x=sites_bio1),color="black")+
  theme_minimal()


ggplot(simplemetricsallyearsperreplicate)+
  geom_histogram(aes(x=opt_clim-sites_bio1))+
  theme_minimal()


test <-
  e %>% 
  dplyr::group_by(site,rep, latitude) %>% 
  dplyr::mutate(relative= freq/sum(freq)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(site, year, sites_bio1,latitude) %>% 
  dplyr::summarise(avg_freq = mean(freq),
                   cv_freq = var(freq)/mean(freq),
                   var_freq = var(freq),
                   opt_clim = weighted.mean(sites_bio1-ecotypes_bio1,relative,na.rm=T),
                   dif_clim = mean(sites_bio1-ecotypes_bio1,na.rm=T)
  ) 

ggplot(test)+
  geom_histogram( aes(x=opt_clim - dif_clim))
ggplot(test)+
  geom_histogram( aes(x=opt_clim), color="black")+
  geom_histogram(aes(x=dif_clim), color="green")

weighted_optimal_ecotypes_bio1<-
ggplot(test)+
  geom_point( aes(y=opt_clim, x=dif_clim, color=latitude))+
  geom_abline(intercept = 0,slope = 1)+
  xlab("Average distance ecotype-to-garden (Temp C)")+
  ylab("Weighted frequency distance ecotype-to-garden (Temp C)")+
  labs(color="Latitude (N)")+
  theme_minimal()
weighted_optimal_ecotypes_bio1
save_plot("figs/fig-simple-metrics-ecotype-bio1-origin-garden-weighted.pdf",
          weighted_optimal_ecotypes_bio1,
          base_width = 6,base_height = 5)
save_plot("figs/fig-simple-metrics-ecotype-bio1-origin-garden-weighted.png",
          weighted_optimal_ecotypes_bio1,
          base_width = 6,base_height = 5)
