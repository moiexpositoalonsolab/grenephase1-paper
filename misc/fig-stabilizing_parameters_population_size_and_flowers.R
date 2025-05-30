################################################################################
### Goal
### See if Vs per site estimates extinction of populations or num. flowers

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

# load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate.rda")
# load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize.rda")

load(file="data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize_flowers.rda")
load(file="data-intermediate/ecotype_terminal_frequencies_long_raw_climate_popsize.rda")


e<-enew <- ecotype_allyears_frequencies_long_raw_climate_popsize_flowers
et<-enew <- ecotype_terminal_frequencies_long_raw_climate_popsize %>% 
  rename(year=max_year, freq=maxfreq) # Use the terminal
et[1:5,1:10]

################################################################################

################################################################################

# This is the population size created by Tati
survival<-read.csv("data/survival.csv")

survival<-
survival %>% 
  # head %>% 
  dplyr::select(site,plot,ends_with("survival")) %>% 
  dplyr::select(-ends_with("comments")) %>% 
  pivot_longer(cols = ends_with("survival"),
               # names_pattern = "X(.)_(.)",
               names_prefix = "X",
               # names_suffix = "_survival",
               names_to="year",
               values_to="survival"
               ) %>% 
  dplyr::mutate(generation=gsub("_survival","",year)) %>% 
  dplyr::mutate(generation=as.numeric(generation)) %>%   
  dplyr::mutate(survival=as.numeric(survival))%>% 
  dplyr::mutate(year=generation)  %>%  # so that then you merge with the year after survival
  dplyr::mutate(rep=plot)
survival$survival[survival$survival < 0] <- NA
survival$survival[is.na(survival$survival)] <- NA

################################################################################
# Relationship of summary simple statistics

load(file = "data-intermediate/simplesummaries-persite-allyears.rda")
load(file = "data-intermediate/simplesummaries-persite-perreplicate-allyears.rda")
sim<-simplemetricsallyears
sim<-simplemetricsallyearsperreplicate

survival2<-
  survival %>% 
  dplyr::filter(year<=3) %>% 
  dplyr::group_by(site,plot) %>% 
  dplyr::mutate(totalsurvivalyears=sum(survival,na.rm=T))

# Merge
sim %>% 
  dplyr::filter(year==1) %>% 
  merge(.,
        survival2,
        by=c("site","year","rep")
  )->
  sumsur


plot(sumsur$totalsurvivalyears ~ sumsur$var_freq)
plot(sumsur$totalsurvivalyears ~ sumsur$avg_freq)

################################################################################
# Get the relationship of Vs per site and per replicate with the survival of that replicate
survival2<-
  survival %>% 
  dplyr::filter(year<=3) %>% 
  dplyr::group_by(site,plot) %>% 
  dplyr::mutate(totalsurvivalyears=sum(survival,na.rm=T))

# Merge
vsur<-
  merge(vssitesyearsperreps,
        survival2,
        by=c("site","year","rep")
        )

vsur %>% dim
vsur %>% head
vsur$survival %>% table
vsur$totalsurvivalyears %>% table

mod<-
vsur %>% 
  dplyr::filter(year==1) %>% 
  lm(data=., 
      totalsurvivalyears ~ flowerscorectedmean * Vssite #,
      # family = binomial(link = "logit")
      )
mod %>% summary()

plot(vsur$totalsurvivalyears ~ vsur$flowerscorectedmean)
plot(vsur$totalsurvivalyears ~ vsur$Vssite)

################################################################################
# Get the relationship of Vs per site (no per replicate) with the survival of the entire site

survival3<-
  survival %>% 
  dplyr::filter(year<=3) %>% 
  dplyr::group_by(site,plot) %>% 
  dplyr::mutate(totalsurvivalyears=sum(survival,na.rm=T)) %>% 
  ungroup() %>% 
  dplyr::group_by(site, year) %>% 
  dplyr::summarise(meantotalsurvivalyears=mean(totalsurvivalyears,na.rm=T))

survival3$meantotalsurvivalyears %>% hist

# Merge
vsur<-
  merge(vssitesyears,
        survival3,
        by=c("site","year")
  )

vsur %>% dim
vsur %>% head

mod<-
  vsur %>% 
  dplyr::filter(year==1) %>% 
  lm(data=., 
     meantotalsurvivalyears ~  Vssite #,
     # family = binomial(link = "logit")
  )

plot(vsur$meantotalsurvivalyears ~ vsur$Vssite)
cor.test(vsur$meantotalsurvivalyears , vsur$Vssite)
qplot(vsur$meantotalsurvivalyears ,x= vsur$Vssite)+stat_smooth(method = 'glm')


# ################################################################################
# # Is Vs connected with population loss?
# 
# # This is the population size created by Tati
# surv<-read.csv("data/survival.csv")
# 
# 
# # This is the population size created by Tati
# popsize=read.csv(paste0(myfolder, "/data-intermediate//pop_size_estimations.csv")) %>%
#   data.frame %>% 
#   dplyr::select(-X,-longitude,-latitude, -predicted) %>% 
#   dplyr::rename(flowerspred=flowerscollected_corrected, plantspred=totalplantnumber_complete,year=generation) %>% 
#   dplyr::select(site,plot,year,flowerspred, plantspred) %>% 
#   dplyr::group_by(site,year) %>% 
#   summarise(plantspredmean=mean(plantspred), 
#             plantspredsd=sd(plantspred),
#             flowerspredmean=mean(flowerspred),
#             flowerspredsd=sd(flowerspred)
#             )
# 
# # Samples flower numbers
# load("grene/data/samples_data.rda")
# samples_tmp<-samples_data %>% 
#   dplyr::filter(year<2021) %>% 
#   dplyr::mutate(year= year-2017) %>% 
#   dplyr::rename(flowers=flowerscollected) %>% 
#   dplyr::select(site,plot,year, date,month,day,flowers) %>% 
#   # dplyr::mutate(julian_day = yday(ymd(date))) %>% 
#   dplyr::select(site, year, plot, flowers) %>% 
#   ungroup %>% 
#   dplyr::group_by(site, year) %>% 
#   dplyr::summarise(flowersmean=mean(flowers), flowerssd=sd(flowers))
# 
# 
# vs<-
#   vssitesyears %>% 
#   merge(.,popsize, by=c("site","year")) %>% 
#   merge(.,samples_tmp, by=c("site","year")) %>% 
#   mutate(Vssite=as.numeric(Vssite))
# vs
# 
# ggplot(vs,
#        aes(y=Vssite, x=plantspredmean)
#        )+
#   geom_point()+
#   stat_smooth(method='glm')+
#   geom_hline(yintercept = 0,lty='dotted')+
#   ylab("1/Vs garden")
# 
# ggplot(vs,
#        aes(y=Vssite, x=flowerspredmean)
# )+
#   geom_point()+
#   stat_smooth(method='glm')+
#   geom_hline(yintercept = 0,lty='dotted')+
#   ylab("1/Vs garden")
# 
# ggplot(vs,
#        aes(y=Vssite, x=plantspredsd)
# )+
#   geom_point()+
#   stat_smooth(method='glm')+
#   geom_hline(yintercept = 0,lty='dotted')+
#   ylab("1/Vs garden")
# 
# 
# 
# ggplot(vs,
#        aes(y=Vssite, x=flowerspredsd)
# )+
#   geom_point()+
#   stat_smooth(method='glm')+
#   geom_hline(yintercept = 0,lty='dotted')+
#   ylab("1/Vs garden")
# 
# 
# ggplot(vs,
#        aes(y=Vssite, x=flowersmean)
# )+
#   geom_point()+
#   stat_smooth(method='glm')+
#   geom_hline(yintercept = 0,lty='dotted')+
#   ylab("1/Vs garden")
# 
# 
# ggplot(vs,
#        aes(y=Vssite, x=flowerssd/flowersmean)
# )+
#   geom_point()+
#   stat_smooth(method='glm')+
#   geom_hline(yintercept = 0,lty='dotted')+
#   ylab("1/Vs garden")

