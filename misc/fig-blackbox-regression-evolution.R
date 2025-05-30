################################################################################
### Goal
### Regression 

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
# et<-enew <- ecotype_terminal_frequencies_long_raw_climate_popsize %>% 
#   rename(year=max_year, freq=maxfreq) # Use the terminal
# et[1:5,1:10]

################################################################################
################################################################################
# This is the population size created by Tati
# survival<-read.csv("data/survival.csv")
load("data/flowers_survival_long.rda")


################################################################################
# Random forest of ecotype frequency based on climate of origina and climate of site

train<-e %>% 
  dplyr::select(id, flowers,freq,year,starts_with("sites_"),starts_with("ecotypes_"), latitude,longitude)  %>% 
  group_by(latitude,longitude,id, year, 
           sites_bio1, ecotypes_bio1,
           sites_bio2, ecotypes_bio2,
           sites_bio3, ecotypes_bio3,
           sites_bio4, ecotypes_bio4,
           sites_bio5, ecotypes_bio5,
           sites_bio6, ecotypes_bio6,
           sites_bio7, ecotypes_bio7,
           sites_bio8, ecotypes_bio8,
           sites_bio9, ecotypes_bio9,
           sites_bio10, ecotypes_bio10,
           sites_bio11, ecotypes_bio11,
           sites_bio12, ecotypes_bio12,
           sites_bio13, ecotypes_bio13,
           sites_bio14, ecotypes_bio14,
           sites_bio15, ecotypes_bio15,
           sites_bio16, ecotypes_bio16,
           sites_bio17, ecotypes_bio17,
           sites_bio18, ecotypes_bio18,
           sites_bio19, ecotypes_bio19,
  ) %>% 
  summarise(freq=weighted.mean(freq,flowers),
            totalflowers=sum(flowers)
  ) %>% 
  mutate(individuals=freq*totalflowers) %>% 
  na.omit
dim(train)

# Train the lm
train %>% 
  ungroup() %>% 
  dplyr::select(individuals,freq,starts_with("sites_"),starts_with("ecotypes_"), year ) %>% 
  dplyr::mutate(log10_ind=log10(feq+1) ,
                dbio1=(sites_bio1 - ecotypes_bio1)^2,
                dbio2=(sites_bio2 - ecotypes_bio2)^2,
                dbio3=(sites_bio3 - ecotypes_bio3)^2,
                dbio4=(sites_bio4 - ecotypes_bio4)^2,
                dbio5=(sites_bio5 - ecotypes_bio5)^2,
                dbio6=(sites_bio6 - ecotypes_bio6)^2,
                dbio7=(sites_bio7 - ecotypes_bio7)^2,
                dbio8=(sites_bio8 - ecotypes_bio8)^2,
                dbio9=(sites_bio9 - ecotypes_bio9)^2,
                dbio10=(sites_bio10 - ecotypes_bio10)^2,
                dbio11=(sites_bio11 - ecotypes_bio11)^2,
                dbio12=(sites_bio12 - ecotypes_bio12)^2,
                dbio13=(sites_bio13 - ecotypes_bio13)^2,
                dbio14=(sites_bio14 - ecotypes_bio14)^2,
                dbio15=(sites_bio15 - ecotypes_bio15)^2,
                dbio16=(sites_bio16 - ecotypes_bio16)^2,
                dbio17=(sites_bio17 - ecotypes_bio17)^2,
                dbio18=(sites_bio18 - ecotypes_bio18)^2,
                dbio19=(sites_bio19 - ecotypes_bio19)^2
                ) %>% 
  dplyr::select(log10_ind, starts_with("dbio"), starts_with("sites_")) %>% 
  lm(data=.,log10_ind~.)->
  # randomForest(data=.,log10(individuals+1)~.)->
  mainlm

mainlm %>% summary
plot(log10(train$individuals+1) ~ predict(mainlm))
cor(log10(train$individuals+1), predict(mainlm))

################################################################################
# RF: Train the random forest with bioclim of sites adn distances
train %>% 
  ungroup() %>% 
  dplyr::select(individuals,freq,starts_with("sites_"),starts_with("ecotypes_"), year ) %>% 
  dplyr::filter(year==2) %>% 
  dplyr::mutate(
                # y=log10(individuals+1) ,
                y=log10(freq+0.001),
                dbio1=(sites_bio1 - ecotypes_bio1)^2,
                dbio2=(sites_bio2 - ecotypes_bio2)^2,
                dbio3=(sites_bio3 - ecotypes_bio3)^2,
                dbio4=(sites_bio4 - ecotypes_bio4)^2,
                dbio5=(sites_bio5 - ecotypes_bio5)^2,
                dbio6=(sites_bio6 - ecotypes_bio6)^2,
                dbio7=(sites_bio7 - ecotypes_bio7)^2,
                dbio8=(sites_bio8 - ecotypes_bio8)^2,
                dbio9=(sites_bio9 - ecotypes_bio9)^2,
                dbio10=(sites_bio10 - ecotypes_bio10)^2,
                dbio11=(sites_bio11 - ecotypes_bio11)^2,
                dbio12=(sites_bio12 - ecotypes_bio12)^2,
                dbio13=(sites_bio13 - ecotypes_bio13)^2,
                dbio14=(sites_bio14 - ecotypes_bio14)^2,
                dbio15=(sites_bio15 - ecotypes_bio15)^2,
                dbio16=(sites_bio16 - ecotypes_bio16)^2,
                dbio17=(sites_bio17 - ecotypes_bio17)^2,
                dbio18=(sites_bio18 - ecotypes_bio18)^2,
                dbio19=(sites_bio19 - ecotypes_bio19)^2
  ) %>% 
  dplyr::select(y, starts_with("dbio"), starts_with("sites_")) %>% 
  # lm(data=.,log10_ind~.)->
  randomForest(data=.,y~.)->
  mainrf

mainrf

plot(mainrf$y ~ mainrf$predicted)
cor(mainrf$y , mainrf$predicted)

# Save random forest
rf_ecofreq_year2_bioclim_bioclimdistances<-mainrf
sink(file = "models/rf_ecofreq_year2_bioclim_bioclimdistances.info")
print(rf_ecofreq_year2_bioclim_bioclimdistances)
print(rf_ecofreq_year2_bioclim_bioclimdistances$importance)
sink()
save(file = "models/rf_ecofreq_year2_bioclim_bioclimdistances.rda",rf_ecofreq_year2_bioclim_bioclimdistances)

# Train the random forest with years
train %>% 
  ungroup() %>% 
  dplyr::select(individuals,freq,starts_with("sites_"),starts_with("ecotypes_"), year ) %>% 
  # dplyr::filter(year==2) %>% 
  dplyr::mutate(
    # y=log10(individuals+1) ,
    y=log10(freq+0.001),
    dbio1=(sites_bio1 - ecotypes_bio1)^2,
    dbio2=(sites_bio2 - ecotypes_bio2)^2,
    dbio3=(sites_bio3 - ecotypes_bio3)^2,
    dbio4=(sites_bio4 - ecotypes_bio4)^2,
    dbio5=(sites_bio5 - ecotypes_bio5)^2,
    dbio6=(sites_bio6 - ecotypes_bio6)^2,
    dbio7=(sites_bio7 - ecotypes_bio7)^2,
    dbio8=(sites_bio8 - ecotypes_bio8)^2,
    dbio9=(sites_bio9 - ecotypes_bio9)^2,
    dbio10=(sites_bio10 - ecotypes_bio10)^2,
    dbio11=(sites_bio11 - ecotypes_bio11)^2,
    dbio12=(sites_bio12 - ecotypes_bio12)^2,
    dbio13=(sites_bio13 - ecotypes_bio13)^2,
    dbio14=(sites_bio14 - ecotypes_bio14)^2,
    dbio15=(sites_bio15 - ecotypes_bio15)^2,
    dbio16=(sites_bio16 - ecotypes_bio16)^2,
    dbio17=(sites_bio17 - ecotypes_bio17)^2,
    dbio18=(sites_bio18 - ecotypes_bio18)^2,
    dbio19=(sites_bio19 - ecotypes_bio19)^2
  ) %>% 
  dplyr::select(y, starts_with("dbio"), starts_with("sites_"),year) %>% 
  # lm(data=.,log10_ind~.)->
  randomForest(data=.,y~.)->
  mainyearrf
mainyearrf
plot(mainyearrf$y ~ mainyearrf$predicted)
cor(mainyearrf$y , mainyearrf$predicted)

################################################################################
# Train the random forest
train %>% 
  ungroup() %>% 
  dplyr::select(individuals,freq,starts_with("sites_"),starts_with("ecotypes_"), year ) %>% 
  dplyr::filter(year==2) %>% 
  dplyr::mutate(
    # y=log10(individuals+1) ,
    y=log10(freq+0.001),
    dbio1=(sites_bio1 - ecotypes_bio1)^2,
    dbio2=(sites_bio2 - ecotypes_bio2)^2,
    dbio3=(sites_bio3 - ecotypes_bio3)^2,
    dbio4=(sites_bio4 - ecotypes_bio4)^2,
    dbio5=(sites_bio5 - ecotypes_bio5)^2,
    dbio6=(sites_bio6 - ecotypes_bio6)^2,
    dbio7=(sites_bio7 - ecotypes_bio7)^2,
    dbio8=(sites_bio8 - ecotypes_bio8)^2,
    dbio9=(sites_bio9 - ecotypes_bio9)^2,
    dbio10=(sites_bio10 - ecotypes_bio10)^2,
    dbio11=(sites_bio11 - ecotypes_bio11)^2,
    dbio12=(sites_bio12 - ecotypes_bio12)^2,
    dbio13=(sites_bio13 - ecotypes_bio13)^2,
    dbio14=(sites_bio14 - ecotypes_bio14)^2,
    dbio15=(sites_bio15 - ecotypes_bio15)^2,
    dbio16=(sites_bio16 - ecotypes_bio16)^2,
    dbio17=(sites_bio17 - ecotypes_bio17)^2,
    dbio18=(sites_bio18 - ecotypes_bio18)^2,
    dbio19=(sites_bio19 - ecotypes_bio19)^2
  ) %>% 
  dplyr::select(y, starts_with("dbio"), starts_with("sites_"),starts_with("ecotypes_")) %>% 
  # lm(data=.,log10_ind~.)->
  randomForest(data=.,y~.)->
  mainrf

mainrf

plot(mainrf$y ~ mainrf$predicted)
cor(mainrf$y , mainrf$predicted)

################################################################################
# Without distances?
train %>% 
  ungroup() %>% 
  dplyr::select(individuals,freq,starts_with("sites_"),starts_with("ecotypes_"), year ) %>% 
  dplyr::filter(year==2) %>% 
  dplyr::mutate(
    # y=log10(individuals+1) ,
    y=log10(freq+0.001)
  ) %>% 
  dplyr::select(y, starts_with("sites_"),starts_with("ecotypes_")) %>% 
  # lm(data=.,log10_ind~.)->
  randomForest(data=.,y~.)->
  mainrf

mainrf
plot(mainrf$y ~ mainrf$predicted)
cor(mainrf$y , mainrf$predicted)


# Save random forest
rf_ecofreq_year2_bioclimsites_bioclimecotypes<-mainrf
sink(file = "models/rf_ecofreq_year2_bioclimsites_bioclimecotypes.info")
print(mainrf)
print(mainrf$importance)
sink()
save(file = "models/rf_ecofreq_year2_bioclimsites_bioclimecotypes.rda",
     rf_ecofreq_year2_bioclimsites_bioclimecotypes)
