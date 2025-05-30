################################################################################
### Goal
### Classifier of population composition

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

load(file="data-intermediate/ecotype_allyears_frequencies_long_raw.rda")
load(file="data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize_flowers.rda")
load(file="data-intermediate/ecotype_terminal_frequencies_long_raw_climate_popsize.rda")


e<-enew <- ecotype_allyears_frequencies_long_raw_climate_popsize_flowers
# et<-enew <- ecotype_terminal_frequencies_long_raw_climate_popsize %>% 
#   rename(year=max_year, freq=maxfreq) # Use the terminal
# et[1:5,1:10]

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

# Train the random forest
train %>% 
  ungroup() %>% 
  dplyr::select(individuals,starts_with("sites_"),starts_with("ecotypes_"), year ) %>% 
  randomForest(data=.,log10(individuals+1)~.)->
  mainrf

# Study predictrability random forest
mainrf
importance(mainrf) %>% data.frame %>% arrange(IncNodePurity)

# Make prediction to the dataset
train$predicted<-predict(mainrf)
plot(log10(train$individuals+1)~train$predicted) 
lm(log10(train$individuals+1)~train$predicted) %>% summary

blackbox_predictor<-
ggplot(train)+
  geom_hex(aes(y=log10(individuals+1), x= predicted))+
  ylab("log10 (individuals per ecotype)")+
  xlab("Predicted")+
  scale_fill_gradientn("# Obs",colors=rev(brewer.pal(9,"Greys")[-1]), trans="log10")+
  geom_abline(slope = 1,intercept = 0, lty="dotted")+
  theme_minimal()
blackbox_predictor

save_plot("figs/fig-blackbox-predict-ecotype-individuals-from-climate.pdf",
          blackbox_predictor, base_height = 4,base_width = 5)
save_plot("figs/fig-blackbox-predict-ecotype-individuals-from-climate.png",
          blackbox_predictor, base_height = 4,base_width = 5)


################################################################################

# Get present data
load("grene/data/worldclim_ecotypesdata.rda")
load("grene/data/worldclim_sitesdata.rda")
load("grene/data/locations_data.rda")

colnames(worldclim_sitesdata)[2:20]<-paste0("sites_",colnames(worldclim_sitesdata)[2:20])

worldclim_ecotypesdata$id<-worldclim_ecotypesdata$ecotypeid
colnames(worldclim_ecotypesdata)[2:20]<-paste0("ecotypes_",colnames(worldclim_ecotypesdata)[2:20])

### Predictions for the future
load("data-intermediate/ecotypes_future_climate.rda")
load("data-intermediate/sites_future_climate.rda")

myid<-ecotypes_future_climate$id
mysites<-sites_future_climate$site

ecotypes_future_climate<-ecotypes_future_climate %>% data.frame() %>%  dplyr::select(-id)
sites_future_climate<-sites_future_climate %>% data.frame() %>%  dplyr::select(-site)

colnames(ecotypes_future_climate)<-paste0("ecotypes_",colnames(ecotypes_future_climate))
colnames(sites_future_climate)<-paste0("sites_",colnames(sites_future_climate))

ecotypes_future_climate$id<-myid
sites_future_climate$site<-mysites

dataforpredictionpresent<-
  expand.grid(unique(myid), unique(mysites)) %>% 
  rename(id=Var1, site=Var2) %>% 
  merge(.,worldclim_ecotypesdata,by="id",all.x=T) %>% 
  merge(.,worldclim_sitesdata,by="site",all.x=T) %>% 
  mutate(year=1)

dataforprediction<-
  expand.grid(unique(myid), unique(mysites)) %>% 
  rename(id=Var1, site=Var2) %>% 
  merge(.,ecotypes_future_climate,by="id",all.x=T) %>% 
  merge(.,sites_future_climate,by="site",all.x=T) %>% 
  mutate(year=1)

presentprediction<-predict(object = mainrf, dataforpredictionpresent)
futureprediction<-predict(object = mainrf, dataforprediction)

presentprediction %>% hist
futureprediction %>% hist
length(presentprediction)
length(futureprediction)

(futureprediction-presentprediction) %>% hist

# look at predictions per site

pres<-dataforpredictionpresent
pres$presentprediction<-presentprediction

futs<-dataforprediction
futs$futureprediction<-futureprediction

head(10^pres$presentprediction)
hist(10^pres$presentprediction)
hist(10^futs$futureprediction)


# Explore whether sites have a total number of flowers different in space
toplot<-
  data.frame(site=pres$site, present=10^pres$presentprediction, future=10^futs$futureprediction) %>% 
  merge(.,locations_data,by="site",all.x=T) %>% 
  group_by(site,longitude,latitude) %>% 
  summarise(sumflowerp=sum(present, na.rm=T),
            sumflowerf=sum(future, na.rm=T)
            )
ggplot(toplot)+
  geom_point(aes(x=longitude,y=latitude,color=sumflowerp))+
  xlim(c(-15,50))+
  scale_color_gradientn(colors=brewer.pal(9,"RdBu"))+
  theme_minimal()
ggplot(toplot)+
  geom_point(aes(x=longitude,y=latitude,color=sumflowerf))+
  xlim(c(-15,50))+
  scale_color_gradientn(colors=brewer.pal(9,"RdBu"))+
  theme_minimal()

ggplot(toplot)+
  geom_point(aes(x=longitude,y=latitude,color=sumflowerf-sumflowerp))+
  xlim(c(-15,50))+
  scale_color_gradientn(colors=brewer.pal(9,"RdBu"))+
  theme_minimal()
