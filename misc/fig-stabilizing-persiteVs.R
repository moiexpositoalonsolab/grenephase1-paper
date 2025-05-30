################################################################################
### Goal
### Estimate Vs per location 

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
# Get Vs from sites

source("function-regression-utilities.R")
source("function-colors-grenenet.R")


myregression<-function(y,x){
  mylm <- lm(y ~ x)
  p = format(coefficients(summary(mylm))[2, 4] )
             # , digits = 3, 
             # scientific = TRUE)
  r2 = round(summary(mylm)$r.squared, digits = 3)
  b = format(coefficients(summary(mylm))[2, 1] )
             # , digits = 3, 
             # scientific=TRUE)
  return(b)
}


vssitesyearsperreps<-
  e %>% 
  group_by(site,year,rep,
           sites_bio1, longitude,latitude) %>% 
  dplyr::mutate(e2=(sites_bio1-ecotypes_bio1)^2,
                logp1p0=log(freq/startfreq)
  ) %>% 
  dplyr::filter(e2<400) %>% 
  dplyr::filter(freq!=0) %>% 
  summarise(Vssite= myregression(logp1p0, e2 ),
            flowersmean=mean(flowers),
            flowerscorectedmean=mean(flowerscollected_corrected),
            plantsmean=mean(totalplantnumber_complete),
  ) %>% 
  # merge(.,
  #       dplyr::select(e,site,longitude,latitude, starts_with("sites_"),
  #                     starts_with("flower"), starts_with("plant")), 
  #       by="site") %>% 
  mutate(Vssite=as.numeric(Vssite))

save(file = "data-intermediate/vssitesyearsperreps.rda",vssitesyearsperreps)

vssitesyears<-
e %>% 
  group_by(site,year,
           sites_bio1, longitude,latitude) %>% 
  dplyr::mutate(e2=(sites_bio1-ecotypes_bio1)^2,
                logp1p0=log(freq/startfreq)
                ) %>% 
  dplyr::filter(e2<400) %>% 
  dplyr::filter(freq!=0) %>% 
  summarise(Vssite= myregression(logp1p0, e2 ),
            flowersmean=mean(flowers),
            flowerscorectedmean=mean(flowerscollected_corrected),
            plantsmean=mean(totalplantnumber_complete),
            ) %>% 
  # merge(.,
  #       dplyr::select(e,site,longitude,latitude, starts_with("sites_"),
  #                     starts_with("flower"), starts_with("plant")), 
  #       by="site") %>% 
  mutate(Vssite=as.numeric(Vssite))

save(file = "data-intermediate/vssitesyears.rda",vssitesyears)
write.csv(file = "data-intermediate/stabilizing_selection_simple_vspersite.csv",
         x=vssitesyears,
         row.names=F)

vssites<-
  e %>% 
  group_by(site) %>% 
  summarise(Vssite= myregression(log(freq+0.001/startfreq), (sites_bio1-ecotypes_bio1)^2 )) %>% 
  # merge(.,
  #       dplyr::select(e,site,longitude,latitude, starts_with("sites_"),starts_with("flower"), starts_with("plant")), 
  #       by="site")%>% 
  mutate(Vssite=as.numeric(Vssite))


vssitesterminal<-
  et %>% 
  group_by(site,
           sites_bio1, longitude,latitude,
  ) %>% 
  dplyr::mutate(e2=(sites_bio1-ecotypes_bio1)^2,
                logp1p0=log(freq/startfreq)
  ) %>% 
  dplyr::filter(e2<400) %>% 
  dplyr::filter(freq!=0) %>% 
  summarise(Vssite= myregression(logp1p0, e2 ),
            flowerscorectedmean=mean(flowerscollected_corrected),
            plantsmean=mean(totalplantnumber_complete),
            flowerssum=sum(flowerscollected_corrected),
            plantssum=sum(totalplantnumber_complete),
  ) %>% 
  mutate(Vssite=as.numeric(Vssite))
save(file = "data-intermediate/vssitesterminal.rda",vssitesterminal)
  


################################################################################
# Several visualizations of the main pattern

# Filtering out 
persitestabilizing<-
e %>% 
  dplyr::filter(freq!=0) %>% 
  ggplot(., aes(y=log(freq/startfreq),
                x= (sites_bio1-ecotypes_bio1)^2,
                group=site
  ))+
  geom_point(color="grey",alpha=0.25)+
  stat_smooth(method='glm', aes(group=site,color=as.numeric(sites_bio1)), se=F)+
  scale_color_gradientn("Temp. (C)", colors=brewer.pal(5,"Reds"))+
  xlim(c(0,400))+
  ylab("log(p1/p0)")+xlab("Temperature distance (e^2)")+
  theme_minimal()
save_plot("figs/fig-stabilizing-simple-Vssite-allyears-bio1.pdf",
          persitestabilizing,
          base_width = 6,base_height = 5)

save_plot("figs/fig-stabilizing-simple-Vssite-allyears-bio1.png",
          persitestabilizing,
          base_width = 6,base_height = 5)


e %>% 
  dplyr::filter(freq!=0) %>% 
  ggplot(., aes(y=log(freq/startfreq),
                x= (sites_bio1-ecotypes_bio1)^2,
                group=site
  ))+
  geom_point(color="grey",alpha=0.25)+
  stat_smooth(method='glm', aes(group=site,color=as.numeric(sites_bio1)), se=F)+
  scale_color_gradientn("Temp. (C)", colors=brewer.pal(5,"Reds"))+
  xlim(c(200,650))+
  ylab("log(p1/p0)")+xlab("Temperature distance (e^2)")+
  theme_minimal()


et %>% 
  dplyr::filter(freq!=0) %>% 
  ggplot(., aes(y=log(freq/startfreq),
              x= (sites_bio1-ecotypes_bio1)^2,
              group=site
              ))+
  geom_point(color="grey",alpha=0.25)+
  stat_smooth(method='glm', aes(group=site,color=as.numeric(sites_bio1)), se=F)+
  scale_color_gradientn("Temp. (C)", colors=brewer.pal(5,"Reds"))+
  ylab("log(p1/p0)")+xlab("Temperature distance (e^2)")+
  theme_minimal()
  
persitestabilizingterminal<-
et %>% 
  dplyr::filter(freq!=0) %>% 
  ggplot(., aes(y=log(freq/startfreq),
                x= (sites_bio1-ecotypes_bio1)^2,
                group=site
  ))+ 
  geom_point(color="grey",alpha=0.25)+
  stat_smooth(method='glm', aes(group=site,color=as.numeric(sites_bio1)), se=F)+
  scale_color_gradientn("Temp. (C)", colors=brewer.pal(5,"Reds"))+
  ylab("log(p1/p0)")+xlab("Temperature distance (e^2)")+
  xlim(c(0,400))+
  theme_minimal()

save_plot("figs/fig-stabilizing-simple-Vssite-terminal-bio1.pdf",
          persitestabilizingterminal,
          base_width = 6,base_height = 5)

save_plot("figs/fig-stabilizing-simple-Vssite-terminal-bio1.png",
          persitestabilizingterminal,
          base_width = 6,base_height = 5)



################################################################################
# Several visualizations

# Plot on 
ggplot(vssitesyears)+
  geom_point(aes(y=Vssite,x=longitude))
ggplot(vssitesterminal)+
  geom_point(aes(y=Vssite,x=longitude))
# Plot on 


# Plot the strength of selection against number of flowers

  
vssitesterminal %>% 
  # dplyr::filter(longitude > -15) %>% 
ggplot()+
  geom_point(aes(y=Vssite,x=latitude))+
  stat_smooth(aes(y=Vssite,x=latitude),method = "glm",color="grey")+
  stat_smooth(aes(y=Vssite,x=latitude),method = "glm", formula=y~poly(x,2), color="grey")+
  xlab("Latitude (N)")+ylab("-1/Vs per site")+
  theme_minimal()->
  vsterminal_latitude
vsterminal_latitude
save_plot("figs/fig-stabilizing-simple-Vssite-terminalyear-and-latitude.pdf",
          vsterminal_latitude,
          base_width = 9,base_height = 9)
save_plot("figs/fig-stabilizing-simple-Vssite-terminalyear-and-latitude.png",
          vsterminal_latitude,
          base_width = 9,base_height = 9)

vssitesterminal %>% 
  # dplyr::filter(longitude > -15) %>% 
ggplot(.)+
  geom_point(aes(y=Vssite,x=flowerscorectedmean))+
  # stat_smooth(aes(y=Vssite,x=flowerscorectedmean),method = "glm")+
  stat_smooth(aes(y=Vssite,x=flowerscorectedmean),method = "glm", formula=y~poly(x,2),color="grey")+
  xlab("Num. flowers")+ylab("-1/Vs per site")+
  theme_minimal()->
  vsterminal_flowers
vsterminal_flowers
save_plot("figs/fig-stabilizing-simple-Vssite-terminalyear-and-meanflowers.pdf",
          vsterminal_flowers,
          base_width = 9,base_height = 9)
save_plot("figs/fig-stabilizing-simple-Vssite-terminalyear-and-meanflowers.png",
          vsterminal_flowers,
          base_width = 9,base_height = 9)

vssitesyears %>% 
dplyr::filter(longitude > -15) %>% 
ggplot(.)+
  geom_point(aes(y=Vssite,x=latitude,color=as.numeric(year)))+
  stat_smooth(aes(y=Vssite,x=latitude,color=as.numeric(year), group=year),method = "glm", formula=y~poly(x,2))+
  scale_color_continuous("year",type = 'viridis')+
  xlab("Latitude (N)")+ylab("-1/Vs per site")+
  theme_minimal()->
  vsterminal_flowersperyear
vsterminal_flowersperyear

save_plot("figs/fig-stabilizing-simple-Vssite-allyears-and-latitude.pdf",
          vsterminal_flowersperyear,
          base_width = 9,base_height = 9)
save_plot("figs/fig-stabilizing-simple-Vssite-allyears-and-latitude.png",
          vsterminal_flowersperyear,
          base_width = 9,base_height = 9)

vssitesyears %>% 
  dplyr::filter(longitude > -15) %>% 
ggplot(.)+
  geom_point(aes(y=Vssite,x=flowerscorectedmean,color=as.numeric(year)))+
  stat_smooth(aes(y=Vssite,x=flowerscorectedmean,color=as.numeric(year),group=year),method = "glm", formula=y~poly(x,2))+
  scale_color_continuous("year",type = 'viridis')+
  xlab("Num. flowers")+ylab("-1/Vs per site")+
  theme_minimal()->
  vsterminal_latitudeperyear
vsterminal_latitudeperyear
save_plot("figs/fig-stabilizing-simple-Vssite-allyears-and-meanflowers.png",vsterminal_latitudeperyear,
          base_width = 9,base_height = 9)
save_plot("figs/fig-stabilizing-simple-Vssite-allyears-and-meanflowers.png",vsterminal_latitudeperyear,
          base_width = 9,base_height = 9)


################################################################################
#### Visualize how well we could predict trajectories

SITETEST=4

train<-
  et %>% 
  dplyr::filter(freq!=0) %>%
  dplyr::filter(site!=SITETEST) %>% 
  mutate(y=log(freq/startfreq),
       x= (sites_bio1-ecotypes_bio1)^2) %>% 
  dplyr::filter(x<400) 

test<-
  et %>% 
  dplyr::filter(freq!=0) %>%
  dplyr::filter(site==SITETEST) %>% 
  mutate(y=log(freq/startfreq),
       x= (sites_bio1-ecotypes_bio1)^2)%>% 
  dplyr::filter(x<400) 

lmod<-
  lm(data=train, y~x)
  
lpred<-
  predict(lmod,test)

mytest<-
  data.frame(y=test$y,
             ypred=lpred)

plot(test$y ~ lpred )
cor.test(test$y ,lpred , method='s')

sit4predictability<-
ggplot(mytest,
       aes(y=(y),x=(ypred)))+
  geom_point(alpha=0.5)+
  # xlim(-2.7,-1.6)+
  stat_smooth(method='glm',
              color="#41AB5D")+
  labs(x="predicted log(p1/p0)", y="observed log(p1/p0)")+
  theme_minimal()
sit4predictability

save_plot(
  paste0("figs/fig-stabilizing-simple-Vssite-terminal-predictability-site",SITETEST,".pdf"),
          sit4predictability,
          base_width = 5,base_height = 5)
save_plot(
  paste0("figs/fig-stabilizing-simple-Vssite-terminal-predictability-site",SITETEST,".png"),
          sit4predictability,
          base_width = 5,base_height = 5)


