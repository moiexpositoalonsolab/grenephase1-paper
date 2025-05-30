################################################################################
### Goal
### Compare population size with the parallelism
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

## Climate 
load("grene/data/worldclim_sitesdata.rda")
load("grene/data/locations_data.rda")

### Parallelism
para<-read.table("data-intermediate/generation_1_parallelism.txt",header = T) %>% 
  dplyr::filter(source=="snp")
  

### Calculate quadratic terms
quad<-data.frame(site=unique(popsize$site),
                 a=NA,
                 b=NA,
                 c=NA
                 )
counter=1
for(counter in 1:length(unique(popsize$site))){
  i=quad[counter,1]
  message(i)
  tmpsub<-dplyr::filter(popsize,site==i)
  if(length(unique(tmpsub$generation))==1){
    # only one generation
    tmpmod<-coefficients(lm(data=tmpsub,flowerscollected_corrected ~ 1))
  }else{
    # 2 or 3 generations
  tmpmod<-tryCatch(
    tmpmod<-coefficients(lm(data=tmpsub,flowerscollected_corrected ~ poly(generation,2))),
    error= function(e) tmpmod<-coefficients(lm(data=tmpsub,flowerscollected_corrected ~ poly(generation,1)))
  )
  }
  if(length(tmpmod)==3){
    a=tmpmod[1]
    b=tmpmod[2]
    c=tmpmod[3]
  }else if(length(tmpmod)==2){
    a=tmpmod[1]
    b=tmpmod[2]
    c=0
  }else if(length(tmpmod)==1){
    a=tmpmod[1]
    b=0
    c=0
  }else{
    stop("something happened")
  }
  message(a," ", b, " ", c)
  quad[counter,"a"]<-a
  quad[counter,"b"]<-b
  quad[counter,"c"]<-c
}      

quad %>% head

################################################################################
#### Merge
m<-
  merge(popsize,worldclim_sitesdata,by="site") %>% 
  merge(.,para, by="site") %>% 
  merge(.,quad,by="site")

#### Plot
quadraticplot<-
ggplot(m)+
  geom_point(aes(y=c,x=mean,size=a,color=bio1))+
  scale_color_gradientn(colors = redblue)+
  labs(y="Pop. quadratic term (c)", 
       x="Evolutionary parallelism", 
       size="Pop. intercept (a)",
       color="Temperature (C)"
       )
quadraticplot

ggsave("figs/fig-quadratic-population-term-cliamte-evo-parallelism.pdf",
       quadraticplot,
       height = 5, width = 6, units="in")
ggsave("figs/fig-quadratic-population-term-cliamte-evo-parallelism.png",
       quadraticplot,
       height = 5, width = 6, units="in")
