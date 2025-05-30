################################################################################
### Goal
### Estimate some basic metrics of ecotypes increasing in frequency across
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


#####*********************************************************************######
### Load data
load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize_flowers.rda")
e = ecotype_allyears_frequencies_long_raw_climate_popsize_flowers


#####*********************************************************************######
#####* Quantify heritability with lmer
# Check that it works
esub<-e %>%
  dplyr::filter(year==1)

REDO=F
if(REDO){
# Test LME4
library(lme4)
mod<-
  lmer(
    formula= log((freq+0.001)/startfreq) ~ 1 + (1|id),
    data = esub
  )
id<-as.data.frame(VarCorr(mod),comp="Variance")[1,"vcov"]
res<-as.data.frame(VarCorr(mod),comp="Variance")[2,"vcov"]
id/(id+res)

geth2freq<-function(freq,startfreq,id){
  mydata<-data.frame(
    freq,startfreq,id
  )
  mod<-
    lmer(
      formula= log((freq+0.001)/startfreq) ~ 1 + (1|id),
      data = mydata
    )
  id<-as.data.frame(VarCorr(mod),comp="Variance")[1,"vcov"]
  res<-as.data.frame(VarCorr(mod),comp="Variance")[2,"vcov"]
  id/(id+res)
}

# Now per garden
# myres<-
#   esub %>%
#     group_by(site) %>%
#     summarise(h2freq=geth2freq(freq,startfreq,id))

myres<-c()
for(i in unique(esub$site)){
  esubsub<-
    esub %>% dplyr::filter(site==i)
  # esubsub$y<-log((esubsub$freq+0.001)/esubsub$startfreq)
  if(length(unique(esubsub$rep))>1){
    mod<-
      lmer(
        formula= log((freq+0.001)/startfreq) ~ 1 + (1|id),
        # formula= y ~ 1 + (1|id),
        data = esubsub
      )
    id<-as.data.frame(VarCorr(mod),comp="Variance")[1,"vcov"]
    res<-as.data.frame(VarCorr(mod),comp="Variance")[2,"vcov"]
    myres<-c(myres,id/(id+res))
  }else{
    myres<-c(myres,NA)
  }
}
myva<-data.frame(site=unique(esub$site),pseudoVa=myres)
write.table(x = myva,
            file="data-intermediate/generation_1_heritability_ecotypefreq.txt",
            quote = F,row.names = F
            )
}else{
  myva<-read.table("data-intermediate/generation_1_heritability_ecotypefreq.txt",header = T)
}
#####*********************************************************************######
#####* Quantify repeatability

para<-read.table("data-intermediate/generation_1_parallelism.txt", header = T) %>%
  dplyr::filter(source=="ecotype")

head(para)
tail(para)


#####*********************************************************************######
##### Repeatability vs heritability

m<-merge(para,myva,by="site")

qplot(data=m, x=mean, y=pseudoVa)

cor.test(m$mean, m$pseudoVa, method='p')
