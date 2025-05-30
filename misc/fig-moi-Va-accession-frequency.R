## Goal
## Try to quantify Va(w) based on heritability of accession frequency

library(dplyr)
library(ggplot2)
library(RColorBrewer)
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)
source("function-colors-grenenet.R")

# prefix <- "FUTURE_PREDICTION_"
# PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
# PATH <-  "~/Shareddrives/MOI-LAB/PROJECTS/grenenet//GRENE-net_PHASE1/grenephase1-analyses/"
# setwd(PATH)

#####*********************************************************************######
# Predictability as Va (w)
# save(file = "data-intermediate/simplemetrics_persite.rda",simplemetrics_persite)


# load("data-intermediate/ecotype_terminal_frequencies_long_raw_climate.rda")
load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize_flowers.rda")

e = ecotype_allyears_frequencies_long_raw_climate_popsize_flowers

# Check that it works
esub<-e %>%
  dplyr::filter(year==1)

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
      formula= (freq-startfreq) ~ 1 + (1|id), # FOR RAW FREA
      # formula= log((freq+0.00001)/startfreq) ~ 1 + (1|id), # FOR LOG RATIO
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

hist(myva$pseudoVa)
summary(myva$pseudoVa)

# vasur<-merge(myva,sur,by="site")
# vasur<-merge(myva,mersurclimate,by="site")
#
# vasurmean<-
#   vasur %>%
#   # group_by(site, pseudoVa) %>%
#   group_by(site, pseudoVa,bio1) %>%
#   summarise(
#             predmean=mean(sp_r),
#             survival1=sum(X1_survival>0),
#             survival2=sum(X2_survival>0)
#             )
#
#
# ggplot(vasurmean)+
#   geom_point(aes(x=pseudoVa,y=predmean,color=bio1),shape=16,size=3)+
#   stat_smooth(aes(x=pseudoVa,y=predmean), method="glm", formula=y~poly(x,2), color="grey",se=F)+
#   scale_color_gradientn("",colours = rev(redblue))+
#   #   xlab("Predictability year 1 (stabilizing r2 LOO)")+
#   # ylab("Survival at year 3 (survival/death)")+
#   theme_minimal()+
#   theme(legend.position = "none")
#
# ggplot(vasur)+
#   geom_point(aes(x=pseudoVa,y=sp_r))
#
# cor.test(vasur$sp_r, vasur$pseudoVa)
# cor.test(vasur$sp_r, vasur$pseudoVa, method="s")
#
# ggplot(vasur)+
#   geom_point(aes(x=pseudoVa,y=X2_survival))
# ggplot(vasur)+
#   geom_point(aes(x=pseudoVa,y=X2_survival))
#
# ggplot(vasurmean)+
#   geom_point(aes(x=pseudoVa,y=survival1))
# ggplot(vasurmean)+
#   geom_point(aes(x=pseudoVa,y=survival2))
