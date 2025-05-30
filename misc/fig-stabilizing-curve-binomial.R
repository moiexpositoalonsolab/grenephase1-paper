################################################################################
### Goal
### Estimate some basic metrics of ecotypes increasing in frequency across
### sites, means, variances, optimal places. Including poulation size

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

load(file=paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize.rda"))
load(file=paste0(myfolder,"/data-intermediate/ecotype_terminal_frequencies_long_raw_climate_popsize.rda"))


e<-enew <- ecotype_allyears_frequencies_long_raw_climate_popsize
# e<-enew <- ecotype_terminal_frequencies_long_raw_climate_popsize %>% 
  # rename(year=max_year, freq=maxfreq) # Use the terminal
enew[1:5,1:12]


################################################################################
# Make an easy model

# We want to do the stabilizing selection but with flower count

ggplot(enew)+
  geom_point(aes(y=freq*flowerscollected_corrected, x=sites_bio1-ecotypes_bio1))+
  geom_vline(xintercept = 0, lty='dotted')+
  # geom_vline(xintercept = mean(em$sites_latitude-em$ecotypes_latitude), lty='dashed')+
  # geom_vline(xintercept = weighted.mean(em$sites_latitude-em$ecotypes_latitude, w = em$freq), lty='solid')+
  xlab("Latitude transplant (site - ecotype)")+
  ylab("Change in ecotype freq")+
  theme_minimal()

ggplot(enew)+
  geom_point(aes(y=log(freq/startfreq), x=(sites_bio1-ecotypes_bio1)^2 ))+
  geom_vline(xintercept = 0, lty='dotted')+
  # geom_vline(xintercept = mean(em$sites_latitude-em$ecotypes_latitude), lty='dashed')+
  # geom_vline(xintercept = weighted.mean(em$sites_latitude-em$ecotypes_latitude, w = em$freq), lty='solid')+
  xlab("Latitude transplant (site - ecotype)")+
  ylab("Change in ecotype freq")+
  theme_minimal()


# Regression of log p1 p0
# get temp data
tmpdata<-
  enew %>% 
  dplyr::mutate(e2=(sites_bio1-ecotypes_bio1)^2,
                logp1p0=log(freq/startfreq),
                ecocounts=round(freq*flowerscollected_corrected)
                
  ) %>% 
  dplyr::select(ecocounts,logp1p0,flowerscollected_corrected,year,e2, sites_bio1,ecotypes_bio1) %>% 
  na.omit %>% 
  dplyr::filter(!is.infinite(logp1p0))
lm(data=tmpdata,
     logp1p0~e2
     ) ->
  lmmod
  
summary(lmmod)

# Predict the probabilities
tmpdata$predicted_logp1p0<-
  predict(lmmod,tmpdata,type="response")
summary(tmpdata$predicted_logp1p0)

ggplot(tmpdata)+
  geom_point(aes(y=exp(predicted_logp1p0),x=e2))


cor.test(tmpdata$logp1p0, tmpdata$predicted_logp1p0,method = 's')

# tmpdata$predicted_counts<-
#   tmpdata$predicted_probability* tmpdata$flowerscollected_corrected
# summary(tmpdata$predicted_counts)
# summary(tmpdata$predicted_probability)


################################################################################
# Binomial regression of log p1 p0
# get temp data
tmpdata<-
enew %>% 
  dplyr::mutate(e2=(sites_bio1-ecotypes_bio1)^2,
                logp1p0=log(freq/startfreq),
                ecocounts=round(freq*flowerscollected_corrected)
                
  ) %>% 
  dplyr::select(ecocounts,logp1p0,flowerscollected_corrected,year,e2, sites_bio1,ecotypes_bio1) %>% 
  na.omit
# fit glm
glm(data=tmpdata,
    cbind(ecocounts, flowerscollected_corrected - ecocounts) ~ year + e2,
    family = quasibinomial(link = "logit")) ->
  glmmod
summary(glmmod)

# Predict the probabilities
tmpdata$predicted_probability<-
  predict(glmmod,tmpdata,type="response")
tmpdata$predicted_counts<-
  tmpdata$predicted_probability* tmpdata$flowerscollected_corrected
summary(tmpdata$predicted_counts)
summary(tmpdata$predicted_probability)


# Visualize the decay in probability with distance
ggplot(tmpdata)+
  geom_point(aes(x=y=ecocounts,x=e2))
ggplot(tmpdata)+
  geom_point(aes(y=predicted_counts,y=e2))
ggplot(tmpdata)+
  geom_point(aes(y=predicted_probability,x=e2, color=year))

cor.test(tmpdata$ecocounts, tmpdata$predicted_counts,method = 's')
cor.test(tmpdata$logp1p0, tmpdata$predicted_probability,method = 's')

install.packages("pscl")
library(pscl)
pR2(glmmod)
install.packages("rcompanion")
library(rcompanion)
nagelkerke(glmmod)


################################################################################
#### Multinomial fitness #####

glmfunctionpbinomialP<-function(testfreq, testcount,testyear){
  lmmod<-glm(cbind(testfreq,testcount) ~ testyear + testbio,
             family = quasibinomial(link = logit))
  p=drop1(lmmod, test="Chisq")
  p=as.data.frame(p)[3,4]
  return(p)
}
glmfunctionpbinomialB<-function(testfreq, testcount,testyear){
  lmmod<-glm(cbind(testfreq,testcount) ~ testyear ,
             family = quasibinomial(link = logit))
  # p=drop1(lmmod, test="Chisq")
  # b=as.data.frame(p)[3,4]
  b=coef(lmmod)[2]
  return(b)
}

testrun<-
enew %>%
  dplyr::group_by(site,id) %>% 
  dplyr::summarise(b=glmfunctionpbinomialB(freq, 
                                           flowerscollected_corrected,
                                           as.numeric(year))
  )

load("grene/data/worldclim_ecotypesdata.rda")
load("grene/data/worldclim_sitesdata.rda")
colnames(worldclim_sitesdata)[-1] <- paste0("site_",colnames(worldclim_sitesdata)[-1])
colnames(worldclim_ecotypesdata)[-1] <- paste0("ecotype_",colnames(worldclim_ecotypesdata)[-1])
testrun<- 
  testrun %>% 
  merge(., 
        worldclim_sitesdata,
        by="site") %>% 
  merge(.,
        worldclim_ecotypesdata,
        by.x="id", by.y="ecotypeid"
        )

ggplot(testrun)+
  # geom_point(aes(y=b, x= ecotype_bio1))+
  stat_summary(aes(y=b, x= ecotype_bio1))
ggplot(testrun)+
  # geom_point(aes(y=b, x= ecotype_bio1))+
  stat_summary(aes(y=b, x= ecotype_bio18))


ggplot(testrun)+
  # geom_point(aes(y=b, x= ecotype_bio1))+
  stat_summary(aes(y=b, x= site_bio1, color=ecotype_bio1>15))


ggplot(testrun)+
  geom_jitter(aes(y=site_bio1, x= ecotype_bio1, color=b>0))+
  # scale_color_gradientn(colors = brewer.pal(9,"RdBu"))+
  theme_minimal()
ggplot(testrun)+
  geom_jitter(aes(y=site_bio18, x= ecotype_bio18, color=b>0))+
  # scale_color_gradientn(colors = brewer.pal(9,"RdBu"))+
  theme_minimal()



# testrun
# testfreq=enew$freq
# testcount=enew$flowerscollected_corrected
# testyear=as.numeric(enew$year)
# testbio=enew$sites_bio1

# 
# ecount<-
#   enew %>% 
#   data.frame() %>% 
#   dplyr::filter(site==4) %>% 
#   dplyr::mutate(ecocount=round(freq*flowerscollected_corrected)) %>% 
#   dplyr::select(site,rep,year,id, ecocount) %>% 
#   dplyr::mutate(id=paste0("X",id)) %>% 
#   pivot_wider(names_from = id, values_from = ecocount)
# ecount_long <- ecount %>%
#   pivot_longer(cols = starts_with("X"), 
#                names_to = "id", values_to = "ecocount")
# head(ecount)
# ecount %>% dim
# 
# 
# model <- glm(ecocount ~ year + id, data = ecount_long, family = quasibinomial(link = "logit"))
# model <- glm(ecocount ~ year + id, data = ecount_long, family = poisson(link = "log"))
# model <- glm(ecocount ~ year:id, data = ecount_long, family = binomial(link = "logit"))
# model %>% summary
# 
# library(MCMCglmm)
# modmc<-
#   MCMCglmm(
#     data=ecount, family = "zipoisson",
#     formula=cbind(X9764,X9723,X9579) ~ year,
#   )
# 
# # formula=cbind(X9764,X9723,X9579,X9596,X9518,X7217,X7164,X9910,X9526,X5811,X9784,X6987,X9629,X100002,X9992,X7346,X7103,X9920,X7125,X9736,X9758,X9653,X9595,X7255,X7008,X8387,X9506,X9879,X9537,X9574,X9534,X9619,X9837,X6979,X9743,X9803,X7378,X6929,X6898,X9941,X7077,X9542,X9726,X6150,X9562,X9602,X7143,X7025,X9933,X10006,X9940,X9814,X9597,X7347,X9741,X9699,X9637,X7126,X9925,X9555,X7031,X7316,X9977,X159,X6961,X9539,X7186,X9766,X8311,X7372,X10013,X6243,X9774,X9804,X7218,X8354,X8249,X9564,X9737,X9529,X6911,X9587,X9820,X9522,X766,X9594,X9657,X9528,X7282,X6195,X7323,X8230,X8247,X9557,X7288,X9947,X768,X9586,X9598,X6073,X7298,X9058,X9606,X9559,X7127,X9716,X9978,X6958,X772,X9507,X9527,X6188,X9371,X9481,X9394,X7002,X7067,X9640,X7333,X7165,X6945,X8312,X7003,X9769,X8376,X7028,X5151,X7036,X9544,X9510,X8231,X6939,X9759,X9416,X765,X9632,X9591,X9659,X9779,X7092,X7384,X9761,X9517,X9857,X9057,X8214,X7244,X8351,X9634,X6216,X9782,X9985,X6244,X9749,X6209,X7000,X265,X6177,X7063,X7521,X7071,X9521,X6184,X7013,X9649,X9781,X9813,X9549,X9584,X10011,X6938,X6074,X5165,X8240,X100001,X9547,X9612,X6013,X9713,X7106,X9548,X7411,X7209,X9524,X7273,X6180,X10002,X6963,X9470,X7394,X7353,X7062,X8357,X9323,X5768,X7268,X7203,X763,X5784,X9921,X9697,X9625,X7296,X9565,X7276,X9748,X9966,X9523,X9577,X9427,X9698,X9600,X6940,X9512,X9944,X6040,X9560,X9739,X10014,X6932,X9535,X9719,X9775,X6025,X7287,X6108,X6915,X9965,X5772,X9643,X7404) ~ year,
# library(nnet)
# 
# # Fit the multinom model
# model <- multinom(ecocount ~ year:id, data = ecount_long)
# 
# # View the model summary
# summary(model)
# 
# # library(nnet) # For multinomial regression
# # glm(data = ecount,
# #     c(X9506, X9649 ,X6074 ,X6958 ,X768) ~ year, 
# #     family = quasibinomial(link = logit))
# #     )
# 
# library(VGAM)
# model <- vglm(ecocount ~ year + id, family = multinomial, data = ecount_long)
# 
#   
# # ecount %>% head %>% 
# #   pivot_wider(names_from = id, values_from = ecocount)