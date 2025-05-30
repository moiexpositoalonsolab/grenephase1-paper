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

################################################################################
#### Location
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
# myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

p0
load("data-intermediate/freq")

########Ç################################################################
## Check that the change in frequency is inverse to number of flowers

# ## How does the change in frequency
# # Do the plot
# frequency_and_flowers<-
# ggplot(e)+
#   geom_point(aes(y=freq-startfreq , x=flowerscollected_corrected, color=as.numeric(year)), alpha=0.5)+
#   scale_color_continuous("year",type = 'viridis')+
#   ylab("Frequency change (p1-p0)")+
#   xlab("Number of flowers in pool")+
#   theme_minimal()
# 
# frequency_and_flowes_log<-
# ggplot(e)+
#   geom_point(aes(y=freq-startfreq , x=flowerscollected_corrected, color=as.numeric(year)), alpha=0.5)+
#   scale_y_log10()+
#   scale_x_log10()+
#   stat_smooth(aes(y=freq-startfreq , x=flowerscollected_corrected, color=as.numeric(year)), method='glm', color="grey")+
#   scale_color_continuous("year",type = 'viridis')+
#   ylab("Frequency change (p1-p0)")+
#   xlab("Number of flowers in pool")+
#   theme_minimal()
# 
# freq_plot_grid<-
#   plot_grid(frequency_and_flowers,frequency_and_flowers_log,ncol=2)
# freq_plot_grid
# 
# save_plot("figs/fig-allele-frequency-inverse-relationship-flower-count.pdf",freq_plot_grid ,base_height = 6,base_width = 12)
# save_plot("figs/fig-allele-frequency-inverse-relationship-flower-count.png",freq_plot_grid ,base_height = 6,base_width = 12)
# 
# 
# ################################################################################
# ## Variance in frequency
# evariance<-
#   group_by(site, replicate)+
#   
#   
#   
# ggplot(e)+
#   geom_point(aes(y=var(freq-startfreq) , x=flowerscollected_corrected, color=as.numeric(year)), alpha=0.5)+
#   scale_color_continuous("year",type = 'viridis')+
#   ylab("Frequency change (p1-p0)")+
#   xlab("Number of flowers in pool")+
#   theme_minimal()
# 
# 
# 
# 
# # moiR::addggregression
# # lm_eq<-function (y, x, tex = TRUE) 
# # {
# #   mylm <- lm(y ~ x)
# #   p = format(coefficients(summary(mylm))[2, 4], digits = 3, 
# #              scientific = TRUE)
# #   r2 = round(summary(mylm)$r.squared, digits = 3)
# #   b = round(coefficients(summary(mylm))[2, 1], digits = 3)
# #   if (tex == FALSE) {
# #     return(sprintf("R2= %s, b= %s, p= %s", r2, b, p))
# #   }
# #   else {
# #     return(paste0("$R^2 = $", r2, ", $\\beta = $", b, ", $ p = $", 
# #                   p))
# #   }
# # }
# # lm_eq(y=log10(e$freq-e$startfreq), x=log10(e$flowerscollected_corrected))
# # 


# ################################################################################
# ### Get some variance
# 
# avar<-
#   a %>% 
#   dplyr::filter(year==1) %>% 
#   group_by(site, rep, longitude,latitude, sites_bio1, flowers) %>% 
#   dplyr::rename(p1=freq, p0=startfreq) %>% 
#   summarise(varp= var( (p1-p0)/ sqrt(p0*(1-p0) ) ) ) 
# 
# 
# ggplot(avar)+
#   geom_point(aes(y=varp, x=1/flowers))+
#   scale_x_log10()+scale_y_log10()+
#   geom_abline(slope = 1)+
#   theme_minimal()
# 
#   
#   
#   avar<-
#     a %>% 
#     dplyr::filter(year==1) %>% 
#     group_by(site, rep, longitude,latitude, sites_bio1, flowers) %>% 
#     dplyr::rename(p1=freq, p0=startfreq) %>% 
#     summarise(varp= var( (p1-p0)/ (p0*(1-p0) ) ) ) 
#   
#   
#   ggplot(avar)+
#     geom_point(aes(y=varp, x=flowers ))+
#     scale_x_log10()+scale_y_log10()+
#     # geom_abline(slope = 1)+
#     theme_minimal()
  



#####*********************************************************************######
#####* Alleles 

load(file="data-intermediate/allele_allyears_frequencies_long_raw_climate_flowers.rda")
a<-allele_allyears_frequencies_long_raw_climate_flowers

a[1:5,1:10]



## How does the change in frequency
a %>% 
  sample_n(10000) %>% 
  ggplot(.)+
  geom_point(aes(y=(freq-startfreq)^2 , x=flowers), color="white",size=2)+
  geom_point(aes(y=(freq-startfreq)^2 , x=flowers), color="black",size=1.5)+
  geom_point(aes(y=(freq-startfreq)^2 , x=flowers, color=as.numeric(year)))+
  # geom_point(
  #   data=expectation,
  #   aes(y=(freq-startfreq)^2 , x=freq, color=as.numeric(year)))+
  # 
  scale_y_log10()+
  scale_x_log10()+
  scale_color_gradientn("year",colors = brewer.pal(8,"Greens")[-1])+
  ylab("Frequency change (p1-p0)^2")+
  xlab("Number of flowers in pool")+
  theme_minimal()

## How frequency change depends on starting
a %>% 
  dplyr::filter(flowers==100) %>% 
  sample_n(1000) %>% 
  # dplyr::mutate(freq=ifelse(freq>0.5,1-freq,freq)) %>% 
  ggplot(.)+
  geom_point(aes(y=(freq-startfreq)^2 , x=freq), color="white",size=2)+
  geom_point(aes(y=(freq-startfreq)^2 , x=freq), color="black",size=1.5)+
  geom_point(aes(y=(freq-startfreq)^2 , x=freq, color=flowers))+
  # geom_point(aes(y=(1/flowers)^2, x=freq, color=flowers), color="green")+
  geom_point(aes(y=(0.5-freq)^2, x=freq, color=flowers), color="green")+
  # geom_point(aes(y=ceiling(100*freq) , x=freq))+
  # geom_line(data=expectation,aes(y=Varp10 , x=p))+
  # geom_line(data=expectation,aes(y=Varp100 , x=p))+
  # geom_line(data=expectation,aes(y=Varp1000 , x=p))+
  # geom_line(data=expectation,aes(y=ceiling(freq*100) , x=p))+
  # scale_y_log10()+
  # scale_x_log10()+
  scale_color_gradientn("Year",colors = brewer.pal(9,"Greys")[-1])+
  ylab("Frequency change (p1-p0)^2")+
  xlab("Starting frequency (p0)")+
  # facet_wrap(~year)+
  theme_minimal()



#####*********************************************************************######
##### *Generation 1 ####
######
# Get the observed
test<-
  a %>% 
  dplyr::filter(year==1) %>% 
  # sample_n(10000) %>% 
  # dplyr::filter(flowers==100) %>% 
  dplyr::filter(flowers>3) %>% 
  # dplyr::mutate(freq=ifelse(freq>0.5,1-freq,freq)) %>% 
  dplyr::mutate(p0=round(startfreq*100)/100)  %>% 
  dplyr::mutate(N=ceiling(flowers/10)*10)  %>% 
  dplyr::mutate(N=ifelse(N>100,100, N))  %>% 
  dplyr::filter(N>1) %>% 
  dplyr::group_by(p0,N) %>% 
  dplyr::summarize(meandeltap=mean((freq-startfreq)^2)) 
head(test)
dim(test)


# colorful
p10exp<-
  ggplot(test)+
  geom_point(aes(y=meandeltap,x=p0, color=factor(N)))+
  # geom_line(data=expectation,aes(y=Varp/200 , x=p0))+
  # geom_line(data=expectation,aes(y=Varp/5 , x=p0))+
  geom_line(aes(y=(p0*(1-p0)/N),x=p0, group=N, color=factor(N)))+
  ylab(TeX("Mean squared change $(p_1-p_0)^2$"))+
  xlab(TeX("Starting frequency $(p_0)$"))+
  theme_minimal()
p10exp  

# gradient
p10exp<-
  ggplot(test)+
  geom_point(aes(y=meandeltap,x=p0, color=N))+
  # geom_line(data=expectation,aes(y=Varp/200 , x=p0))+
  # geom_line(data=expectation,aes(y=Varp/5 , x=p0))+
  geom_line(aes(y=(p0*(1-p0)/N),x=p0, group=N, color=N))+
  ylab(TeX("Mean squared change $(p_1-p_0)^2$"))+
  xlab(TeX("Starting frequency $(p_0)$"))+
  theme_minimal()
p10exp  

save_plot("figs/fig-allele-frequency-change-expection.pdf",
          p10exp,base_height = 5,base_width = 6)
save_plot("figs/fig-allele-frequency-change-expection.png",
          p10exp,base_height = 5,base_width = 6)



ggplot(test)+
  geom_point(aes(y=meandeltap,x=p0, color=N))+
  # geom_line(data=expectation,aes(y=Varp/200 , x=p0))+
  # geom_line(data=expectation,aes(y=Varp/5 , x=p0))+
  geom_line(aes(y=(p0*(1-p0)/N),x=p0, group=N, color=N))+
  ylab(TeX("Mean squared change $(p_1-p_0)^2$"))+
  xlab(TeX("Starting frequency $(p_0)$"))+
  scale_y_log10()+
  theme_minimal()

test %>% 
  dplyr::filter(N>20) %>% 
ggplot(.)+
  geom_point(aes(y=meandeltap/(p0*(1-p0)/N),x=p0, color=N))+
  # geom_line(data=expectation,aes(y=Varp/200 , x=p0))+
  # geom_line(data=expectation,aes(y=Varp/5 , x=p0))+
  # geom_line(aes(y=(p0*(1-p0)/N),x=p0, group=N, color=N))+
  geom_hline(yintercept = 1, color="darkgrey")+
  ylab(TeX("Variance excess over neutrality $(p_1-p_0)^2$  / $p_0(1-p_0)/N$"))+
  xlab(TeX("Starting frequency $(p_0)$"))+
  theme_minimal()->
p10ratio
p10ratio

save_plot("figs/fig-allele-frequency-change-expection-ratio-over-neutrality.pdf",
          p10ratio,base_height = 5,base_width = 6)
save_plot("figs/fig-allele-frequency-change-expection-ratio-over-neutrality.png",
          p10ratio,base_height = 5,base_width = 6)


sink("figs/test-allele-frequency-change-expection-ratio-over-neutrality.txt")
attach(test)
hist(meandeltap/(p0*(1-p0)/N))
summary(meandeltap-(p0*(1-p0)/N))
summary(meandeltap/(p0*(1-p0)/N))
t.test(meandeltap/(p0*(1-p0)/N))
sink()

#####*********************************************************************######
#####* Generation 2 ####


test2<-
  a %>% 
  dplyr::filter(year==2) %>% 
  dplyr::filter(flowers>3) %>% 
  # dplyr::mutate(freq=ifelse(freq>0.5,1-freq,freq)) %>% 
  dplyr::mutate(p0=round(startfreq*100)/100)  %>% 
  dplyr::mutate(N=ceiling(flowers/10)*10)  %>% 
  dplyr::mutate(N=ifelse(N>100,100, N))  %>% 
  dplyr::filter(N>1) %>% 
  dplyr::group_by(p0,N) %>% 
  dplyr::summarize(meandeltap=mean((freq-startfreq)^2)) 
head(test)
dim(test)

ggplot(test2)+
  geom_point(aes(y=meandeltap,x=p0, color=factor(N)))+
  # geom_line(data=expectation,aes(y=Varp/200 , x=p0))+
  # geom_line(data=expectation,aes(y=Varp/5 , x=p0))+
  geom_line(aes(y=2*(p0*(1-p0)/N),x=p0, group=N, color=factor(N)))+
  ylab(TeX("Mean squared change $(p_1-p_0)^2$ vs   $t \\times p_0(1-p_0)/N$"))+
  xlab(TeX("Starting frequency $(p_0)$"))+
  labs(color="N")+
  theme_minimal()

test2 %>% 
  dplyr::filter(meandeltap / (3*(p0*(1-p0)/N)) < 5) %>%
  ggplot(.)+
  geom_point(aes(y=meandeltap / (3*(p0*(1-p0)/N)),x=p0, color=factor(N)))+
  # geom_line(data=expectation,aes(y=Varp/200 , x=p0))+
  # geom_line(data=expectation,aes(y=Varp/5 , x=p0))+
  # geom_line(aes(y=3*(p0*(1-p0)/N),x=p0, group=N, color=factor(N)))+
  ylab(TeX("Mean squared change $(p_1-p_0)^2$ /  $t \\times p_0(1-p_0)/N$"))+
  xlab(TeX("Starting frequency $(p_0)$"))+
  geom_hline(yintercept = 1, color="darkgrey")+
  labs(color="N")+
  theme_minimal()


attach(test2)
hist(meandeltap/(2*p0*(1-p0)/N))
summary(meandeltap/(2*p0*(1-p0)/N))


attach(test3)
hist(meandeltap/(3*p0*(1-p0)/N))
summary(meandeltap/(3*p0*(1-p0)/N))


#####*********************************************************************######
#####* Generation 3 ####

test3<-
  a %>% 
  dplyr::filter(year==3) %>% 
  dplyr::filter(flowers>3) %>% 
  # dplyr::mutate(freq=ifelse(freq>0.5,1-freq,freq)) %>% 
  dplyr::mutate(p0=round(startfreq*100)/100)  %>% 
  dplyr::mutate(N=ceiling(flowers/10)*10)  %>% 
  dplyr::mutate(N=ifelse(N>100,100, N))  %>% 
  dplyr::filter(N>1) %>% 
  dplyr::group_by(p0,N) %>% 
  dplyr::summarize(meandeltap=mean((freq-startfreq)^2)) 
head(test)
dim(test)


# ggplot(test)+
#   geom_point(aes(y=meandeltap,x=p0, color=factor(N)))+
#   # geom_line(data=expectation,aes(y=Varp/200 , x=p0))+
#   # geom_line(data=expectation,aes(y=Varp/5 , x=p0))+
#   geom_line(aes(y=1*(p0*(1-p0)/(2*N)),x=p0, group=N, color=factor(N)))+
#   ylab(TeX("Mean squared change $(p_1-p_0)^2$ vs   $p_0(1-p_0)/N$"))+
#   xlab(TeX("Starting frequency $(p_0)$"))+
#   labs(color="N")+
#   theme_minimal()

ggplot(test)+
  geom_point(aes(y=meandeltap,x=p0, color=factor(N)))+
  # geom_line(data=expectation,aes(y=Varp/200 , x=p0))+
  # geom_line(data=expectation,aes(y=Varp/5 , x=p0))+
  geom_line(aes(y=1*(p0*(1-p0)/(N)),x=p0, group=N, color=factor(N)))+
  ylab(TeX("Mean squared change $(p_1-p_0)^2$ vs   $p_0(1-p_0)/N$"))+
  xlab(TeX("Starting frequency $(p_0)$"))+
  labs(color="N")+
  theme_minimal()

p3<-
ggplot(test3)+
  geom_point(aes(y=meandeltap,x=p0, color=factor(N)))+
  # geom_line(data=expectation,aes(y=Varp/200 , x=p0))+
  # geom_line(data=expectation,aes(y=Varp/5 , x=p0))+
  geom_line(aes(y=3*(p0*(1-p0)/N),x=p0, group=N, color=factor(N)))+
  ylab(TeX("Mean squared change $(p_1-p_0)^2$ vs   $t \\times p_0(1-p_0)/N$"))+
  xlab(TeX("Starting frequency $(p_0)$"))+
  labs(color="N")+
  theme_minimal()
p3

p3over<-
  test3 %>% 
  dplyr::filter(meandeltap / (3*(p0*(1-p0)/N)) < 5) %>% 
ggplot(.)+
  geom_point(aes(y=meandeltap / (3*(p0*(1-p0)/N)),x=p0, color=factor(N)))+
  # geom_line(data=expectation,aes(y=Varp/200 , x=p0))+
  # geom_line(data=expectation,aes(y=Varp/5 , x=p0))+
  # geom_line(aes(y=3*(p0*(1-p0)/N),x=p0, group=N, color=factor(N)))+
  ylab(TeX("Mean squared change $(p_1-p_0)^2$ /  $t \\times p_0(1-p_0)/N$"))+
  xlab(TeX("Starting frequency $(p_0)$"))+
  geom_hline(yintercept = 1, color="darkgrey")+
  labs(color="N")+
  theme_minimal()
p3over

save_plot("figs/fig-allele-frequency-change-expection-gen3.pdf",
          p3,base_height = 5,base_width = 6)
save_plot("figs/fig-allele-frequency-change-expection-gen3.png",
          p3,base_height = 5,base_width = 6)

save_plot("figs/fig-allele-frequency-change-expection-ratio-over-neutrality-gen3.pdf",
          p3over,base_height = 5,base_width = 6)
save_plot("figs/fig-allele-frequency-change-expection-ratio-over-neutrality-gen3.png",
          p3over,base_height = 5,base_width = 6)

