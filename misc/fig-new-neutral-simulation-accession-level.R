################################################################################
### Goal
### Create a realistic accession sorting simulation and compare Var(delta p)
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
# Load the genome matrix of LD prunned SNPs
g<-read.table("data/genotype_reference_dosage.txt",header = F)
head(g)
# X=g[,1:5] %>% as.matrix()
p0<-apply(g,1,mean)*0.5
hist(p0)
# p<-apply(g,1,function(x)sum(x==TRUE))/(2*231)
# table(p)
# (apply(X,1,sum)) %>% hist
# p=(apply(X,1,sum)/(2*231))
# hist(p)

# Samples of flowers
load("grene/data/samples_data.rda")
head(samples_data)

flowers<-
  samples_data %>%
  dplyr::filter(year %in% c("2018","2019","2020")) %>%
  dplyr::mutate(year=year-2017) %>%
  dplyr::group_by(site,plot,year) %>%
  dplyr::summarise(flowersum=sum(flowerscollected))

flowers %>% head

# check how many data points
flowers %>%
  dplyr::filter(year==1, flowersum>3) %>% nrow


#####*********************************************************************######
#####* Load frequencies  ######

load(file="data-intermediate/allele_allyears_frequencies_long_raw_climate_flowers.rda")
a<-allele_allyears_frequencies_long_raw_climate_flowers
a[1:5,1:10]

#####*********************************************************************######
#####* Overall variance  ######

obsmean<-
  a %>%
  dplyr::filter(year==1) %>%
  dplyr::filter(flowers>3) %>%
  # dplyr::mutate(freq=ifelse(freq>0.5,1-freq,freq)) %>%
  dplyr::mutate(p0=round(startfreq*100)/100)  %>%
  dplyr::mutate(N=ceiling(flowers/10)*10)  %>%
  dplyr::mutate(N=ifelse(N>100,100, N))  %>%
  dplyr::filter(N>1) %>%
  # dplyr::group_by(site,rep) %>% # Group by site
  dplyr::group_by(p0,N) %>% # Group by starting frequency and sample size
  dplyr::summarize(
    wfdeltap=mean((startfreq*(1-startfreq))/(2*N)),
    meandeltap=mean((freq-startfreq)^2)
  )

obsmean %>% head

sink("tables/ratio_obs_deltap_wfneutral.txt")
summary(obsmean$meandeltap/obsmean$wfdeltap)
hist(obsmean$meandeltap/obsmean$wfdeltap)
wilcox.test(obsmean$meandeltap/obsmean$wfdeltap)
t.test(obsmean$meandeltap/obsmean$wfdeltap)
sink()

#####*********************************************************************######
#####* Observed Variance (delta P)  ######


obs<-
  a %>%
  dplyr::filter(year==1) %>%
  dplyr::filter(flowers>3) %>%
  # dplyr::mutate(freq=ifelse(freq>0.5,1-freq,freq)) %>%
  dplyr::mutate(p0=round(startfreq*100)/100)  %>%
  dplyr::mutate(N=ceiling(flowers/10)*10)  %>%
  dplyr::mutate(N=ifelse(N>100,100, N))  %>%
  dplyr::filter(N>1) %>%
  dplyr::group_by(p0,N) %>% # Group by starting frequency and sample size
  dplyr::summarize(meandeltap=mean((freq-startfreq)^2))
head(obs)
dim(obs)
ggplot(obs)+geom_histogram(aes(x=meandeltap))

#####* Create neutarl expectations Var (delta P) ######

### Ecotype filtering uniform

x=1
neutraluniform<-
  lapply(1:length(allflowers),
  # lapply(1:10,#DEBUG
         FUN= function(x){
           afake<-flowers2plus[x,]
           acc<-sample(1:231,afake$flowersum, replace = T)
           p<-apply(g[,acc],1,mean)*0.5

           afake <- afake %>%
             slice(rep(1, length(p))) %>%         # replicate row by length(p)
             mutate(freq = p, startfreq = p0)     # assign columns

           return(afake)
         }
  ) %>%
  do.call(rbind,.)

neuuni<-
  # pipe the sampled dataset
  neutraluniform %>%
  rename(flowers=flowersum) %>%
  # repeat the same as observation
  dplyr::filter(year==1) %>%
  dplyr::filter(flowers>3) %>%
  dplyr::mutate(p0=round(startfreq*100)/100)  %>%
  dplyr::mutate(N=ceiling(flowers/10)*10)  %>%
  dplyr::mutate(N=ifelse(N>100,100, N))  %>%
  dplyr::filter(N>1) %>%
  dplyr::group_by(p0,N) %>%
  dplyr::summarize(meandeltap=mean((freq-startfreq)^2))
head(neuuni)


### Ecotype filtering poisson
lambdapoisson=7.2
x=1
neutralpois<-
  lapply(1:length(allflowers),
         # lapply(1:10,#DEBUG
         FUN= function(x){
           # Create dataset of all data points
           afake<-flowers2plus[x,]
           # Sample accessions in every sample
           # acc<-sample(1:231,afake$flowersum, replace = T)
           # p<-apply(g[,acc],1,mean)*0.5
           accfreq<-rpois(231,lambda = lambdapoisson)
           acc<-sample(1:231,allflowers[x], replace = T,prob = accfreq)
           p<-apply(g[,acc],1,mean)*0.5
           # Fill the dataset
           afake <- afake %>%
             slice(rep(1, length(p))) %>%         # replicate row by length(p)
             mutate(freq = p, startfreq = p0)     # assign columns

           return(afake)
         }
  ) %>%
  do.call(rbind,.)

neupois<-
  # pipe the sampled dataset
  neutralpois %>%
  rename(flowers=flowersum) %>%
  # repeat the same as observation
  dplyr::filter(year==1) %>%
  dplyr::filter(flowers>3) %>%
  dplyr::mutate(p0=round(startfreq*100)/100)  %>%
  dplyr::mutate(N=ceiling(flowers/10)*10)  %>%
  dplyr::mutate(N=ifelse(N>100,100, N))  %>%
  dplyr::filter(N>1) %>%
  dplyr::group_by(p0,N) %>%
  dplyr::summarize(meandeltap=mean((freq-startfreq)^2))
head(neupois)


### PLOTS ###


# Observation vs theory
ggplot(obs)+
  geom_point(aes(y=meandeltap,x=p0, color=factor(N)))+
  # geom_line(data=expectation,aes(y=Varp/200 , x=p0))+
  # geom_line(data=expectation,aes(y=Varp/5 , x=p0))+
  geom_line(aes(y=1*(p0*(1-p0)/(N)),x=p0, group=N, color=factor(N)))+
  ylab(TeX("Mean squared change $(p_1-p_0)^2$ vs   $p_0(1-p_0)/N$"))+
  xlab(TeX("Starting frequency $(p_0)$"))+
  labs(color="N")+
  theme_minimal()

# Uniform vs theory

ggplot(neuuni)+
  geom_point(aes(y=meandeltap,x=p0, color=factor(N)))+
  # geom_line(data=expectation,aes(y=Varp/200 , x=p0))+
  # geom_line(data=expectation,aes(y=Varp/5 , x=p0))+
  geom_line(aes(y=1*(p0*(1-p0)/(N)),x=p0, group=N, color=factor(N)))+
  ylab(TeX("Mean squared change $(p_1-p_0)^2$ vs   $p_0(1-p_0)/N$"))+
  xlab(TeX("Starting frequency $(p_0)$"))+
  labs(color="N")+
  theme_minimal()

# Poisson vs theory


ggplot(neupois)+
  geom_point(aes(y=meandeltap,x=p0, color=factor(N)))+
  # geom_line(data=expectation,aes(y=Varp/200 , x=p0))+
  # geom_line(data=expectation,aes(y=Varp/5 , x=p0))+
  geom_line(aes(y=1*(p0*(1-p0)/(N)),x=p0, group=N, color=factor(N)))+
  ylab(TeX("Mean squared change $(p_1-p_0)^2$ vs   $p_0(1-p_0)/N$"))+
  xlab(TeX("Starting frequency $(p_0)$"))+
  labs(color="N")+
  theme_minimal()

# All together
library(ggnewscale)

ggplot()+
  geom_point(data=neupois,aes(y=meandeltap/1*(p0*(1-p0)/(N)),x=p0), color="red")+
  # geom_point(data=neuuni,aes(y=meandeltap/1*(p0*(1-p0)/(N)),x=p0), color="blue")+
  geom_point(data=obs,aes(y=meandeltap/1*(p0*(1-p0)/(N)),x=p0), color="black")+
  # geom_line(data=expectation,aes(y=Varp/200 , x=p0))+
  # geom_line(data=expectation,aes(y=Varp/5 , x=p0))+
  # geom_line(data=neupois,aes(y=1*(p0*(1-p0)/(N)),x=p0, group=N, color=(N)))+
  ylab(TeX("Mean squared change $(p_1-p_0)^2$ vs   $p_0(1-p_0)/N$"))+
  xlab(TeX("Starting frequency $(p_0)$"))+
  scale_y_log10()+
  labs(color="N")+
  theme_minimal()

ggplot()+
  geom_point(data=neupois,aes(y=meandeltap/1*(p0*(1-p0)/(N)),x=p0,color=N))+
  scale_color_gradientn("Poisson",colours = brewer.pal(5,"Reds"))+
  new_scale_color()+
  geom_point(data=neuuni,aes(y=meandeltap/1*(p0*(1-p0)/(N)),x=p0,color=N))+
  scale_color_gradientn("Uniform",colours = brewer.pal(5,"Blues"))+
  new_scale_color()+
  geom_point(data=obs,aes(y=meandeltap/1*(p0*(1-p0)/(N)),x=p0,color=N))+
  scale_color_gradientn("Observed",colours = brewer.pal(5,"Greys"))+
  # geom_line(data=neupois,aes(y=1*(p0*(1-p0)/(N)),x=p0, group=N, color=(N)))+
  ylab(TeX("Mean squared change $(p_1-p_0)^2$ vs   $p_0(1-p_0)/N$"))+
  xlab(TeX("Starting frequency $(p_0)$"))+
  scale_y_log10()+
  labs(color="N")+
  theme_minimal()

ggplot()+
  # Filtering poisson
  geom_point(data=neupois,aes(y=meandeltap,x=p0,color=N))+
  scale_color_gradientn("Poisson",colours = brewer.pal(5,"Greens")[-1])+
  # Filtering uniform
  new_scale_color()+
  geom_point(data=neuuni,aes(y=meandeltap,x=p0,color=N))+
  scale_color_gradientn("Uniform",colours = brewer.pal(5,"Blues")[-1])+
  # Wright fisher
  new_scale_color()+
  geom_line(data=neupois,aes(y=1*(p0*(1-p0)/(N)),x=p0, group=N, color=N))+
  scale_color_gradientn("Uniform",colours = brewer.pal(5,"Greys")[-1])+
  # Observed
  new_scale_color()+
  geom_point(data=obs,aes(y=meandeltap,x=p0,color=N))+
  scale_color_gradientn("Observed",colours = brewer.pal(9,"Reds")[-1])+
  # Labs
  ylab(TeX("Mean squared change $(p_1-p_0)^2$ vs   $p_0(1-p_0)/N$"))+
  xlab(TeX("Starting frequency $(p_0)$"))+
  scale_y_log10()+
  labs(color="N")+
  theme_minimal()

#WF diploid
finalplot<-
ggplot()+
  # Filtering poisson
  geom_jitter(data=neupois,aes(y=meandeltap/(p0*(1-p0)/(2*N)),x=factor(N),color=N),height = 0)+
  geom_boxplot(data=neupois,aes(y=meandeltap/(p0*(1-p0)/(2*N)),x=factor(N),color=N),fill="white",alpha=0.5)+
  scale_color_gradientn("Simulation\nsorting\nPoisson",colours = brewer.pal(5,"Greens")[-1])+
  # Filtering uniform
  new_scale_color()+
  geom_jitter(data=neupois,aes(y=meandeltap/(p0*(1-p0)/(2*N)),x=factor(N),color=N),height = 0)+
  geom_boxplot(data=neupois,aes(y=meandeltap/(p0*(1-p0)/(2*N)),x=factor(N),color=N),fill="white",alpha=0.5)+
  scale_color_gradientn("Simulation\nsorting\nuniform",colours = brewer.pal(5,"Blues")[-1])+
  # Wright fisher
  # new_scale_color()+
  # geom_line(data=neupois,aes(y=1*(p0*(1-p0)/(N)),x=p0, group=N, color=N))+
  # scale_color_gradientn("Uniform",colours = brewer.pal(5,"Greys")[-1])+
  # Observed
  new_scale_color()+
  geom_jitter(data=obs,aes(y=meandeltap/(p0*(1-p0)/(2*N)),x=factor(N),color=N),height = 0)+
  geom_boxplot(data=obs,aes(y=meandeltap/(p0*(1-p0)/(2*N)),x=factor(N),color=N),fill="white",alpha=0.5)+
  scale_color_gradientn("Observed",colours = brewer.pal(9,"Reds")[-1])+
  # Line
  geom_hline(yintercept = 1,color="grey")+
  annotate(geom = "text", x=factor(90),y=0.9,label="Neutral WF expectation", color="darkgrey")+
  # Labs
  ylab(TeX("Ratio mean squared change $(p_1-p_0)^2$ /  $p_0(1-p_0)/2N$"))+
  xlab(TeX("Sample size $(N)$"))+
  # scale_y_log10()+
  labs(color="N")+
  theme_minimal()


finalplot

save_plot("figs/fig-allele-frequency-change-expection-ratio-over-ecotype-sorting-simulation-and-theoreticalWFneutrality-gen1.pdf",
          finalplot,base_height = 6,base_width = 8)
save_plot("figs/fig-allele-frequency-change-expection-ratio-over-ecotype-sorting-simulation-and-theoreticalWFneutrality-gen1.png",
          finalplot,base_height = 6,base_width = 8)

