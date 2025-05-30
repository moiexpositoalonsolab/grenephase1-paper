################################################################################
### Goal
### Get a similar plot as the Diffusion map for those non-believers on
### things that are not a PCA for decomposition of variable
################################################################################
# 
# #### If in google drive
# from google.colab import drive
# drive.mount('/content/drive')
# %load_ext rpy2.ipython


################################################################################
# %%R

# detach(ggplot2)
# detach("package:ggplot2", unload=TRUE)

# library(grene)
library(tidyverse)
library(dplyr)

library(ggplot2)
library(RColorBrewer)
library(tidyr)
# theme_set(theme_cowplot())
# theme_set(theme_classic)

################################################################################
# %%R
# Set location

# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE

setwd(myfolder)

################################################################################
# Load files

file_plot_pca12="/fig-pca12-bio1-ecotypes.pdf"

# Load ecotypes and their ID
eco0<-read.table("data/founder_ecotype_frequency.txt", header=F) %>% 
  dplyr::rename(id=V1,e0=V2)
ecofreq<-read.table("data/merged_ecotype_frequency.txt", header=T) 
ecodelta<-read.table("data/delta_ecotype_freq.txt", header=T)


df_wide<-ecodelta
df_wide<-ecofreq


################################################################################
#### PCA analysis
# run pca
pcamod<-prcomp(df_wide[,-c(1)])
pcamod<-prcomp(df_wide)
# plot(pcamod)
pdat<-pcamod$x
dim(pdat)

pdat<-data.frame(pdat)

# calculate variance explained
percentageesplained<- pcamod$sdev^2 / sum(pcamod$sdev^2)
percentageesplained<-round(percentageesplained*100,1)
head(percentageesplained)
percentageesplainedlabels<- paste0("PC",1:length(percentageesplained)," (", percentageesplained,"%)")

# Quick plot
normalize<- function(x) {(x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm=T))}


ggplot(pdat) +
  geom_point(aes(x=scale(PC1)+0.1, y=scale(PC2)+0.01))+
  scale_y_log10()+
  scale_x_log10()+
  theme_minimal()

################################################################################
### Add the climate to the 
##############
# Load climates
load("grene/data/worldclim_sitesdata.rda")
load("grene/data/locations_data.rda")
# worldclim_sitesdata<-apply(worldclim_sitesdata,2,as.numeric)
# worldclim_sitesdata$Latitude<-sites.clim$LATITUDE

worldclim_sitesdata<-worldclim_sitesdata %>% data.frame %>% mutate_all(as.numeric)
locations_data<-locations_data %>% data.frame %>%  mutate_all(as.numeric)

##############
# prepare for merge
pdat<-data.frame(pdat)
df_wide<-data.frame(df_wide)
pdat$site_rep<-df_wide$site_rep
pdat$site<-unlist(lapply(pdat$site_rep, function(i) strsplit(i,split="_")[[1]][1]))
tail(pdat)

ggplot(pdat) +
  geom_point(aes(x=PC1, y=PC2, color=site))

ggplot(pdat) +
  geom_point(aes(x=PC1, y=PC3, color=site))

# Merge
pdat %>% head
pmer<-merge(pdat,by.x="site", worldclim_sitesdata, by.y="site")

# Plot
pcaplot<-
ggplot(pmer) +
  geom_point(aes(x=PC1, y=PC2, color=bio1))+
  scale_color_gradientn("Temp. (C)",colours =  brewer.pal(9,"Reds"))+
  xlab(percentageesplainedlabels[1])+
  ylab(percentageesplainedlabels[2])+
  theme_minimal()
save_plot(pcaplot,
          filename=paste0(myfolder,file_plot_pca12),base_height = 4.5,base_width = 5)
