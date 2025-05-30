rm(list=ls())

## This script plots the Vs Wmax and stabilizing selection curves
library(dplyr)
library(ggplot2)
library(RColorBrewer)


prefix <- "LOCAL_ADAPTATION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

## draw the stabilizing selection curves for each ecotype

ecotype_climate <- read.delim("data-external/bioclimvars_ecotypes_era5.csv",sep=",")
site_climate <- read.delim("data-external/bioclimvars_sites_era5_year_2018.csv",sep=",")
meta <- read.delim("data/merged_sample_table.csv",sep=",")

## load Vs Wmax 
df<-read.table("data-intermediate/LOCAL_ADAPTATION_ecotype_Wmax_Vs.txt",header = T)

## load ecotype_info
ecotype_climate <- read.delim("data-external/bioclimvars_ecotypes_era5.csv",sep=",")
colnames(ecotype_climate) <- c("ecotype",paste0("bio",1:19,"_ecotype"))

ecotype_info <- read.delim("data/1001g_regmap_grenet_ecotype_info_corrected_2024May16.csv",sep = ",")
colnames(ecotype_info)[c(1,5,6)] <- c("ecotype","latitude_ecotype","longitude_ecotype")

ecotype_info <- ecotype_info %>%
  filter(GrENE.net == TRUE) %>%
  arrange(match(ecotype,Vsresults$ecotype)) %>%
  left_join(.,Vsresults,by="ecotype") %>%
  left_join(.,ecotype_climate,by="ecotype")

site_climate <- site_climate %>%
  filter(site %in% meta$site) %>%
  arrange(match(site,unique(meta$site)))

Vs_0.9 <- quantile(df$Vs,0.9)
wmax_0.9 <- quantile(df$Wmax,0.9)

stabilizing_selection <- function(x,z0,wmax,vs){ wmax * exp( - (x-z0)^2 / vs)} 

pdf("stabilizing_curve_bio1.pdf",width = 7,height = 5)
plot(0,0,xlim=c(0,25),ylim=c(0,2.5),type="n",bty = "n",yaxt="n",ylab=NA,xlab=NA)
for(i in 1:231){
  if(df$Wmax[i] < wmax_0.9 & df$Vs[i] < Vs_0.9){
    curve(stabilizing_selection(x,ecotype_info$bio1_ecotype[i],df$Wmax[i],df$Vs[i]),from = -30,to = 50,col=adjustcolor("grey80",alpha.f = 0.5),add = T,lwd=1)
  }
}

for(i in 1:231){
  if(df$Vs[i] > Vs_0.9){
    curve(stabilizing_selection(x,ecotype_info$bio1_ecotype[i],df$Wmax[i],df$Vs[i]),from = -30,to = 50,col=adjustcolor("#CB181D",alpha.f = 0.5),add = T)
  }
}
for(i in 1:231){
  if(df$Wmax[i] > wmax_0.9){
    curve(stabilizing_selection(x,ecotype_info$bio1_ecotype[i],df$Wmax[i],df$Vs[i]),from = -30,to = 50,col=adjustcolor("#31A354",alpha.f = 0.5),add = T)
  }
}

#for(i in 1:31){
#  abline(v=site_climate[i,2],col="grey40",lty=3)
#  #text(x = site_climate[i,2],y=2.5,label=site_climate$site[i],cex = 0.5,col="red")
#}
dev.off()



################################################################################
### Adapted from MOI
### Plot the Gaussian model lines
################################################################################

library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(scales)

PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

################################################################################
## Raw data
file_edelta = 'data/delta_ecotype_freq.txt'
file_e0 = 'data/merged_ecotype_frequency.txt'
file_eid = 'data/greneNet_final_v1.0.fam'
file_climate_sites= "grene/data/worldclim_sitesdata.rda"
file_info_sites=  "grene/data/sites.clim.rda"
files_climate_ecotypes="grene/data/worldclim_ecotypesdata.rda"
files_info_ecotypes="grene/data/ecotypes_data.rda"

load(files_climate_ecotypes)
load(file_climate_sites)
load(files_info_ecotypes)


## New results Sep 2024
Vsresults=read.table("data-intermediate/LOCAL_ADAPTATION_ecotype_Wmax_Vs.txt",header=T)
head(Vsresults)
Vsresults<-Vsresults %>% dplyr::rename(wmax=Wmax,Vs=Vs)


# Merge
tmp= merge(ecotypes_data, by.x="ecotypeid", worldclim_ecotypesdata, by.y="ecotypeid")
dfeco= merge(Vsresults, by.x="ecotype", tmp, by.y="ecotypeid")


################################################################################
# Plot a model example

lineprojection<-
  expand.grid(unique(dfeco$ecotype), seq(0,25,0.5)) %>% 
  rename(id=Var1, bio1transplant=Var2) %>% 
  merge(.,dfeco, by.x="id",by.y="ecotype") %>% 
  # mutate(Westimated=wmax-(1/Vs)*(bio1-bio1transplant)^2) # Gaussian model. incorrect wmax already in the real space
  mutate(Westimated=wmax * exp(( -(bio1-bio1transplant)^2) / Vs)) 
# dfeco$vulnerability <- dfeco$wmax * exp(-(dfeco$temperaturechange) / dfeco$Vs)

#redgreen<- c("#CB181D", "#FB6A4A", "#FCAE91","#FEE5D9","#EDF8E9","#A1D99B","#74C476","#006D2C")

# [1] "#F7FCF5" "#E5F5E0" "#C7E9C0" "#A1D99B" "#74C476" "#41AB5D" "#238B45" "#006D2C" "#00441B"
ggplot(lineprojection)+
  geom_line(aes(y=Westimated,x=bio1-bio1transplant,color=1/Vs, group=id))+
  scale_color_gradientn(colors=c(rev(RColorBrewer::brewer.pal(5,"Reds")),RColorBrewer::brewer.pal(5,"Greens")),
                        values = rescale(c(0,0.001,0.002,0.004,0.006,0.008,0.011,0.015,0.02)))+
  xlab("Temperature (C)")+
  ylab("W modeled")+
  theme_minimal()

fig_stabilizing_modellines<-
  ggplot(lineprojection)+
  geom_line(aes(y=Westimated,x=bio1transplant,color=1/Vs, group=id),alpha=0.5)+
  scale_color_gradientn(colors=c(rev(RColorBrewer::brewer.pal(5,"Reds")),RColorBrewer::brewer.pal(5,"Greens")),
                        values = rescale(c(0,0.001,0.002,0.004,0.006,0.008,0.011,0.015,0.02)))+
  xlab("Temperature (C)")+
  ylab("W modeled")+
  theme_minimal(base_size = 18)
fig_stabilizing_modellines

ggsave(plot = fig_stabilizing_modellines,filename = "figs/LOCAL_ADAPTATION_stabilizing_selection_curve.pdf",device = "pdf", height = 5,width = 6)

