## This script plots the Vs Wmax and stabilizing selection curves
library(dplyr)
library(ggplot2)
library(RColorBrewer)


setwd("~/Dropbox/Carnegie_Stanford/Projects/GrENE-net/stabilizing_selection_model/")
rm(list = ls())

df <- read.delim("stabilizing_selection_wmax_vs.txt")

ggplot(df,aes(x=wmax_mean,y=vs_mean,color=vs_sd))+
  geom_point(shape=16,cex=2.5)+
  scale_color_gradientn(colors = rev(brewer.pal(4,"Reds")))+
  theme_classic()

ggsave("bio1_wmax_vs.pdf",device = "pdf",units = "in",width = 5,height = 5)
ggsave("bio1_wmax_vs.png",device = "png",dpi = 800,width = 5,height = 5)


## draw the stabilizing selection curves for each ecotype

ecotype_climate <- read.delim("../metadata/worldclim_ecotypesdata_sorted_20240517.csv",sep="\t")
#predicted_ecotype_climate <- read.delim("../meta_information/predicted_climate_grenet_accession.txt",header=F)
site_climate <- read.delim("../metadata/bioclimvars_experimental_sites_era5.csv",sep=",")
meta <- read.delim("../merged_frequency/merged_sample_table.csv",sep=",")

site_climate <- site_climate %>%
  filter(site %in% meta$site) %>%
  arrange(match(site,unique(meta$site)))

Vs_0.95 <- quantile(df$vs_mean,0.95)
wmax_0.95 <- quantile(df$wmax_mean,0.95)

stabilizing_selection <- function(x,z0,wmax,vs){ wmax * exp( - (x-z0)^2 / vs)} 

pdf("stabilizing_curve_bio1.pdf",width = 7,height = 5)
plot(0,0,xlim=c(-30,50),ylim=c(0,2.5),type="n")
for(i in 1:231){
  if(df$wmax_mean[i]< wmax_0.95 & df$vs_mean[i] < Vs_0.95){
    curve(stabilizing_selection(x,ecotype_climate$bio1[i],df$wmax_mean[i],df$vs_mean[i]),from = -30,to = 50,col=adjustcolor("grey80",alpha.f = 0.3),add = T)
  }
}

for(i in 1:231){
  if(df$vs_mean[i] > Vs_0.95){
    curve(stabilizing_selection(x,ecotype_climate$bio1[i],df$wmax_mean[i],df$vs_mean[i]),from = -30,to = 50,col=adjustcolor("#A50F15",alpha.f = 0.5),add = T)
  }
}
for(i in 1:231){
  if(df$wmax_mean[i] > wmax_0.95){
    curve(stabilizing_selection(x,ecotype_climate$bio1[i],df$wmax_mean[i],df$vs_mean[i]),from = -30,to = 50,col=adjustcolor("#006D2C",alpha.f = 0.5),add = T)
  }
}

for(i in 1:31){
  abline(v=site_climate[i,2],col="grey40",lty=3)
  #text(x = site_climate[i,2],y=2.5,label=site_climate$site[i],cex = 0.5,col="red")
}
dev.off()

png("stabilizing_curve_bio1.png",width = 7,height = 5,units = "in",res = 800)
plot(0,0,xlim=c(-30,50),ylim=c(0,2.5),type="n")
for(i in 1:231){
  if(df$wmax_mean[i]< wmax_0.95 & df$vs_mean[i] < Vs_0.95){
    curve(stabilizing_selection(x,ecotype_climate$bio1[i],df$wmax_mean[i],df$vs_mean[i]),from = -30,to = 50,col=adjustcolor("grey80",alpha.f = 0.3),add = T)
  }
}

for(i in 1:231){
  if(df$vs_mean[i] > Vs_0.95){
    curve(stabilizing_selection(x,ecotype_climate$bio1[i],df$wmax_mean[i],df$vs_mean[i]),from = -30,to = 50,col=adjustcolor("#A50F15",alpha.f = 0.5),add = T)
  }
}
for(i in 1:231){
  if(df$wmax_mean[i] > wmax_0.95){
    curve(stabilizing_selection(x,ecotype_climate$bio1[i],df$wmax_mean[i],df$vs_mean[i]),from = -30,to = 50,col=adjustcolor("#006D2C",alpha.f = 0.5),add = T)
  }
}

for(i in 1:31){
  abline(v=site_climate[i,2],col="grey40",lty=3)
  #text(x = site_climate[i,2],y=2.5,label=site_climate$site[i],cex = 0.5,col="red")
}
dev.off()

