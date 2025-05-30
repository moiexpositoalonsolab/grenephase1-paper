rm(list=ls())

library(ggplot2)
library(dplyr)

## local adaptation

## This script will use the estimated Vs to draw illustrative figure of higher Vs slower fitness decline when transplanting


prefix <- "LOCAL_ADAPTATION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

## load the dataframe
mydata <- read.delim("data-intermediate/LOCAL_ADAPTATION_stabilizing_selection_data_2024Jun18_collectionsite_era5_allyears_pt1_pt0.txt",header = T) %>%
  filter(generation == 1) %>%
  filter(site != 33)

stabilizing_parameters <- read.delim("data-intermediate/LOCAL_ADAPTATION_ecotype_Wmax_Vs.txt")
site_wavg <- read.delim("data-intermediate/LOCAL_ADAPTATION_site_wavg.txt")
df <- tibble(log_p_ratio = log(mydata$pt1/mydata$pt0),
             p_diff = mydata$pt1 - mydata$pt0,
             env_2 = (mydata$bio1_ecotype - mydata$bio1_site)^2,
             ecotype = mydata$ecotype,
             site = mydata$site) %>% 
  left_join(.,stabilizing_parameters,by="ecotype") %>%
  left_join(.,site_wavg,by="site") %>%
  filter(!is.infinite(log_p_ratio))

RColorBrewer::brewer.pal(9,"Greens")

ggplot(df,aes(x=env_2,y=log_p_ratio))+
  geom_point(color="#C7E9C0",pch=16,cex=3)+
  #population level
  geom_abline(slope = - 1/mean(df$Vs),intercept = mean(log(df$Wmax/df$wavg)),colour = "grey50",linetype = 1)+
  #ecotype with large Vs # 9542
  geom_abline(slope = - 1/df$Vs[df$ecotype==772],intercept = log(unique(df$Wmax[df$ecotype==772]) / unique(df$wavg[df$site==5])),colour = "salmon",linetype = 1)+
  #ecotype with small Vs #6961
  geom_abline(slope = - 1/df$Vs[df$ecotype==6961],intercept = log(unique(df$Wmax[df$ecotype==6961])/ unique(df$wavg[df$site==5])),colour = "#00441B",linetype = 1)+
  annotate(geom = "text",x = 900,y=0.5,label="GrENE-net founder j",color="salmon",cex=5)+
  annotate(geom = "text",x = 900,y=-2,label="GrENE-net population avg",color="grey50",cex=5)+
  annotate(geom = "text",x = 900,y=-5,label="GrENE-net founder i",color="#00441B",cex=5)+
  theme_minimal(base_size = 18)
ggsave(file = paste0("figs/",prefix,"ecotype_Vs_line_illustration.pdf"),device = "pdf",width = 6,height = 5)

