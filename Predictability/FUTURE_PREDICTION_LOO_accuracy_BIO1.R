rm(list = ls())
setwd("~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/data-intermediate/Leave_one_out_prediction/")
library(dplyr)

era5 <- read.delim("stabilizing_selection_era5/Genomic_offset_stabilizing_selection_loo_era5_prediction_summary.txt")
site_climate <- read.delim("../../data-external/bioclimvars_sites_era5_year_2018.csv",sep = ",")
#site_climate <- read.delim("../../data-external/bioclimvars_experimental_sites_era5.csv",sep = ",")

ecotype_climate <- read.delim("../../data-external/bioclimvars_ecotypes_era5.csv",sep=",")

population_center <- median(ecotype_climate$bio1_new)

era5 <- era5 %>%
  left_join(.,site_climate,by="site")


ggplot()+
  geom_point(data=ecotype_climate,mapping = aes(x=bio1_new,y=bio12_new),col="green")+
  geom_point(data=era5,mapping = aes(x=bio1,y=bio12,cex=r2),col="red2")
  

era5$site <- factor(era5$site,levels = unique(era5$site[order(era5$bio1)]))

summary(lm(r2~poly(bio1,2,raw = T),data = era5))

ggplot(era5)+
  geom_point(aes(x=bio1,y=r2),cex=2)+
  geom_smooth(aes(x=bio1,y=r2),method="lm",formula = y~poly(x,2,raw = T))+
  #geom_smooth(data=era5[era5$bio1<=17,],aes(x=bio1,y=r2),method="lm",formula = y~x)+
  #geom_smooth(data=era5[era5$bio1>=17,],aes(x=bio1,y=r2),method="lm",formula = y~x)+
  geom_vline(xintercept = 9.544582,colour = "red",linetype = 2,linewidth = 1)+
  xlab("Experimental Garden BIO1")+
  ylab(expression(LOO~prediction~accuracy~(r^2)))+
  annotate(geom = "text",x = 10,y=0.25,label= "R^2==0.126 ~ P==1.356e-10",parse=T,cex=5)+
  theme_minimal(base_size=18)

ggsave("../../figs/FUTURE_PREDICTION_LOO_BIO1.pdf",device = "pdf",height = 7,width = 7,units = "in")
ggsave("../../figs/FUTURE_PREDICTION_LOO_BIO1.png",device = "png",height = 7,width = 7,units = "in",dpi = 300)




ggplot(era5)+
  geom_point(aes(x=bio1,y=sp_r),cex=2)+
  geom_smooth(aes(x=bio1,y=sp_r),method="lm",formula = y~poly(x,2,raw = T))+
  #geom_smooth(data=era5[era5$bio1<=17,],aes(x=bio1,y=r2),method="lm",formula = y~x)+
  #geom_smooth(data=era5[era5$bio1>=17,],aes(x=bio1,y=r2),method="lm",formula = y~x)+
  geom_vline(xintercept = 9.544582,colour = "red",linetype = 2,linewidth = 1)+
  xlab("Experimental Garden BIO1")+
  ylab(expression(LOO~prediction~accuracy~(spr^2)))+
  #annotate(geom = "text",x = 10,y=0.25,label= "R^2==0.125 ~ P==1.356e-6",parse=T,cex=5)+
  theme_minimal(base_size=18)
summary(lm(era5$sp_r~poly(era5$bio1,2,raw = T)))
