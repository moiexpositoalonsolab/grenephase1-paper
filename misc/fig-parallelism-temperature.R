##### Goal replot the parallelism  plot from Xing for publication

#####
# library(grene)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(tidyr)
theme_set(theme_cowplot())

myfolder<-"~/grenephase1-analyses/"
# setwd("~/Google Drive/RESEARCH/MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/figures-raw/plot_ecotype_trend/")
setwd(myfolder)

################################################################################
# Load parallelism
para<-read.table("data-intermediate/generation_1_parallelism.txt", header = T)
head(para)

# Load climate of sites
load("grene/data/worldclim_sitesdata.rda")
worldclim_sitesdata

# Merged

m<- merge(para, worldclim_sitesdata,by="site")
head(m)

################################################################################
gen1parallel<-
m %>% 
  filter(generation==1) %>% 
ggplot(.)+
  geom_point(aes(x = bio1, y=mean, color=source, shape=lower<0.01), size=5)+
  geom_segment(aes(x = bio1, xend=bio1,y=lower, yend=upper, color=source, shape=lower<0.01), size=1)+
  
  scale_color_manual(values = c("snp"="#CCCCCC","ecotype"="#000000"))+
  scale_shape_manual(values = c("TRUE"=1,"FALSE"=16))+
  stat_smooth(aes(x = bio1, y=mean, color=source, shape=lower<0), method="glm",se=F)+
  xlab("Annual temperature (C)") +
  ylab("Evolutionary parallelism (r)")+
  theme_classic()

gen1parallel

save_plot(
  "figs/fig-parallelism-temperature.pdf",
  gen1parallel,
  base_width = 5,base_height = 4 )

tmp<-
m %>% 
  filter(generation==1) %>% 
  filter(source=="ecotype") %>% 
  filter(lower > 0) 
cor.test(tmp$mean, tmp$bio1)
cor.test(tmp$mean, tmp$bio1, method = "s")

tmp<-
  m %>% 
  filter(generation==1) %>% 
  filter(source=="snp") %>% 
  filter(lower > 0) 
cor.test(tmp$mean, tmp$bio1)
cor.test(tmp$mean, tmp$bio1, method = "s")


m %>% 
  filter(generation==1) %>% 
  filter(lower<0.01) %>% 
  select(site, source)
