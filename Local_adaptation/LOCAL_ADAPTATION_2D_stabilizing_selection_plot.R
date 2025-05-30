rm(list=ls())

library(corrplot)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(gridExtra)
library(cowplot)
## local adaptation

## This script local adaptation 2D plot using terminal generation data

prefix <- "LOCAL_ADAPTATION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

## load the dataframe
p0 <- read.delim("data/founder_ecotype_frequency.txt",header=F)
colnames(p0) <- c("ecotype","freq")
mydata <- read.delim("data-intermediate/LOCAL_ADAPTATION_stabilizing_selection_data_2024Jun18_collectionsite_era5_allyears_pt1_pt0.txt",header = T)
#mydata <- read.delim("data-intermediate/stabilizing_selection_data_2024Jun18_collectionsite_era5_terminal_generation_pt_p0.txt",header = T)
mydata <- mydata %>%
  filter(site != 33) %>%
  left_join(.,p0,by="ecotype")

mydata$ecotype <- factor(mydata$ecotype,levels = unique(mydata$ecotype))
mydata$site <- factor(mydata$site,levels = unique(mydata$site))
unique_sites <- levels(mydata$site)


## Bio1 and Bio18
df <- mydata %>%
  filter(site != 33) %>%
  mutate(p_diff = mydata$pt1 - freq) %>%
  filter(!is.infinite(p_diff) & !is.na(p_diff)) %>%
  mutate(bio1_diff = as.numeric(bio1_ecotype - bio1_site)) %>%
  mutate(bio18_diff = as.numeric(bio18_ecotype - bio18_site)) %>%
  dplyr::select(p_diff,site,ecotype,bio1_diff,bio18_diff)


ggplot(df)+
  geom_point(aes(x=bio1_diff,y=p_diff,color=p_diff,size = p_diff))+
  scale_color_gradientn(colors=brewer.pal(9,"Greens")[-1])+
  scale_size_continuous(range =c(1,3))+
  geom_vline(xintercept = 0,linetype = "dotted")+
  theme_classic(base_size = 18)

stabilizing_2d <- df %>%
  arrange(c(p_diff)) %>%
  ggplot(.) +
  geom_point(aes(x=bio1_diff,y=bio18_diff,color=p_diff,size = p_diff),shape = 16)+
  scale_color_gradientn(colors=brewer.pal(9,"Greens")[-1])+
  scale_size_continuous(range =c(1,4))+
  geom_vline(xintercept = 0,linetype = "dotted")+
  geom_hline(yintercept = 0,linetype = "dotted")+
  xlab("Annual Mean Temp (C) \n (ecotype - site) ") +
  ylab("Precipitation in the Warmest Quarter (MM) \n (ecotype - site) ")+
  theme_classic(base_size = 18)
stabilizing_2d
bio1_hist <- ggplot(df)+
  geom_point(aes(x=bio1_diff,y=p_diff,color=p_diff,size = p_diff))+
  scale_color_gradientn(colors=brewer.pal(9,"Greens")[-1])+
  scale_size_continuous(range =c(1,4))+
  geom_vline(xintercept = 0,linetype = "dotted")+
  xlab(NULL)+
  ylab("delta ecotype frequency")+
  theme_classic(base_size = 18)
bio1_hist

bio18_hist <- ggplot(df)+
  geom_point(aes(x=bio18_diff,y=p_diff,color=p_diff,size = p_diff))+
  scale_color_gradientn(colors=brewer.pal(9,"Greens")[-1])+
  scale_size_continuous(range =c(1,4))+
  geom_vline(xintercept = 0,linetype = "dotted")+
  coord_flip()+
  xlab(NULL)+
  ylab("delta ecotype frequency")+
  theme_classic(base_size = 18)
bio18_hist
legend <- get_legend(stabilizing_2d)
stabilizing_2d <- stabilizing_2d + theme(legend.position = "none")
bio1_hist <- bio1_hist + theme(legend.position = "none")
bio18_hist <- bio18_hist + theme(legend.position = "none")

combined_plot <- grid.arrange(bio1_hist,legend,stabilizing_2d, bio18_hist, ncol = 2, nrow = 2)
combined_plot

ggsave(filename = paste0("figs/",prefix,"2D_stabilizing_selection_bio1_bio18_terminal_generation.pdf"),combined_plot, device = "pdf",width = 11,height = 10,units = "in")



## plot the 2D stabilizing selection for Bio1 and Bio12
df <- mydata %>%
  filter(site != 33) %>%
  mutate(p_diff = pt - p0) %>%
  filter(!is.infinite(p_diff) & !is.na(p_diff)) %>%
  mutate(bio1_diff = as.numeric(bio1_ecotype - bio1_site)) %>%
  mutate(bio12_diff = as.numeric(bio12_ecotype - bio12_site)) %>%
  dplyr::select(p_diff,site,ecotype,bio1_diff,bio12_diff)


ggplot(df)+
  geom_point(aes(x=bio1_diff,y=p_diff,color=p_diff))+
  scale_color_gradientn(colors=brewer.pal(9,"Greens")[-1])+
  scale_size_continuous(range =c(1,3))+
  geom_vline(xintercept = 0,linetype = "dotted")+
  theme_classic(base_size = 18)

stabilizing_2d <- df %>%
  arrange(c(p_diff)) %>%
  ggplot(.) +
  geom_point(aes(x=bio1_diff,y=bio12_diff,color=p_diff),shape = 16,cex=3)+
  scale_color_gradientn(colors=brewer.pal(9,"Greens")[-1])+
  #scale_size_continuous(range =c(1,4))+
  geom_vline(xintercept = 0,linetype = "dotted")+
  geom_hline(yintercept = 0,linetype = "dotted")+
  xlab(NULL)+
  ylab(NULL)+
  xlim(c(-33,20))+
  ylim(c(-2700,2350))+
  #xlab("Annual Mean Temp (C) \n (ecotype - site) ") +
  #ylab("Annual Percipitation (MM) \n (ecotype - site) ")+
  theme_classic(base_size = 18)
stabilizing_2d
bio1_hist <- ggplot(df)+
  geom_point(aes(x=bio1_diff,y=p_diff,color=p_diff),cex=3)+
  scale_color_gradientn(colors=brewer.pal(9,"Greens")[-1])+
  scale_size_continuous(range =c(1,4))+
  geom_vline(xintercept = 0,linetype = "dotted")+
  xlab(NULL)+
  ylab(NULL)+
  xlim(c(-33,20))+
  theme_classic(base_size = 18)
bio1_hist

bio12_hist <- ggplot(df)+
  geom_point(aes(x=bio12_diff,y=p_diff,color=p_diff),cex=3)+
  scale_color_gradientn(colors=brewer.pal(9,"Greens")[-1])+
  #scale_size_continuous(range =c(1,4))+
  geom_vline(xintercept = 0,linetype = "dotted")+
  coord_flip()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(c(-2700,2350))+
  theme_classic(base_size = 18)
bio12_hist
legend <- get_legend(stabilizing_2d)
stabilizing_2d <- stabilizing_2d + theme(legend.position = "none")
bio1_hist <- bio1_hist + theme(legend.position = "none")
bio12_hist <- bio12_hist + theme(legend.position = "none")

combined_plot <- grid.arrange(bio1_hist,legend,stabilizing_2d, bio12_hist, ncol = 2, nrow = 2)
combined_plot

ggsave(filename = paste0("figs/",prefix,"2D_stabilizing_selection_bio1_bio12_terminal_generation.pdf"),combined_plot, device = "pdf",width = 9,height = 8,units = "in")


