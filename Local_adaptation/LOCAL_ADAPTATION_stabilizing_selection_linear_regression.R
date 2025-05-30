rm(list=ls())

library(ggplot2)
library(corrplot)
library(MCMCglmm)
library(r2glmm)
library(genio)
library(MASS)
library(grDevices)
library(dplyr)
library(patchwork)

## local adaptation

## This script fit the stabilizing selection model, a simple linear regression using first year data
## to see which bioclim variables are more explanatory 

prefix <- "LOCAL_ADAPTATION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

## load the era5 stabilizing selection table
mydata <- read.delim("data-intermediate/LOCAL_ADAPTATION_stabilizing_selection_data_2024Jun18_collectionsite_era5_allyears_pt1_pt0.txt",header = T)
mydata$ecotype <- factor(mydata$ecotype,levels = unique(mydata$ecotype))
mydata$site <- factor(mydata$site,levels = unique(mydata$site))
unique_sites <- levels(mydata$site)

df_bio <- mydata %>%
  filter(generation==1) %>% 
  filter(site != 33) %>%
  mutate(log_p_ratio = log(pt1/pt0),
         bio1_diff = (bio1_ecotype -bio1_site)^2,
         bio2_diff = (bio2_ecotype -bio2_site)^2,
         bio3_diff = (bio3_ecotype -bio3_site)^2,
         bio4_diff = (bio4_ecotype -bio4_site)^2,
         bio5_diff = (bio5_ecotype -bio5_site)^2,
         bio6_diff = (bio6_ecotype -bio6_site)^2,
         bio7_diff = (bio7_ecotype -bio7_site)^2,
         bio8_diff = (bio8_ecotype -bio8_site)^2,
         bio9_diff = (bio9_ecotype -bio9_site)^2,
         bio10_diff = (bio10_ecotype -bio10_site)^2,
         bio11_diff = (bio11_ecotype -bio11_site)^2,
         bio12_diff = (bio12_ecotype -bio12_site)^2,
         bio13_diff = (bio13_ecotype -bio13_site)^2,
         bio14_diff = (bio14_ecotype -bio14_site)^2,
         bio15_diff = (bio15_ecotype -bio15_site)^2,
         bio16_diff = (bio16_ecotype -bio16_site)^2,
         bio17_diff = (bio17_ecotype -bio17_site)^2,
         bio18_diff = (bio18_ecotype -bio18_site)^2,
         bio19_diff = (bio19_ecotype -bio19_site)^2) %>% 
  filter(!is.infinite(log_p_ratio) & !is.na(log_p_ratio))
df_bio <- df_bio[,49:68]
r2 <- c()
p <- c()
for(i in 2:20){
  a <- summary(lm(df_bio[,1]~df_bio[,i]))
  r2 <- c(r2, a$adj.r.squared)
  p <- c(p,a$coefficients[2,4])
}

df <- tibble(bioclim = paste0("bio",1:19),
       r2 = r2)
df$bioclim <- factor(df$bioclim,levels = paste0("bio",1:19))
ggplot(df,aes(x=bioclim,y=r2))+
  geom_bar( stat="identity" ) + theme_minimal(base_size = 15) +
  ylab("R2 \n (y ~ bio_distance)")+xlab("")

ggsave(filename =  paste0("figs/",prefix,"stabilizing_selection_bioclim_linear_regression.pdf"),device = "pdf",width = 10,height = 6,units = "in")



