rm(list=ls())

library(ggplot2)
library(patchwork)
library(dplyr)

## EXPERIMENTAL_INTRO

## This script plot the climate PCA of 1TG ecotypes, 231 founders and 43 sites

prefix <- "EXPERIMENTAL_INTRO_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

samples <- read.delim("data/samples_data_fix57.csv",sep=",")

sum(samples$usesample==TRUE)

hist(samples$flowerscollected,breaks=50)
hist(unique(samples$weighted_mean_coverage),breaks=50)

p1 <- ggplot(samples)+
  geom_histogram(aes(flowerscollected),bins = 100,fill = "grey70",color="grey20")+
  xlab("Number of flowers collected")+
  theme_minimal()
p2 <- samples %>%
  filter(usesample==T) %>%
  ggplot(.)+
  geom_histogram(aes(weighted_mean_coverage),bins = 50,binwidth = 0.5,fill = "grey70",color="grey20")+
  xlab("Average sequencing coverage per population")+
  theme_minimal()
combined <- p1 | p2

ggsave(combined,file=paste0("figs/",prefix,"GrENE_net_flower_samples_sequencing_historgram.pdf"),device = "pdf",width = 8,height = 4)

table(samples$code)
