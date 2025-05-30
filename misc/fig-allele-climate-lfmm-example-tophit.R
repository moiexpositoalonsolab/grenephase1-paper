################################################################################
### Goal
### plot of allele frequency by climate 

################################################################################
library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(patchwork)
# Define the custom theme function
theme_grenenet <- function(mysize=12) {
  theme_minimal() +
    theme(
      axis.title = element_text(size = mysize, color = "black"),
      axis.text = element_text(size = mysize, color = "black"),
      legend.text = element_text(size = mysize, color = "black"),
      legend.title = element_text(size = mysize, color = "black")
    )
}

################################################################################
#### Location
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
# myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

#### Allele data
# The intermediate dataset is produced by this Colab
# https://colab.research.google.com/drive/1u4nQh6-f7U_VcjBxQ2JgtJDEfyj5MVjA?authuser=1#scrollTo=3wn60Z9cyAbq

a<-read.csv("data-intermediate/allele_lfmm_bio1_intermediate_merge_climate.csv")
dim(a)
head(a)

################################################################################
# Assuming 'a' is your data frame and it has 1000 columns
set.seed(123)  # For reproducibility
random_columns <- sample(names(a), 10)  # Select 10 random columns

# Initial plot with the specified column X2_9902238
allele_lfmm_adaptive <- ggplot(a, aes(x = bio1)) +
  geom_point(aes(y = X2_9902238), color = "#74C476", size = 3, alpha = 0.5) +
  stat_summary(aes(y=X2_9902238,x=bio1),color="#74C476", size=2)+
  stat_smooth(aes(y = X2_9902238), method = 'glm', color = "#41AB5D") +
  labs(y = "Allele frequency change (p1-p0)", x = "Temperature (C)") +
  theme_grenenet(15)

# Add random columns to the plot
allele_lfmm_neutral<-
  ggplot(a, aes(x = bio1)) +
  labs(y = "Allele frequency change (p1-p0)", x = "Temperature (C)") +
  theme_grenenet(15)
for (col in random_columns) {
  allele_lfmm_neutral <- 
    allele_lfmm_neutral +
    geom_point(aes_string(y = col), color = "lightgrey", size = 3, alpha = 0.2) +
    stat_summary(aes_string(y = col), color = "lightgrey", size = 2, alpha = 0.2) +
    stat_smooth(aes_string(y = col), method = 'glm', color = "lightgrey")
}
allele_lfmm_neutral


save_plot(paste0("figs/fig-allele_climate_lfmm_neutral.pdf"),allele_lfmm_neutral,base_width = 5,base_height = 5)
save_plot(paste0("figs/fig-allele_climate_lfmm_adaptive_X2_9902238.pdf"),allele_lfmm_adaptive,base_width = 5,base_height = 5)
save_plot(paste0("figs/fig-allele_climate_lfmm_neutral.png"),allele_lfmm_neutral,base_width = 5,base_height = 5)
save_plot(paste0("figs/fig-allele_climate_lfmm_adaptive_X2_9902238.png"),allele_lfmm_adaptive,base_width = 5,base_height = 5)


