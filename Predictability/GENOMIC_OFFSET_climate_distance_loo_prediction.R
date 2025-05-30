rm(list=ls())

library(corrplot)
library(MCMCglmm)
library(r2glmm)
library(genio)
library(MASS)
library(grDevices)
library(dplyr)
library(patchwork)

## Genomic offset

## This script will compare the LOO prediction accuracy using the climate distance


prefix <- "Genomic_offset_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

## load the stabilizing selection table (era5)
mydata <- read.delim("data-intermediate/LOCAL_ADAPTATION_stabilizing_selection_data_2024Jun18_collectionsite_era5_allyears_pt1_pt0.txt",header = T) %>%
  filter(!site==33) %>%
  filter(generation == 1)
mydata$ecotype <- factor(mydata$ecotype,levels = unique(mydata$ecotype))
mydata$site <- factor(mydata$site,levels = unique(mydata$site))
unique_sites <- levels(mydata$site)

site_plot <- unique(paste(mydata$site,mydata$plot,sep="_"))

prediction_summary <- as.data.frame(matrix(nrow=length(site_plot),ncol=6))
colnames(prediction_summary) <- c("site","plot","pearson_r","sp_r","r2","source")
prediction_summary$site <- gsub("_.*","",site_plot)
prediction_summary$plot <- gsub(".*_","",site_plot)
prediction_summary$source <- "climate_distance"


for(i in 1:30){

  site_ <- as.numeric(unique_sites[i])
  index <- which(prediction_summary$site==site_)
  plot_ <- prediction_summary$plot[index]

  for(j in 1:length(plot_)){

    df <- mydata %>%
      filter(site == site_) %>%
      filter(plot == plot_[j]) %>%
      mutate(bio1_diff_sq = -(bio1_ecotype - bio1_site)^2)

    # p1/p0
    true <- (df$pt1)/df$pt0
    summary_index <- which(prediction_summary$site==site_ & prediction_summary$plot == plot_[j])

    ## pearson r
    prediction_summary[summary_index,3] <- cor(true,df$bio1_diff_sq)

    ## spearman rho
    prediction_summary[summary_index,4] <- cor(true,df$bio1_diff_sq,method = "spearman")

    ## R2
    prediction_summary[summary_index,5] <- summary(lm(true~df$bio1_diff_sq))$r.squared

  }
}

write.table(prediction_summary,file = paste0("data-intermediate/Leave_one_out_prediction/climate_distance/",prefix,"climate_distance_loo_prediction_summary.txt"),
            append = F,quote = F,sep = "\t",row.names = F)




