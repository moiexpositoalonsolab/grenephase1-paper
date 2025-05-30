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

## This script will compare the prediction accuracy using the mean frequency change
## across all the plots in the same site

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
prediction_summary$source <- "plot_mean"


for(i in 1:30){

  site_ <- as.numeric(unique_sites[i])
  index <- which(prediction_summary$site==site_)
  plot_ <- prediction_summary$plot[index]


  for(j in 1:length(plot_)){

    random_plot <- sample(plot_[plot_!=plot_[j]],size = 10,replace = T)
    mean_frequency_change <- mydata %>%
      filter(site == site_) %>%
      filter(plot %in% random_plot) %>%
      mutate(p_ratio = pt1/pt0) %>%
      group_by(ecotype) %>%
      summarise(mean = mean(p_ratio))

    df <- mydata %>%
      filter(site == site_) %>%
      filter(plot == plot_[j])

    true <- df$pt1/df$pt0
    summary_index <- which(prediction_summary$site==site_ & prediction_summary$plot == plot_[j])

    ## pearson r
    prediction_summary[summary_index,3] <- cor(true,mean_frequency_change$mean)

    ## spearman rho
    prediction_summary[summary_index,4] <- cor(true,mean_frequency_change$mean,method = "spearman")

    ## R2
    prediction_summary[summary_index,5] <- summary(lm(true~mean_frequency_change$mean))$r.squared

    # ## top10 accuracy
    # prediction_summary[summary_index,6] <- sum(order(true,decreasing = T)[1] %in% order(mean_frequency_change$mean,decreasing = T)[1:10] )
    #
    # ## bottom10 accuracy
    # prediction_summary[summary_index,7] <- sum(order(true,decreasing = F)[1] %in% order(mean_frequency_change$mean,decreasing = F)[1:10] )
    #
  }
}

write.table(prediction_summary,file = paste0("data-intermediate/Leave_one_out_prediction/replicates/",prefix,"plot_replication_prediction_summary.txt"),
            append = F,quote = F,sep = "\t",row.names = F)




