## This script calculate the genomic offset for 30 sites using
## pearson, spearman and r2 using the climatic distance

rm(list=ls())

library(ggplot2)
library(dplyr)

setwd("~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/")

mydata <- read.delim("data-intermediate/LOCAL_ADAPTATION_stabilizing_selection_data_2024Jun18_collectionsite_era5_allyears_pt1_pt0.txt")
mydata <- mydata %>%
  filter(generation ==1) %>%
  filter(site != 33 )

unique_sites <- unique(mydata$site)
mydata$ecotype <- factor(mydata$ecotype,levels = unique(mydata$ecotype))

site_plot <- unique(paste(mydata$site,mydata$plot,sep="_"))

go_results <- as.data.frame(matrix(nrow=length(site_plot),ncol=6))
colnames(go_results) <- c("site","plot","pearson_r","sp_r","r2","source")
go_results$plot <- gsub(".*_","",site_plot)
go_results$site <- gsub("_.*","",site_plot)
go_results$source <- "climate_distance"

k=1

for(i in 1:length(unique_sites)){
  
  tmp_site <- mydata %>%
    filter(site==unique_sites[i])
  
  unique_plots <- unique(tmp_site$plot)
  
  for(j in 1:length(unique_plots)){
    
    tmp_plot <- tmp_site %>%
      filter(plot==unique_plots[j]) %>%
      mutate(bio1_sq_neg = -(bio1_ecotype - bio1_site)^2,
             log_p1_p0 = log(pt1 / pt0))
    
    go_results[k,3] <- cor(tmp_plot$bio1_sq_neg,exp(tmp_plot$log_p1_p0),method = "pearson")
    go_results[k,4] <- cor(tmp_plot$bio1_sq_neg,exp(tmp_plot$log_p1_p0),method = "spearman")
    go_results[k,5] <- summary(lm(exp(tmp_plot$log_p1_p0)~tmp_plot$bio1_sq_neg))$adj.r.squared
    
    k = k + 1
  }
}

write.table(go_results,"data-intermediate/Leave_one_out_prediction/climate_distance/Genomic_offset_site_climate_distance.txt",append = F,quote = F,sep = "\t",row.names = F)


