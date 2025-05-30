################################################################################
### Goal
### Plot replicate leave one out 

################################################################################
library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(patchwork)

################################################################################
#### Location
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
# myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

dat<-read_tsv(file="tables/leave_one_out_comparison_all.txt")
load("grene/data/worldclim_sitesdata.rda")
worldclim_sitesdata
sites_simple_names<-read.csv("grene/data/sites_simple_names.csv")

################################################################################
# Summarise replicate r

library(readr)
res<-
  dat %>% 
  dplyr::filter(source=='replicates') %>% 
  group_by(site) %>% 
  summarise(mean_r=mean(pearson_r),
            sd_r=sd(pearson_r)
            ) %>% 
  merge(.,worldclim_sitesdata,by="site") %>% 
  merge(.,sites_simple_names, by="site") 
  
head(res,2)

# Add a column for the rank of bio1
res <- res %>%
  mutate(rank_bio1 = rank(bio1))

# Create the plot
plot_rep<-
ggplot(res) +
  geom_point(aes(y = mean_r, x = reorder(city, bio1)), size=2,alpha=0.5) +
  theme_classic() +
  labs(x = "City (ranked bio1)", y = "Mean r (replicates)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot_rep

save_plot("figs/fig-predictability-replicate-leaveoneout-bycity.pdf",
          plot_rep, base_height = 2.5,base_width = 12)
save_plot("figs/fig-predictability-replicate-leaveoneout-bycity.png",
          plot_rep, base_height = 2.5,base_width = 12)

