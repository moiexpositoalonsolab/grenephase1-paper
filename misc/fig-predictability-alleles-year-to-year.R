################################################################################
### Goal
### Test predictability to next year


################################################################################
#### If transformed to collab can use the stuff below
# from google.colab import drive
# drive.mount('/content/drive')
# %load_ext rpy2.ipython
# %%R

################################################################################
#### Libraries
library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(patchwork)

################################################################################
#### Location if run locally
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
# myfolder<-"~/grenephase1-analyses"

myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
# write.csv( file= paste0(myfolder,"/data-intermediate/allele_allyears_frequencies_semilong_raw.csv"), 
#            allele_allyears_frequencies_semilong_raw)
load(file=paste0("data-intermediate/allele_allyears_frequencies_semilong_raw.rda"))

plot(a$y3 ~ a$y2)
cor.test(a$y3 , a$y2)


table(a$y2 >0, a$y1 >0)
table(a$y2 >0, a$y1 >0) %>% fisher.test()

asub<-a %>% 
  dplyr::filter(site==4)


table(asub$y2 >0, asub$y1 >0) %>% fisher.test()
plot(asub$y2, asub$y1 )

table(asub$y3 >0, asub$y2 >0) %>% fisher.test()
plot(asub$y3, asub$y2 )


table((asub$y3-asub$y2) >0, (asub$y2-asub$y1) >0) %>% fisher.test()

