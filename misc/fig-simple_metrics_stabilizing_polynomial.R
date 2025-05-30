################################################################################
### Goal
### Estimate some basic metrics of ecotypes increasing in frequency across
### sites, means, variances, optimal places. Including poulation size
### Intead of MCMCglmm do spatial interpolation of frequency change

################################################################################
library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(cowplot)

transparent<-function (col, alpha = 0.5)
{
  res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255,
                                                c[3]/255, alpha))
  return(res)
}

################################################################################
#### Location
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
# myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
# Load the long format with ecotype frequencies and climates

load(file=paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize.rda"))
load(file=paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize.rda"))
load(file="data-intermediate/ecotype_terminal_frequencies_long_raw_climate_popsize.rda")
ls()

e<-enew <- ecotype_terminal_frequencies_long_raw_climate_popsize %>%
  rename(freq=maxfreq, year=maxyear)
e<-enew <- ecotype_allyears_frequencies_long_raw_climate_popsize
enew[1:5,1:10]
e[1:5,1:10]

# Load also Vs analysis


# Quantify simple metrics equivalent to Vs analysis
# Get some general info, like what is the per ecotype average
simplemetrics <-
  e %>%
  # dplyr::mutate(count=freq*flowerscollected_corrected) %>%
  dplyr::mutate(count=freq) %>%
  # dplyr::group_by(id,year) %>%
  dplyr::group_by(id) %>%
  summarise(avg_freq = mean(count),
            var_freq = var(count)/mean(count),
            opt_clim = weighted.mean(sites_bio1,count/sum(count))
  )
# add back all the other variables
simplemetricse<-
  # merge(simplemetrics,e,by=c("id","year")) #
  merge(simplemetrics,e,by=c("id")) #
################################################################################

# Plot
fig_optimal_bio1_lm<-
  simplemetricse %>%
  arrange(-opt_clim) %>%
  ggplot(.) +
    geom_point(aes(y=log(freq/startfreq),x=sites_bio1, group=id),
               shape=16, alpha=0.5,
               size=0.1
               )+
    stat_smooth(
        aes(y=log(freq/startfreq),x=sites_bio1, group=opt_clim,color=opt_clim),
        method = "glm",
        alpha=0.8,
        lwd=0.25,
        #color = transparent("darkgrey"),
        formula=y~poly(x,1), se=F)+
    scale_color_gradientn("Optimal climate (C)",colors = rev(brewer.pal(9,"RdBu")))+
    theme_minimal()+
    xlab("Temperature (C)")+ylab("log(p1/p0)")
fig_optimal_bio1_lm

fig_optimal_bio1_polylm<-
  ggplot(simplemetricse) +
  geom_point(aes(y=log(freq/startfreq),x=sites_bio1, group=id),
             shape=16, alpha=0.5,
             size=0.1
  )+
  stat_smooth(
    aes(y=log(freq/startfreq),x=sites_bio1, group=opt_clim,color=opt_clim),
    method = "glm",
    alpha=0.8,
    lwd=0.25,
    #color = transparent("darkgrey"),
    formula=y~poly(x,2), se=F)+
  scale_color_gradientn("Optimal climate (C)",colors = rev(brewer.pal(9,"RdBu")))+
  theme_minimal()+
  xlab("Temperature (C)")+ylab("log(p1/p0)")
fig_optimal_bio1_polylm

save_plot("figs/fig_ecotype_weighted_optimum_freq_by_bio1_lm.pdf",
          fig_optimal_bio1_lm, base_width = 7,base_height =5)
save_plot("figs/fig_ecotype_weighted_optimum_freq_by_bio1_lm.png",
          fig_optimal_bio1_lm, base_width = 7,base_height =5)
save_plot("figs/fig_ecotype_weighted_optimum_freq_by_bio1_polynomial.pdf",
          fig_optimal_bio1_polylm, base_width = 7,base_height =5)
save_plot("figs/fig_ecotype_weighted_optimum_freq_by_bio1_polynomial.png",
          fig_optimal_bio1_polylm, base_width = 7,base_height =5)


