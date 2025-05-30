################################################################################
### Goal
### Estimate some basic metrics of ecotypes increasing in frequency across
### sites, means, variances, optimal places. Including poulation size

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

################################################################################
# Load the long format with ecotype frequencies and climates

# load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate.rda")
# load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize.rda")

load(file=paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize.rda"))
load(file=paste0(myfolder,"/data-intermediate/ecotype_terminal_frequencies_long_raw_climate_popsize.rda"))


e<-enew <- ecotype_terminal_frequencies_long_raw_climate_popsize %>%
            rename(year=max_year, freq=maxfreq) # Use the terminal
e<-enew <- ecotype_allyears_frequencies_long_raw_climate_popsize
enew[1:5,1:10]

#####*********************************************************************######
# Quick correlation
cor.test(log(e$freq/e$startfreq), (e$sites_bio1-e$ecotypes_bio1)^2,
         method='s')


################################################################################
# 1D stabilizing selection to check
bio1stabilizing<-
  e %>%
  arrange((freq-startfreq)) %>%
  ggplot(. ) +
  geom_point(aes(x=sites_bio1-ecotypes_bio1,
                 y=(freq-startfreq),
                 color=(freq-startfreq), size=(freq-startfreq)),
             alpha=0.9, shape=16)+
  scale_color_gradientn("freq. change (p1-p0)",colors=brewer.pal(9,"Greens")[-1])+
  scale_size_continuous("freq. change (p1-p0)", range = c(1,4))+
  geom_hline(yintercept = 0,lty='dotted')+
  geom_vline(xintercept = 0,lty='dotted')+
  theme_minimal()+
  labs(x = "Temperature transplant (C)", y = "Change in frequency (p1-p0)")
bio1stabilizing
save_plot("figs/fig_stabilizing_bio1.pdf",bio1stabilizing, base_width = 7,base_height =5)
save_plot("figs/fig_stabilizing_bio1.png",bio1stabilizing, base_width = 7,base_height =5)

bio1stabilizing_blank<-
  e %>%
  arrange((freq-startfreq)) %>%
  ggplot(. ) +
  geom_point(aes(x=sites_bio1-ecotypes_bio1,
                 y=(freq-startfreq),
                 color=(freq-startfreq), size=(freq-startfreq)),
             alpha=0, shape=16)+
  scale_color_gradientn("freq. change (p1-p0)",colors=brewer.pal(9,"Greens")[-1])+
  scale_size_continuous("freq. change (p1-p0)", range = c(1,4))+
  geom_hline(yintercept = 0,lty='dotted')+
  geom_vline(xintercept = 0,lty='dotted')+
  theme_minimal()+
  labs(x = "Temperature transplant (C)", y = "Change in frequency (p1-p0)")
save_plot("figs/fig_stabilizing_bio1-blank.pdf",bio1stabilizing_blank, base_width = 7,base_height =5)
save_plot("figs/fig_stabilizing_bio1-blank.png",bio1stabilizing_blank, base_width = 7,base_height =5)
################################################################################
# bio 18

bio18stabilizing<-
e %>%
  arrange((freq-startfreq)) %>%
  ggplot(. ) +
  geom_point(aes(x=sites_bio18-ecotypes_bio18,
                 y=(freq-startfreq)),
             color="black",
             alpha=0.9, shape=16)+
  geom_point(aes(x=sites_bio18-ecotypes_bio18,
                 y=(freq-startfreq),
                 color=(freq-startfreq), size=(freq-startfreq)),
             alpha=0.9, shape=16)+
  # scale_color_gradientn("freq. change (p1-p0)",colors=brewer.pal(9,"Greens")[-1])+
  scale_color_gradientn("freq. change (p1-p0)",colors=redgreen[-1])+
  scale_size_continuous("freq. change (p1-p0)", range = c(1,3))+
  geom_hline(yintercept = 0,lty='dotted')+
  geom_vline(xintercept = 0,lty='dotted')+
  theme_minimal()+
  labs(x = "Precipitation transplant (mm)", y = "Change in frequency (p1-p0)")
bio18stabilizing
save_plot("figs/fig_stabilizing_bio18.pdf",bio18stabilizing, base_width = 7,base_height =6)
save_plot("figs/fig_stabilizing_bio18.png",bio18stabilizing, base_width = 7,base_height =6)


e %>%
  arrange((freq-startfreq)) %>%
  ggplot(. ) +
  geom_point(aes(x=(sites_bio18-ecotypes_bio18)^2,
                 y=(freq-startfreq)),
             color="black",
             alpha=0.9, shape=16)+
  geom_point(aes(x=(sites_bio18-ecotypes_bio18)^2,
                 y=(freq-startfreq),
                 color=(freq-startfreq), size=(freq-startfreq)),
             alpha=0.9, shape=16)+
  # scale_color_gradientn("freq. change (p1-p0)",colors=brewer.pal(9,"Greens")[-1])+
  scale_color_gradientn("freq. change (p1-p0)",colors=redgreen[-1])+
  scale_size_continuous("freq. change (p1-p0)", range = c(1,2))+
  geom_hline(yintercept = 0,lty='dotted')+
  geom_vline(xintercept = 0,lty='dotted')+
  theme_minimal()+
  labs(x = "Precipitation transplant (mm)", y = "Change in frequency (p1-p0)")

################################################################################
# 2D stabilizing selection
twodimensionstabilizing<-
  e %>%
    arrange((freq-startfreq)) %>%
    ggplot(. ) +
    # geom_point(aes(x=sites_bio12-ecotypes_bio12,
    #                y=sites_bio1-ecotypes_bio1,
    #                size=freq-startfreq),
    #            color='black',shape=16)+
    geom_point(aes(x=sites_bio12-ecotypes_bio12,
                   y=sites_bio1-ecotypes_bio1,
                   color=(freq-startfreq), size=(freq-startfreq)),
               alpha=0.9, shape=16)+
    scale_color_gradientn("freq. change (p1-p0)",colors=brewer.pal(9,"Greens")[-1])+
  scale_size_continuous("freq. change (p1-p0)", range = c(1,4))+
    geom_hline(yintercept = 0,lty='dotted')+
    geom_vline(xintercept = 0,lty='dotted')+
    theme_minimal()+
    labs(x = "Precipitation distance (mm)", y = "Temperature distance (°C)")
twodimensionstabilizing
save_plot("figs/fig_stabilizing_D2_bio1_bio12.pdf",twodimensionstabilizing, base_width = 6,base_height =5)
save_plot("figs/fig_stabilizing_D2_bio1_bio12.png",twodimensionstabilizing, base_width = 6,base_height =5)

twodimensionstabilizing<-
  e %>%
  arrange((freq-startfreq)) %>%
  ggplot(. ) +
  # geom_point(aes(x=sites_bio12-ecotypes_bio12,
  #                y=sites_bio1-ecotypes_bio1,
  #                size=freq-startfreq),
  #            color='black',shape=16)+
  geom_point(aes(x=sites_bio18-ecotypes_bio18,
                 y=sites_bio1-ecotypes_bio1,
                 color=(freq-startfreq), size=(freq-startfreq)),
             alpha=0.9, shape=16)+
  scale_color_gradientn("freq. change (p1-p0)",colors=brewer.pal(9,"Greens")[-1])+
  scale_size_continuous("freq. change (p1-p0)", range = c(1,4))+
  geom_hline(yintercept = 0,lty='dotted')+
  geom_vline(xintercept = 0,lty='dotted')+
  theme_minimal()+
  labs(x = "Summer precipitation distance (mm)", y = "Temperature distance (°C)")
twodimensionstabilizing
save_plot("figs/fig_stabilizing_D2_bio1_bio18.pdf",twodimensionstabilizing, base_width = 6,base_height =5)
save_plot("figs/fig_stabilizing_D2_bio1_bio18.png",twodimensionstabilizing, base_width = 6,base_height =5)

# Compose the 2D plot with

main_plot <- e %>%
  arrange((freq - startfreq)) %>%
  ggplot() +
  # geom_point(aes(x = sites_bio12 - ecotypes_bio12,
  #                y = sites_bio1 - ecotypes_bio1),
  #            color = 'black', shape = 16) +
  geom_point(aes(x = sites_bio12 - ecotypes_bio12,
                 y = sites_bio1 - ecotypes_bio1,
                 color = (freq - startfreq), size = (freq - startfreq)),
             alpha = 0.9, shape = 16) +
  scale_color_gradientn(colors = brewer.pal(9, "Greens")[-1]) +
  scale_size_continuous(range =  c(1,3))+
  geom_hline(yintercept = 0, lty = 'dotted') +
  geom_vline(xintercept = 0, lty = 'dotted') +
  theme_minimal()
main_plot_nolegend<-
  main_plot+
  labs(x = "", y = "")+
  theme(legend.position = "none")
# main_plot

# Side plot 1
side_plot1 <- ggplot(e, aes(x = sites_bio1 - ecotypes_bio1, y = freq - startfreq)) +
  geom_point(aes(color = (freq - startfreq), size = (freq - startfreq))) +
  theme_minimal() +
  scale_size_continuous(range =  c(1,3))+
  scale_color_gradientn(colors = brewer.pal(9, "Greens")[-1]) +
  # geom_hline(yintercept = 0,lty='dotted')+
  geom_vline(xintercept = 0,lty='dotted')+
  # labs(x = "sites_bio1 - ecotypes_bio1", y = "freq - startfreq")
  labs(x = "", y = "")
side_plot1<-
  side_plot1+ coord_flip()+
  theme(legend.position = "none")

# Side plot 2
side_plot2 <- ggplot(e, aes(x = sites_bio12 - ecotypes_bio12, y = freq - startfreq)) +
  geom_point(aes(color = (freq - startfreq), size = (freq - startfreq))) +
  theme_minimal() +
  scale_size_continuous(range =  c(1,3))+
  scale_color_gradientn(colors = brewer.pal(9, "Greens")[-1]) +
  # geom_hline(yintercept = 0,lty='dotted')+
  geom_vline(xintercept = 0,lty='dotted')+
  # labs(x = "sites_bio12 - ecotypes_bio12", y = "freq - startfreq")+
  labs(x = "", y = "")
side_plot2<-
  side_plot2+ #scale_y_reverse()+
  theme(legend.position = "none")

# Equation explanation plot
side_plot3<-
  ggplot(e, aes(x = (sites_bio1 - ecotypes_bio1)^2 ,
                y = log10(freq / startfreq),
                # color=log10(freq / startfreq)
                color=freq - startfreq
                )) +
  geom_point()+
  # geom_point(aes(color = (freq - startfreq), size = (freq - startfreq))) +
  theme_minimal() +
  scale_size_continuous(range =  c(1,3))+
  scale_color_gradientn(colors = brewer.pal(9, "Greens")[-1]) +
  stat_smooth(method='glm',color='grey')+
  # labs(x = "sites_bio12 - ecotypes_bio12", y = "freq - startfreq")
  labs(x = "", y = "")
side_plot3<-
  side_plot3+theme(legend.position = "none")

combined_plot<-
(side_plot2 + side_plot3 ) / ( main_plot_nolegend + side_plot1)


# combined_plot
save_plot("figs/fig_stabilizing_D2_bio1_bio12_multipanel.pdf",combined_plot, base_width = 9,base_height =9)
save_plot("figs/fig_stabilizing_D2_bio1_bio12_multipanel.png",combined_plot, base_width = 9,base_height =9)

save_plot("figs/fig_stabilizing_D2_bio1_bio12.pdf",main_plot, base_width = 9,base_height =9)
save_plot("figs/fig_stabilizing_D2_bio1_bio12.png",main_plot, base_width = 9,base_height =9)


################################################################################
# Do an enclidean distance
euclideandist <- lapply(1:19, function(i) {
  ( scale(e[,paste0("sites_bio",i)]) - scale(e[,paste0("ecotypes_bio",i)]))^2
  }) %>% do.call(cbind,.) %>% apply(.,1,sum)
e$euclideandist<-euclideandist

fig_euclidean<-
e %>%
  # dplyr::filter(year==3) %>%
ggplot(., aes(x = euclideandist,
              y = freq / startfreq)) +
  geom_point(aes(color = (freq / startfreq), size = (freq / startfreq))) +
  theme_minimal() +
  scale_size_continuous(range =  c(1,3))+
  scale_color_gradientn(colors = brewer.pal(9, "Greens")[-1]) +
  # labs(x = "sites_bio12 - ecotypes_bio12", y = "freq - startfreq")+
  labs(y = "Frequency change (p1/p0)", x = "Euclidean climatic distances (bio1-19)")

save_plot("figs/fig_stabilizing_euclidean_bio1-19.pdf",fig_euclidean, base_width = 7,base_height =5)
save_plot("figs/fig_stabilizing_euclidean_bio1-19.png",fig_euclidean, base_width = 7,base_height =5)

fig_euclidean+scale_y_log10()

# Add the equation
source("function-regression-utilities.R")
eq_euclidean<-
  lm_eq(log10(e$freq+0.001 / e$startfreq) , e$euclideandist )
eq_euclidean
save_plot("figs/fig_stabilizing_euclidean_bio1-19-equation.pdf",fig_euclidean, base_width = 7,base_height =5)
################################################################################
# Add the Vs to the modeling plot

tmpfig<-
  ggplot(e, aes(x = (sites_bio1 - ecotypes_bio1)^2 ,
                y = log10(freq / startfreq),
                # color=log10(freq / startfreq)
                color=freq - startfreq
  )) +
    geom_point()+
    # geom_point(aes(color = (freq - startfreq), size = (freq - startfreq))) +
    theme_minimal() +
    scale_size_continuous(range =  c(1,3))+
    scale_color_gradientn("freq change (pq-p0)", colors = brewer.pal(9, "Greens")[-1]) +
    stat_smooth(method='glm',color='grey')+
    # labs(x = "sites_bio12 - ecotypes_bio12", y = "freq - startfreq")
    labs(y = "log(p1/p0)", x = "Temperature distance (e^2)")

tmpfig2<-tmpfig


# Load the Vs results
Vsresults<-read.table("data-intermediate/stabilizing_selection_wmax_vs.txt",header = T) %>%
  rename(Vs=vs_mean, wmax=wmax_mean) %>%
  arrange(Vs)
Vsresults$rank<-round(rank(Vsresults$Vs))
colpal<-c(rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ) %>% rev
colpalet<-colorRampPalette(colpal)(nrow(Vsresults))

for(i in 1:nrow(Vsresults)){
  tmpfig2<-
  tmpfig2+
    geom_abline(intercept = log(Vsresults$wmax[i]),
                slope = - 1/Vsresults$Vs[i],
                col=colpalet[Vsresults$rank[i]], alpha=0.24)
}
tmpfig2

save_plot("figs/fig_stabilizing_bio1_model_logp1p0_distancesquare.pdf",tmpfig, base_width = 7,base_height =5)
save_plot("figs/fig_stabilizing_bio1_model_logp1p0_distancesquare.png",tmpfig, base_width = 7,base_height =5)
save_plot("figs/fig_stabilizing_bio1_model_logp1p0_distancesquare_1vs_slopes.pdf",tmpfig2, base_width = 7,base_height =5)
save_plot("figs/fig_stabilizing_bio1_model_logp1p0_distancesquare_1vs_slopes.png",tmpfig2, base_width = 7,base_height =5)
