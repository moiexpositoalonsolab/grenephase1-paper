################################################################################
### Goal
### Estimate some basic metrics of ecotypes increasing in frequency across
### sites, means, variances, optimal places

################################################################################
library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(cowplot)

################################################################################
#### Location
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
# Load the long format with ecotype frequencies and climates
load("data-intermediate/ecotype_terminal_frequencies_long_raw_climate.rda")
load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize_flowers.rda")
load("grene/data/worldclim_sitesdata.rda")

e = ecotype_allyears_frequencies_long_raw_climate_popsize_flowers
et = ecotype_terminal_frequencies_long_raw_climate
et$freq=et$maxfreq
head(e)

################################################################################
# Check classic local adaptation as quality check

ggplot(e ) +
  geom_point(aes(x=sites_bio18-ecotypes_bio18,y=freq))+
  # stat_smooth(aes(x=avg_freq,y=var_freq), color = "darkgrey")+
  # stat_smooth(aes(x=avg_freq,y=var_freq),method = "glm",color = "darkgrey")+
  # xlab("Mean in frequency GrENE-net (mean)")+  ylab("CV in frequency GrENE-net")+
  theme_minimal()

ggplot(e ) +
  geom_point(aes(x=sites_bio18-ecotypes_bio18,
                 y=sites_bio1-ecotypes_bio1,
                 color=log10(freq+0.001/startfreq)))+
  theme_minimal()

################################################################################
# There are 3 things in these ecotypes,
# first the distance to origin, does it fit with where they perform?
# are the frequency values typically high or typically low?
# are the variance between accessions large or small? like fitness variance

SAVEAGAIN=F

# Create a subset to get things easier
es= et %>%
  dplyr::select(id,site,longitude,latitude,
                ecotypes_bio1, sites_bio1,
                ecotypes_bio12, sites_bio12,
                freq, startfreq)
head(es)

# Get some general info, like what is the per ecotype average
simplemetricsterminal <-
  es %>%
  dplyr::group_by(id,ecotypes_bio1,ecotypes_bio12) %>%
  summarise(avg_freq = mean(freq),
            var_freq = var(freq)/mean(freq),
            opt_clim = weighted.mean(sites_bio1,freq/sum(freq))
            )
if(SAVEAGAIN)save(file = "data-intermediate/simplesummaries-terminal.rda",simplemetricsterminal)

# Create a subset to get things easier
es= e %>%
  dplyr::select(id,site,longitude,latitude, ecotypes_bio1,ecotypes_bio12, sites_bio1,sites_bio12, freq, startfreq)
head(es)

# Get some general info, like what is the per ecotype average
simplemetricsallyears <-
  es %>%
  dplyr::group_by(id,ecotypes_bio1,ecotypes_bio12) %>%
  summarise(avg_freq = mean(freq),
            var_freq = var(freq)/mean(freq),
            opt_clim = weighted.mean(sites_bio1,freq/sum(freq))
  )
if(SAVEAGAIN)save(file = "data-intermediate/simplesummaries-allyears.rda",simplemetricsallyears)

# to easy writing
simplemetrics=simplemetricsallyears

################################################################################
### Simple summaries for each site
simplemetrics_persite <-
  e %>%
  dplyr::group_by(site,rep,year, sites_bio1) %>%
  reframe(avg_freq = mean(freq),
            var_freq = var(freq)/mean(freq),
          opt_clim = weighted.mean(ecotypes_bio1,(freq/sum(freq, na.rm=T)),na.rm=T),
          opt_clim_permuted = weighted.mean(ecotypes_bio1,(sample(freq)/sum(freq, na.rm=T)),na.rm=T),
          top_ecotype = ecotypes_bio1[freq==max(freq)]
  )
if(SAVEAGAIN)save(file = "data-intermediate/simplemetrics_persite.rda",simplemetrics_persite)

#####

ggplot(simplemetrics_persite)+
  geom_abline()+
  geom_jitter(aes(y=opt_clim, x= sites_bio1))

ggplot(simplemetrics_persite)+
  geom_hline(yintercept = 0)+
  # geom_point(aes(y=opt_clim-opt_clim_permuted, x= sites_bio1))+
  stat_summary(aes(y=opt_clim-opt_clim_permuted, x= sites_bio1))



################################################################################
# Visualize averages vs coeffiecient variation, and where they come from


# ggplot(simplemetrics ) +
#   geom_point(aes(x=log10(avg_freq),y=log(var_freq)) )+
#   # stat_smooth(aes(x=avg_freq,y=var_freq), color = "darkgrey")+
#   # stat_smooth(aes(x=avg_freq,y=var_freq),method = "glm",color = "darkgrey")+
#   xlab("Mean in frequency GrENE-net (mean)")+  ylab("CV in frequency GrENE-net")+
#   theme_minimal()

fig_ecotype_average_freq_by_CV_freq<-
ggplot(simplemetrics ) +
  geom_point(aes(x=avg_freq,y=var_freq))+
  stat_smooth(aes(x=avg_freq,y=var_freq), color = "darkgrey")+
  stat_smooth(aes(x=avg_freq,y=var_freq),method = "glm",color = "darkgrey")+
  xlab("Mean in frequency GrENE-net (mean)")+  ylab("CV in frequency GrENE-net")+
  theme_minimal()
fig_ecotype_average_freq_by_CV_freq
save_plot("figs/fig_ecotype_average_freq_by_CV_freq.pdf",fig_ecotype_average_freq_by_CV_freq, base_width = 9,base_height =7)
save_plot("figs/fig_ecotype_average_freq_by_CV_freq.png",fig_ecotype_average_freq_by_CV_freq, base_width = 9,base_height =7)


fig_ecotype_average_freq_by_CV_freq_log<-
  ggplot(simplemetrics ) +
  geom_point(aes(x=avg_freq,y=var_freq))+
  stat_smooth(aes(x=avg_freq,y=var_freq), color = "darkgrey")+
  stat_smooth(aes(x=avg_freq,y=var_freq),method = "glm",color = "darkgrey")+
  scale_x_log10()+
  scale_y_log10()+
  xlab("Mean in frequency GrENE-net (mean)")+
  ylab("CV in frequency GrENE-net")+
  theme_minimal()
fig_ecotype_average_freq_by_CV_freq_log
save_plot("figs/fig_ecotype_average_freq_by_CV_freq_log.pdf",fig_ecotype_average_freq_by_CV_freq_log, base_width = 9,base_height =7)
save_plot("figs/fig_ecotype_average_freq_by_CV_freq_log.png",fig_ecotype_average_freq_by_CV_freq_log, base_width = 9,base_height =7)



ggplot(simplemetricsallyears ,
       aes(y=avg_freq,x=ecotypes_bio12)) +
  geom_point()+
  stat_smooth(method='glm', color="lightgrey")+
  stat_smooth(method='glm',formula = y~poly(x,2), color="grey40")+
  ylab("Mean in frequency GrENE-net (mean)")+
  xlab("Annual precipitaiton (mm)")+
  theme_minimal()
fig_ecotype_average_freq_by_bio1<-
ggplot(simplemetrics ,
       aes(y=avg_freq,x=ecotypes_bio1)) +
  geom_point()+
  stat_smooth(method='glm', color="lightgrey")+
  stat_smooth(method='glm',formula = y~poly(x,2), color="grey40")+
  ylab("Mean in frequency GrENE-net (mean)")+
  xlab("Temperature of origin (C)")+
  theme_minimal()
fig_ecotype_average_freq_by_bio1

save_plot("figs/fig_ecotype_average_freq_by_bio1.pdf",fig_ecotype_average_freq_by_bio1, base_width = 9,base_height =7)
save_plot("figs/fig_ecotype_average_freq_by_bio1.png",fig_ecotype_average_freq_by_bio1, base_width = 9,base_height =7)

fig_ecotype_coeffvar_freq_by_bio1<-
ggplot(simplemetrics ,
       aes(y=var_freq,x=ecotypes_bio1)) +
  geom_point()+
  stat_smooth(method='glm', color="lightgrey")+
  stat_smooth(method='glm',formula = y~poly(x,2), color="grey40")+
  ylab("CV in frequency GrENE-net")+
  xlab("Temperature of origin (C)")+
  theme_minimal()

fig_ecotype_coeffvar_freq_by_bio1
save_plot("figs/fig_ecotype_coeffvar_freq_by_bio1.pdf",fig_ecotype_coeffvar_freq_by_bio1, base_width = 9,base_height =7)
save_plot("figs/fig_ecotype_coeffvar_freq_by_bio1.png",fig_ecotype_coeffvar_freq_by_bio1, base_width = 9,base_height =7)



fig_ecotype_weightopt_freq_by_bio1<-
ggplot(simplemetrics ,
       aes(y=opt_clim,x=ecotypes_bio1,color=avg_freq,size=avg_freq)) +
  geom_vline(xintercept= c(unique(es$sites_bio1)),lty='dotted' )+
  geom_hline(yintercept= c(unique(es$sites_bio1)),lty='dotted' )+
  geom_point(color='black')+
  geom_point()+
  scale_color_gradientn("Mean frequency \n across sites",colors=brewer.pal(9,"Greens"))+
  stat_smooth(method='glm')+
  stat_smooth(method='glm', formula = y~poly(x,2), col='darkgrey')+
  geom_abline(intercept = 0,slope = 1, col='darkgrey')+
  xlim(c(-5,+24))+
  ylim(c(-5,+24))+
  ylab("Temperature optimum from GrENE-net (C)")+
  xlab("Temperature of origin (C)")+
  guides(size = "none") +
  theme_minimal()

fig_ecotype_weightopt_freq_by_bio1
save_plot("figs/fig_ecotype_weighted_optimum_freq_by_bio1.pdf",fig_ecotype_weightopt_freq_by_bio1, base_width = 9,base_height =7)
save_plot("figs/fig_ecotype_weighted_optimum_freq_by_bio1.png",fig_ecotype_weightopt_freq_by_bio1, base_width = 9,base_height =7)


#### Add Vs
Vsresults<-read.table("data-intermediate/stabilizing_selection_wmax_vs.txt",header = T) %>%
  dplyr::rename(Vs=vs_mean, wmax=wmax_mean)
simplemetricsmerge<-merge(simplemetrics, by.x="id",Vsresults,by.y="ecotype")


  ggplot(simplemetricsmerge ,
         # aes(y=opt_clim,x=ecotypes_bio1,color=avg_freq,size=avg_freq)) +
      aes(y=opt_clim,x=ecotypes_bio1,color=1/Vs,size=wmax)) +
  geom_vline(xintercept= min(c(unique(es$sites_bio1))),lty='dotted' )+
  geom_vline(xintercept= max(c(unique(es$sites_bio1))),lty='dotted' )+
  geom_hline(yintercept= min(c(unique(es$sites_bio1))),lty='dotted' )+
  geom_hline(yintercept= max(c(unique(es$sites_bio1))),lty='dotted' )+
  geom_point()+
  scale_color_gradientn("1/Vs",colors=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ))+
  stat_smooth(method='glm',col='darkgrey')+
  geom_abline(intercept = 0,slope = 1, col='darkgrey', lty='dotted')+
  xlim(c(-5,+24))+
  ylim(c(-5,+24))+
  ylab("Temperature optimum from GrENE-net (C)")+
  xlab("Temperature of origin (C)")+
  theme_minimal()
fig_ecotype_weightopt_freq_by_bio1_simple

fig_ecotype_weightopt_freq_by_bio1_simple<-
ggplot(simplemetricsmerge ,
         aes(y=opt_clim,x=ecotypes_bio1,color=1/Vs,size=wmax, group=1/Vs > 0.015)) +
  geom_vline(xintercept= min(c(unique(es$sites_bio1))),lty='dotted' )+
  geom_vline(xintercept= max(c(unique(es$sites_bio1))),lty='dotted' )+
  geom_hline(yintercept= min(c(unique(es$sites_bio1))),lty='dotted' )+
  geom_hline(yintercept= max(c(unique(es$sites_bio1))),lty='dotted' )+
  geom_point()+
  scale_color_gradientn("1/Vs",colors=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ))+
  stat_smooth(method='glm',col='darkgrey')+
  geom_abline(intercept = 0,slope = 1, col='darkgrey', lty='dotted')+
  xlim(c(-5,+24))+
  ylim(c(-5,+24))+
    # xlim(c(min(c(unique(es$sites_bio1))),max(c(unique(es$sites_bio1)))))+
    # ylim(c(min(c(unique(es$sites_bio1))),max(c(unique(es$sites_bio1)))))+
  ylab("Temperature optimum from GrENE-net (C)")+
  xlab("Temperature of origin (C)")+
  theme_minimal()

ggplot(simplemetricsmerge ,
       aes(y=opt_clim,x=ecotypes_bio1,color=1/Vs,size=wmax, group=wmax > 2.2)) +
  geom_vline(xintercept= min(c(unique(es$sites_bio1))),lty='dotted' )+
  geom_vline(xintercept= max(c(unique(es$sites_bio1))),lty='dotted' )+
  geom_hline(yintercept= min(c(unique(es$sites_bio1))),lty='dotted' )+
  geom_hline(yintercept= max(c(unique(es$sites_bio1))),lty='dotted' )+
  geom_point()+
  scale_color_gradientn("1/Vs",colors=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ))+
  stat_smooth(method='glm',col='darkgrey')+
  geom_abline(intercept = 0,slope = 1, col='darkgrey', lty='dotted')+
  # xlim(c(-5,+24))+
  # ylim(c(-5,+24))+
  xlim(c(min(c(unique(es$sites_bio1))),max(c(unique(es$sites_bio1)))))+
  ylim(c(min(c(unique(es$sites_bio1))),max(c(unique(es$sites_bio1)))))+
  ylab("Temperature optimum from GrENE-net (C)")+
  xlab("Temperature of origin (C)")+
  theme_minimal()



save_plot("figs/fig_ecotype_weighted_optimum_freq_by_bio1-simplified.pdf",fig_ecotype_weightopt_freq_by_bio1_simple, base_width = 7,base_height =6)
save_plot("figs/fig_ecotype_weighted_optimum_freq_by_bio1-simplified.png",fig_ecotype_weightopt_freq_by_bio1_simple, base_width = 7,base_height =6)


cor.test(simplemetrics$ecotypes_bio1,simplemetrics$opt_clim)
sink("tables/test-correlation-ecotype-bio1-vs-weighted-optimum-bio1.txt")
print("Goal: Optimum bio1 per ecotype is calculated as the mean of grene-net site climate weighted by ecotype's terminal frequency at site")
print(cor.test(simplemetrics$ecotypes_bio1,simplemetrics$opt_clim))
sink()

# consider subseting to accessions that are specialist
sink("tables/test-adaptation-lag-weighted-optimum-bio1.txt")
print("Correlation between inferred optimum of accessions that are specialists (1/Vs>0.015) and not")
simplemetricsmerge %>%
  group_by(1/Vs>0.015) %>%
  summarise(r=cor.test(ecotypes_bio1,opt_clim)$estimate,
            p=cor.test(ecotypes_bio1,opt_clim)$p.value
            )
print("Optimum climate from grene-net - climate of origin. if negative, could indicate lag because the accessions performed better in colder climates than they come from")
print("all accessions")
summary(simplemetricsmerge$opt_clim - simplemetricsmerge$ecotypes_bio1)
# hist(simplemetricsmerge$opt_clim - simplemetricsmerge$ecotypes_bio1)
print("separating accessions with higher specialization >0.015")
simplemetricsmerge %>%
  group_by(1/Vs>0.015) %>%
  summarise(r=mean(opt_clim-ecotypes_bio1, na.rm=T))
print("separating accessions with higher specialization >0.13")
simplemetricsmerge %>%
  group_by(1/Vs>
             (max(1/simplemetricsmerge$Vs)+min(1/simplemetricsmerge$Vs))/2
           ) %>%
  summarise(r=mean(opt_clim-ecotypes_bio1, na.rm=T),
            rIQRup=quantile(opt_clim-ecotypes_bio1, 0.75, na.rm=T),
            rIQRlow=quantile(opt_clim-ecotypes_bio1, 0.25, na.rm=T),
  )

sink()


# ################################################################################
# # Group first by climate to get CV of frequency that means more like Vs
# simplemetrics_binclimate <-
#   es %>%
#   # dplyr::mutate(bio1breaks=cut(sites_bio1,10))+
#   dplyr::mutate(bio1breaks=round(sites_bio1/2)*2) %>%
#   dplyr::group_by(id,ecotypes_bio1, sites_bio1) %>%
#   dplyr::summarise(freq=mean(freq)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(id,ecotypes_bio1) %>%
#   summarise(avg_freq = mean(freq),
#             var_freq = var(freq)/mean(freq),
#             opt_clim = weighted.mean(sites_bio1,freq/sum(freq))
#   )
#
#
# ggplot(simplemetrics_binclimate ,
#        aes(y=avg_freq,x=var_freq)) +
#   geom_point()+
#   stat_smooth(method='glm', color="lightgrey")+
#   stat_smooth(method='glm',formula = y~poly(x,2), color="grey40")+
#   ylab("Mean in frequency GrENE-net [mean]")+
#   xlab("Variance in frequency GrENE-net [CV]")
#
#
# ################################################################################

## How much general attributes across sites predict outcomes?
