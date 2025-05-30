#### If in google drive
# from google.colab import drive
# drive.mount('/content/drive')
# %load_ext rpy2.ipython

################################################################################

library(tidyverse)
library(dplyr)

library(ggplot2)
library(RColorBrewer)
library(tidyr)

################################################################################
# %%R
#### Location
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
# Load the long format with ecotype frequencies and climates
load("data-intermediate/ecotype_terminal_frequencies_long_raw.rda")
load("data-intermediate/ecotype_terminal_frequencies_long_raw_climate_popsize.rda")
load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize.rda")
load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize_flowers.rda")
load("data-external//atlas1001_phenotype_matrix_imputed_onlypheno_NEW.rda")

# Get the ecotype frequency long, rename
e = ecotype_terminal_frequencies_long_raw
e = ecotype_terminal_frequencies_long_raw_climate_popsize %>% 
  dplyr::rename(freq=maxfreq)
e = ecotype_allyears_frequencies_long_raw_climate_popsize
e = ecotype_allyears_frequencies_long_raw_climate_popsize_flowers
head(e)

# ###########
# Combine with phenotype data, subset first
phenosub <-
  pheno %>% 
  dplyr::select(id,
         Delta_13C,
         Leaf_Area,
         leafsize,
         FT10,
         FT16,
         RGR,
         DSDS50,
         Root_horizontal_index_avg,
         Root_growth_rate_avg,
         Relative_root_length
         )


colnames(pheno)[grep("root", colnames(pheno))]

# merge frequency and phenotypes and group to correlate per site freq and phenotype
freqpheno<-mmmtmp<-merge(e,phenosub, by="id") 


mmm<-merge(e,phenosub, by="id") %>% 
  group_by(across(c(site, starts_with("sites_"))), 
           latitude, longitude, altitude) %>% 
  dplyr::summarise(
    r_c13= cor(Delta_13C,freq),
    # r_c13= cor(Delta_13C,freq,method = "s"),
    r_ft10= cor(FT10,freq),
    r_ft16= cor(FT16,freq),
    r_la= cor(Leaf_Area,freq),
    r_roothorizontality= cor(Root_horizontal_index_avg,freq),
    r_rootgrowth= cor(Root_growth_rate_avg,freq),
    r_rootlength= cor(Relative_root_length,freq),
    r_dormancy= cor(DSDS50,freq),
    b_c13= coefficients(
      summary(lm(freq~Delta_13C+FT16 + RGR +leafsize +Leaf_Area +DSDS50 )) # correcting for many things
    )[2,1],
    p_c13= coefficients(
      summary(lm(freq~Delta_13C+FT16 + RGR +leafsize +Leaf_Area +DSDS50 )) # correcting for many things
    )[2,4],
    b_FT= coefficients(
      summary(lm(freq~FT16+Delta_13C + RGR +leafsize +Leaf_Area +DSDS50 )) # correcting for many things
    )[2,1],
    p_FT= coefficients(
      summary(lm(freq~FT16+Delta_13C + RGR +leafsize +Leaf_Area +DSDS50 )) # correcting for many things
    )[2,4],
    b_dormancy= coefficients(
      summary(lm(freq~DSDS50+Delta_13C + RGR +leafsize +Leaf_Area +FT16 )) # correcting for many things
    )[2,1],
    p_dormancy= coefficients(
      summary(lm(freq~DSDS50+Delta_13C + RGR +leafsize +Leaf_Area +FT16 )) # correcting for many things
    )[2,4]
    )



################################################################################
# Leaf area
mmm %>% 
  dplyr::filter(longitude > -15) %>% 
  ggplot(.,
         aes(x=sites_bio1, y=r_la, color=r_la)
  )+
  geom_point(size=5)+
  geom_hline(yintercept = 0,lty="dotted")+
  stat_smooth(method='glm',color='grey', se=F, formula = y~poly(x,2))+
  # scale_x_log10()+
  ylab("Leaf area(r)")+
  xlab("Annual temperature (C)")+
  theme_minimal()

# Plot dormancy
mmm %>% 
  dplyr::filter(longitude > -15) %>% 
  ggplot(.,
         aes(x=sites_bio1, y=r_dormancy, color=r_dormancy)
  )+
  geom_point(size=5)+
  scale_colour_gradientn(colors=redgreen)+
  geom_hline(yintercept = 0,lty="dotted")+
  stat_smooth(method='glm',color='grey', se=F, formula = y~poly(x,2))+
  # scale_x_log10()+
  ylab("Dormancy (DSDS50) (r)")+
  xlab("Annual temperature (C)")+
  theme_minimal()

# Plot root horizontality
mmm %>% 
  dplyr::filter(longitude > -15) %>% 
  ggplot(.,
         aes(x=sites_bio1, y=r_rootlength, color=r_rootlength)
  )+
  geom_point(size=5)+
  scale_colour_gradientn(colors=redgreen)+
  geom_hline(yintercept = 0,lty="dotted")+
  stat_smooth(method='glm',color='grey', se=F, formula = y~poly(x,2))+
  # scale_x_log10()+
  ylab("root h (r)")+
  xlab("Precipitation seasonality (mm)")+
  theme_minimal()
mmm %>% 
  dplyr::filter(longitude > -15) %>% 
  ggplot(.,
         aes(x=sites_bio15, y=r_rootgrowth, color=r_rootgrowth)
  )+
  geom_point(size=5)+
  # scale_colour_gradient2("root h (r)",
  #                        low = "#E3AF70",mid = "lightgrey",high = "#7B97CA",
  #                        midpoint = 0
  # )+
  geom_hline(yintercept = 0,lty="dotted")+
  stat_smooth(method='glm',color='grey', se=F, formula = y~poly(x,2))+
  # scale_x_log10()+
  ylab("root h (r)")+
  xlab("Precipitation seasonality (mm)")+
  theme_minimal()
mmm %>% 
  dplyr::filter(longitude > -15) %>% 
  ggplot(.,
         aes(x=sites_bio12, y=r_roothorizontality, color=r_roothorizontality)
  )+
  geom_point(size=5)+
  # scale_colour_gradient2("root h (r)",
  #                        low = "#E3AF70",mid = "lightgrey",high = "#7B97CA",
  #                        midpoint = 0
  # )+
  geom_hline(yintercept = 0,lty="dotted")+
  stat_smooth(method='glm',color='grey', se=F, formula = y~poly(x,2))+
  # scale_x_log10()+
  ylab("root h (r)")+
  xlab("Precipitation seasonality (mm)")+
  theme_minimal()


# Plot WUE for exploration
mmm %>% 
  dplyr::filter(longitude > -15) %>% 
  ggplot(.,
         aes(x=sites_bio18, y=b_c13, color=b_c13)
         )+
  geom_point(size=5)+
  scale_colour_gradient2("WUE (b)",
                         low = "#E3AF70",mid = "lightgrey",high = "#7B97CA",
                         midpoint = 0
  )+
  geom_hline(yintercept = 0,lty="dotted")+
  stat_smooth(method='glm',color='grey', se=F, formula = y~poly(x,2))+
  scale_x_log10()+
  ylab("WUE (b)")+
  xlab("Summer precipitation (mm)")+
  theme_minimal()


mmm %>% 
  dplyr::filter(longitude > -15) %>% 
ggplot(.,
       aes(x=sites_bio18, y=r_c13, color=r_c13)
       )+
  geom_point(size=5)+
  scale_colour_gradient2("WUE (r)",
                         low = "#E3AF70",mid = "lightgrey",high = "#7B97CA",
                         midpoint = 0
  )+
  stat_smooth(method='glm',color='grey', se=F)+
  geom_hline(yintercept = 0,lty="dotted")+
  ylab("WUE ()")+
  xlab("Summer precipitation (mm)")+
  scale_x_log10()+
  theme_minimal()


# Plot WUE for paper

p1wue<-
ggplot(mmm)+
  geom_jitter(aes(x=sites_bio12, y=sites_bio1, color=r_c13),
             size=5)+
  # geom_point(aes(x=sites_bio12, y=sites_bio1, color=r_c13),
  #            size=5)+
  scale_colour_gradient2("WUE (r)",
                        low = "#E3AF70",mid = "lightgrey",high = "#7B97CA",
                         midpoint = 0
                         )+
  # scale_colour_gradientn("WUE (r)",
  #                        colors=colorRampPalette( c("#9E490A","#C97B3A","#C97B3A","#C97B3A","#E3AF70","lightgrey","lightgrey","#7B97CA","#2B5DAD"))(5)
  # )+
  
  ylab("Annual temperature (C)")+
  xlab("Annual precipitation (mm)")+
  # theme_classic()
  theme_minimal()
p1wue


p2ft<-
ggplot(mmm)+
  geom_point(aes(x=sites_bio12, y=sites_bio1, color=r_ft10),
             size=5)+
  scale_colour_gradient2("Flowering (r)",
                         low = "#76AB33",mid = "lightgrey",high = "#798AE9",
                         midpoint = 0
  )+
  ylab("Annual temperature (C)")+
  xlab("Annual precipitation (mm)")+
  theme_classic()
p2ft

p3la<-
ggplot(mmm)+
  geom_point(aes(x=sites_bio12, y=sites_bio1, color=r_la),
             size=5)+
  scale_colour_gradient2("Leaf Area (r)",
                         low = "#4F5E42",mid = "lightgrey",high = "#76AB33",
                         midpoint = 0, position="bottom"
  )+
  ylab("Annual temperature (C)")+
  xlab("Annual precipitation (mm)")+
  theme_classic()
p3la

phenotypetrajectories<-
  plot_grid(p1wue,p2ft,p3la,
          nrow=1
          )
phenotypetrajectories
save_plot("figs/fig-phenotypes-ecotype_frequencies-climate.pdf",phenotypetrajectories,
          base_height = 3, base_width = 12)
save_plot("figs/fig-phenotypes-ecotype_frequencies-climate.png",phenotypetrajectories,
          base_height = 3, base_width = 12)




################################################################################
### Save cool plots potentially for main paper? with latitude

floweringlatitude<-
ggplot(mmm)+
  geom_point(aes(y=latitude, x=sites_bio1),color="white",
             size=5.2)+
  geom_point(aes(y=latitude, x=sites_bio1),color="black",
             size=5.1)+
  geom_point(aes(y=latitude, x=sites_bio1, color=r_ft10),
              size=5)+
  scale_colour_gradient2("FT10 (r)",
                         # low = "#E69D4E",mid = "lightgrey",high = "#728DC4",
                         low = "#76AB33",mid = "lightgrey",high = "#798AE9",
                         # low = "#76AB33",mid = "grey",high = "#798AE9",
                         # low = "#C97B3A",mid = "lightgrey",high = "#2E5EAB",
                         midpoint = 0
  )+
  ylab("Latitude (N)")+
  xlab("Annual temperature (C)")+
  # theme_classic()
  theme_minimal()
floweringlatitude

floweringlatitude
save_plot("figs/fig-phenotypes-ecotype_frequencies-climate-flowering-and-latitude.pdf",floweringlatitude,
          base_height = 4, base_width = 5)
save_plot("figs/fig-phenotypes-ecotype_frequencies-climate-flowering-and-latitude.png",floweringlatitude,
          base_height = 4, base_width = 5)

## Just correlation in xs

floweringcorrelation<-
  ggplot(mmm)+
  geom_point(aes(y=r_ft10, x=sites_bio1),color="white",
             size=5.2)+
  geom_point(aes(y=r_ft10, x=sites_bio1),color="black",
             size=5.1)+
  geom_point(aes(y=r_ft10, x=sites_bio1, color=r_ft10),
             size=5)+
  scale_colour_gradient2("FT10 (r)",
                         # low = "#E69D4E",mid = "lightgrey",high = "#728DC4",
                         low = "#76AB33",mid = "lightgrey",high = "#798AE9",
                         # low = "#76AB33",mid = "grey",high = "#798AE9",
                         # low = "#C97B3A",mid = "lightgrey",high = "#2E5EAB",
                         midpoint = 0
  )+
  ylab("Correlation ecotype change with flowering (r)")+
  xlab("Annual temperature (C)")+
  geom_hline(yintercept = 0,color='grey',lty='dotted')+
  # theme_classic()
  theme_minimal()
floweringcorrelation


save_plot("figs/fig-phenotypes-ecotype_frequencies-climate-flowering-correlation.pdf",floweringcorrelation,
          base_height = 4, base_width = 5)
save_plot("figs/fig-phenotypes-ecotype_frequencies-climate-flowering-correlation.png",floweringcorrelation,
          base_height = 4, base_width = 5)


floweringcorrelationsummmer<-
  ggplot(mmm)+
  geom_point(aes(y=r_ft10, x=sites_bio18),color="white",
             size=5.2)+
  geom_point(aes(y=r_ft10, x=sites_bio18),color="black",
             size=5.1)+
  geom_point(aes(y=r_ft10, x=sites_bio18, color=r_ft10),
             size=5)+
  scale_colour_gradient2("FT10 (r)",
                         # low = "#E69D4E",mid = "lightgrey",high = "#728DC4",
                         low = "#76AB33",mid = "lightgrey",high = "#798AE9",
                         # low = "#76AB33",mid = "grey",high = "#798AE9",
                         # low = "#C97B3A",mid = "lightgrey",high = "#2E5EAB",
                         midpoint = 0
  )+
  ylab("Correlation ecotype change with flowering (r)")+
  xlab("Summer precipitation (mm)")+
  geom_hline(yintercept = 0,color='grey',lty='dotted')+
  # theme_classic()
  theme_minimal()
floweringcorrelationsummmer


save_plot("figs/fig-phenotypes-ecotype_frequencies-climate-summerprecipitation-flowering-correlation.pdf",
          floweringcorrelationsummmer,
          base_height = 4, base_width = 5)
save_plot("figs/fig-phenotypes-ecotype_frequencies-climate-summerprecipitation-flowering-correlation.png",
          floweringcorrelationsummmer,
          base_height = 4, base_width = 5)

# Flowerig time against summer precipitation
mmm %>% 
  dplyr::filter(longitude > -15) %>% 
  ggplot(.,
         aes(x=sites_bio18, y=r_ft16, color=r_c13)
  )+
  geom_point(size=5.2, color='white')+
  geom_point(size=5.2, color='black')+
  geom_point(size=5)+
  scale_colour_gradient2("WUE (r)",
                         low = "#E3AF70",mid = "lightgrey",high = "#7B97CA",
                         midpoint = 0
  )+
  geom_hline(yintercept = 0,lty="dotted")+
  stat_smooth(method='glm',color='grey', se=F, formula = y~poly(x,2))+
  # scale_x_log10()+
  ylab("FT (r)")+
  xlab("Summer precipitation (mm)")+
  theme_minimal()

# Flowerig time against summer precipitation
mmm %>% 
  dplyr::filter(longitude > -15) %>% 
  ggplot(.,
         aes(x=sites_bio18, y=r_ft16, color=b_c13, shape=p_c13<0.05)
  )+
  geom_point(size=5.2, color='white')+
  geom_point(size=5.2, color='black')+
  geom_point(size=5)+
  scale_colour_gradient2("WUE (r)",
                         low = "#E3AF70",mid = "lightgrey",high = "#7B97CA",
                         midpoint = 0
  )+
  geom_hline(yintercept = 0,lty="dotted")+
  stat_smooth(method='glm',color='grey', se=F, formula = y~poly(x,2))+
  # scale_x_log10()+
  ylab("FT (r)")+
  xlab("Summer precipitation (mm)")+
  theme_minimal()


#### Now of WUE with summer precipitation

wueplotgradient<-
mmm %>% 
    dplyr::filter(longitude > -15) %>% 
    ggplot(.,
           aes(x=sites_bio18, y=b_c13, color=b_c13)
    )+
    geom_point(size=5.2, color='white')+
    geom_point(size=5.2, color='black')+
    geom_point(size=5)+
    scale_colour_gradient2("WUE (b)",
                           low = "#E3AF70",mid = "lightgrey",high = "#7B97CA",
                           midpoint = 0
    )+
    geom_hline(yintercept = 0,lty="dotted")+
    stat_smooth(method='glm',color='grey', se=F, formula = y~poly(x,2))+
    # scale_x_log10()+
    ylab("WUE (b)")+
    xlab("Summer precipitation (mm)")+
    theme_minimal()
wueplotgradient
save_plot("figs/fig-phenotypes-ecotype_frequencies-summerprecipitation-WUE-gradient.pdf",
          wueplotgradient,
          base_height = 4, base_width = 5)
save_plot("figs/fig-phenotypes-ecotype_frequencies-summerprecipitation-WUE-gradient.png",
          wueplotgradient,
          base_height = 4, base_width = 5)


wueplotcorrelation<-
  mmm %>% 
  dplyr::filter(longitude > -15) %>% 
  ggplot(.,
         aes(x=sites_bio18, y=r_c13, color=r_c13)
  )+
  geom_point(size=5.2, color='white')+
  geom_point(size=5.2, color='black')+
  geom_point(size=5)+
  scale_colour_gradient2("WUE (r)",
                         low = "#E3AF70",mid = "lightgrey",high = "#7B97CA",
                         midpoint = 0
  )+
  geom_hline(yintercept = 0,lty="dotted")+
  stat_smooth(method='glm',color='grey', se=F, formula = y~poly(x,2))+
  # scale_x_log10(limits=c(1,1000))+
  ylab("WUE (r)")+
  xlab("Summer precipitation (mm)")+
  theme_minimal()
wueplotcorrelation

save_plot("figs/fig-phenotypes-ecotype_frequencies-summerprecipitation-WUE-correlation.pdf",
          wueplotcorrelation,
          base_height = 4, base_width = 5)
save_plot("figs/fig-phenotypes-ecotype_frequencies-summerprecipitation-WUE-correlation.png",
          wueplotcorrelation,
          base_height = 4, base_width = 5)


################################################################################
# Evolution of phenotypes across landscapes
# The idea is that we can do a weighted sum of phenotypes based on the 
# ecotype frequencies

# founderphenotypes<-
#   
# pp<-
#   projected_data_df %>% 
#   mutate(sampleid=rownames(projected_data_df)) %>% 
#   mutate(sampleid=gsub( "X", "", sampleid)) %>% 
#   mutate(sample=sampleid) %>% 
#   separate(sampleid, c("site","year","rep"), "_") %>% 
#   merge(
#     dplyr::select(locations_data,longitude,latitude,altitude, site), 
#     by='site'
#   ) %>% 
#   merge(
#     .,
#     worldclim_sitesdata, by="site") %>% 
#   merge(
#     .,
#     sites_simple_names, by="site") %>% 
#   mutate(name=paste(city, country))

eco0<-read.table("data/founder_ecotype_frequency.txt", header=F) %>% 
  dplyr::rename(id=V1,freq=V2) %>% 
  merge(.,phenosub, by="id") 
  

# Get the average phenotype per site
founderphenotypes<-
  eco0 %>% 
  dplyr::summarise(
    FT_mean=weighted.mean(FT16,freq),
    WUE_mean=weighted.mean(Delta_13C,freq)
  )
  
sitephenotypes<-
freqpheno %>% 
  group_by(across(c(site, starts_with("sites_"))), 
           latitude, longitude, altitude) %>% 
  dplyr::summarise(
    FT_mean=weighted.mean(FT16,freq),
    WUE_mean=weighted.mean(Delta_13C,freq)
  )

# Do the plot of average phenotypes
# sitesplot<-
sitephenotypes %>% 
  dplyr::filter(longitude > -15) %>% 
  ggplot(.)+
  # geom_point(data=pheno,aes(y=FT16, x=Delta_13C),shape=3)+ # the ecotypes
  geom_point(aes(y=FT_mean, x=WUE_mean,  color=sites_bio1), size=3)+ # the sites
  geom_point(data=founderphenotypes,aes(y=FT_mean, x=WUE_mean),color="black", size=3)+
  xlab("WUE (DC13)")+ylab("Flowering time (days)")+
  scale_color_gradientn("Temperature (C)",colors = rev(brewer.pal(5,"RdBu")))+
  theme_minimal()
sitephenotypes %>% 
  dplyr::filter(longitude > -15) %>% 
  ggplot(.)+
  # geom_point(data=pheno,aes(y=FT16, x=Delta_13C),shape=3)+ # the ecotypes
  geom_point(aes(y=FT_mean, x=WUE_mean,  color=sites_bio18),size=3)+ # the sites
  geom_point(data=founderphenotypes,aes(y=FT_mean, x=WUE_mean),color="black", size=3)+
  xlab("WUE (DC13)")+ylab("Flowering time (days)")+
  scale_color_gradientn("Summer precipitation (mm)",colors = brewer.pal(5,"Blues"), trans="log10")+
  theme_minimal()
sitephenotypes %>% 
  dplyr::filter(longitude > -15) %>% 
  ggplot(.)+
  # geom_point(data=pheno,aes(y=FT16, x=Delta_13C),shape=3)+ # the ecotypes
  geom_point(aes(y=FT_mean, x=WUE_mean,  color=sites_bio15),size=5)+ # the sites
  geom_point(data=founderphenotypes,aes(y=FT_mean, x=WUE_mean),color="black", size=5)+
  xlab("WUE (DC13)")+ylab("Flowering time (days)")+
  scale_color_gradientn("Annual precipitation (mm)",colors = brewer.pal(5,"Spectral")[-1])+
  theme_minimal()

# sitesplot
