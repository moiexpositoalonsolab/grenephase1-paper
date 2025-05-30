#### Summarize phenotype evolution in paper

#####*********************************************************************######

library(tidyverse)
library(dplyr)

library(ggplot2)
library(RColorBrewer)
library(tidyr)


#### Location
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

#####*********************************************************************######
#### Read files #####

mmm<-df<-
  read_tsv(file = "tables/table_selection_coefficients.tsv")
head(mmm)



#####*********************************************************************######
##### Exploration of coef with climate #####

source("function-colors-grenenet.R")
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
  xlab("Annual temperature (C)")+
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

#### Correlation plots vs climate ariable ####


phenoplot1<-
  mmm %>%
  dplyr::filter(longitude > -15) %>%
  ggplot(.,
         aes(x=sites_bio18, y=r_ft10, color=r_ft10, shape=year)
  )+
  geom_point(size=5)+
  scale_colour_gradient2("Freq <-> Flowering (r)",
                         low = "#76AB33",mid = "lightgrey",high = "#798AE9",
                         midpoint = 0
  )+
  stat_smooth(method='glm',color='grey', se=F)+
  geom_hline(yintercept = 0,lty="dotted")+
  ylab("Correlation accession frequency with accession flowering (r)")+
  xlab("Summer precipitation (mm)")+
  scale_x_log10()+
  theme_minimal()
phenoplot1
save_plot("figs/fig-phenotypes-ecotype_frequencies-FT10-bio18.pdf",
          phenoplot1,
          base_height = 4, base_width = 4)
save_plot("figs/fig-phenotypes-ecotype_frequencies-FT10-bio18.png",
          phenoplot1,
          base_height = 4, base_width = 4)



#####*********************************************************************######
##### PAPER FIGS #####
#### Correlation plots vs 2D climate for paper ####

mmm %>%
  dplyr::filter(year==1) %>%
ggplot()+
  geom_jitter(aes(x=sites_bio18, y=sites_bio1, color=r_c13, shape=r_p_c13<0.05),
              size=5, height = 2, width = 10)+
  scale_colour_gradient2("Water Use Efficiency \nevolution (r)",
                        low = "#E3AF70",mid = "lightgrey",high = "#7B97CA",
                         midpoint = 0
                         )+
  scale_shape_manual("p<0.05",values = c("TRUE"=16,"FALSE"=1 ))+
  ylab("Annual temperature (C)")+
  xlab("Summer precipitation (mm)")+
  theme_classic()->
  pwue
  # theme_minimal()
pwue


pft<-
  mmm %>%
  dplyr::filter(year==1) %>%
  ggplot()+
  geom_jitter(aes(x=sites_bio18, y=sites_bio1, color=r_ft10,shape=r_p_ft10<0.05),
              size=5, height = 2, width = 10)+
  scale_colour_gradient2("Flowering time \nevolution (r)",
                         low = "#76AB33",mid = "lightgrey",high = "#798AE9",
                         midpoint = 0
  )+
  scale_shape_manual("p<0.05",values = c("TRUE"=16,"FALSE"=18 ))+

  ylab("Annual temperature (C)")+
  xlab("Summer precipitation (mm)")+
  theme_classic()
pft

pla<-
  mmm %>%
  dplyr::filter(year==1) %>%
  ggplot()+
  geom_jitter(aes(x=sites_bio18, y=sites_bio1, color=r_la,shape=r_p_la<0.05),
              size=5, height = 2, width = 10)+
  scale_colour_gradient2("Leaf area \nevolution (r)",
                         low = "black",mid = "lightgrey",high = "#76AB33",
                         midpoint = 0, position="bottom"
  )+
  scale_shape_manual("p<0.05",values = c("TRUE"=16,"FALSE"=18 ))+
  ylab("Annual temperature (C)")+
  xlab("Summer precipitation (mm)")+
  theme_classic()
pla

proot<-
  mmm %>%
  dplyr::filter(year==1) %>%
  ggplot()+
  geom_jitter(aes(x=sites_bio18, y=sites_bio1, color=r_roothorizontality,shape=r_p_roothorizontality<0.05),
              size=5, height = 2, width = 10)+
  scale_colour_gradient2("Root horizontality \nevolution (r)",
                         low = "gold3",mid = "lightgrey",high = "plum3",
                         midpoint = 0, position="bottom"
  )+
  scale_shape_manual("p<0.05",values = c("TRUE"=16,"FALSE"=18 ))+
  ylab("Annual temperature (C)")+
  xlab("Summer precipitation (mm)")+
  theme_classic()
proot

mmm %>%
  dplyr::filter(year==1) %>%
  ggplot()+
  geom_jitter(aes(x=sites_bio18, y=sites_bio1, color=r_rosettearea,shape=r_p_rosettearea<0.05),
              size=5, height = 2, width = 10)+
  scale_colour_gradient2("Leaf rosette area \nevolution (r)",
                         low = "#E3AF70",mid = "lightgrey",high = "#76AB33",
                         midpoint = 0, position="bottom"
  )+
  scale_shape_manual("p<0.05",values = c("TRUE"=16,"FALSE"=18 ))+
  ylab("Annual temperature (C)")+
  xlab("Summer precipitation (mm)")+
  theme_classic()->
  prosette
prosette


pdormancy<-
  mmm %>%
  dplyr::filter(year==1) %>%
  ggplot()+
  geom_jitter(aes(x=sites_bio18, y=sites_bio1,
                 color=r_dormancy,shape=r_p_dormancy<0.05),
              size=5, height = 2, width = 10)+
  scale_colour_gradient2("Seed dormancy \nevolution (r)",
                         low = "#BF812D",mid = "lightgrey",high = "black",
                         midpoint = 0, position="bottom"
  )+
  scale_shape_manual("p<0.05",values = c("TRUE"=16,"FALSE"=18 ))+
  ylab("Annual temperature (C)")+
  xlab("Summer precipitation (mm)")+
  theme_classic()
pdormancy


pstomata<-
  mmm %>%
  dplyr::filter(year==1) %>%
  ggplot()+
  geom_jitter(aes(x=sites_bio18, y=sites_bio1, color=r_stomata_density,
                 shape=r_p_stomata<0.1),
              size=5, height = 2, width = 10)+
  scale_colour_gradient2("Stomata density \nevolution (r)",
                         low = "#8C510A",mid = "lightgrey",high = "#006D2C",
                         midpoint = 0, position="bottom"
  )+
  scale_shape_manual("p<0.05",values = c("TRUE"=16,"FALSE"=18 ))+
  ylab("Annual temperature (C)")+
  xlab("Summer precipitation (mm)")+
  theme_classic()
pstomata


phenotypetrajectories<-
  plot_grid(
  plot_grid(pft,pdormancy, proot, ncol=1, align = "h",labels = c("A","C",'E')),
  plot_grid(pwue, pstomata,prosette, ncol=1,align = "h",labels = c("B","D",'F')),
  ncol=2
  )
phenotypetrajectories

save_plot("figs/fig-phenotypes-ecotype_frequencies-climate.pdf",phenotypetrajectories,
          base_height = 12, base_width = 10)
save_plot("figs/fig-phenotypes-ecotype_frequencies-climate.png",phenotypetrajectories,
          base_height = 12, base_width = 10)

phenotypetrajectories



#### Correlation vs climate  ####
#### Now of WUE with summer precipitation
wueplotgradient<-
  mmm %>%
  dplyr::filter(longitude > -15) %>%
  ggplot(.,
         aes(x=sites_bio18, y=r_c13, color=r_c13,size=year, shape=r_p_c13<0.05)
  )+
  geom_point()+
  # geom_point(size=5.2, color='white')+
  # geom_point(size=5.2, color='black')+
  # geom_point(size=5)+
  scale_colour_gradient2("WUE evolution (r)",
                         low = "#E3AF70",mid = "lightgrey",high = "#7B97CA",
                         midpoint = 0
  )+
  scale_size_continuous(breaks = c(1,2,3),range = c(3,6))+
  geom_hline(yintercept = 0,lty="dotted")+
  scale_shape_manual("p<0.05",values = c("TRUE"=16,"FALSE"=1 ))+
  stat_smooth(method='glm',color='grey', se=F, formula = y~poly(x,2),
              data = mmm %>% filter(longitude > -15, r_p_c13 < 0.05) # Filtering for p < 0.0
              )+
  # scale_x_log10()+
  ylab("WUE evolution (b)")+
  xlab("Summer precipitation (mm)")+
  theme_minimal()+
  guides(size = "none") # Removes the shape legend

wueplotgradient
save_plot("figs/fig-phenotypes-ecotype_frequencies-summerprecipitation-WUE-gradient.pdf",
          wueplotgradient,
          base_height = 5, base_width = 6)
save_plot("figs/fig-phenotypes-ecotype_frequencies-summerprecipitation-WUE-gradient.png",
          wueplotgradient,
          base_height = 5, base_width = 6)




floweringplotgradient<-
  mmm %>%
  dplyr::filter(longitude > -15) %>%
  # dplyr::filter(year == 3) %>%
  ggplot(.,
         aes(x=sites_bio1, y=r_ft10, color=r_ft10, size=year,shape=r_p_ft10<0.05)
  )+
  # geom_point(size=5.2, color='white')+
  # geom_point(size=5.2, color='black')+
  geom_point()+
  scale_size_continuous(breaks = c(1,2,3),range = c(3,6))+
  scale_colour_gradient2("Flowering time evolution (r)",
                         low = "#76AB33",mid = "lightgrey",high = "#798AE9",
                         midpoint = 0
  )+
  geom_hline(yintercept = 0,lty="dotted")+
  scale_shape_manual("p<0.05",values = c("TRUE"=16,"FALSE"=1))+
  stat_smooth(method='glm',color='grey', se=F, formula = y~poly(x,2),
              data = mmm %>% filter(longitude > -15, r_p_ft10 < 0.05) # Filtering for p < 0.05
  )+
  # scale_x_log10()+
  ylab("Flowering evolution (r)")+
  xlab("Annual temperature (C)")+
  theme_minimal()+
  guides(size = "none") # Removes the shape legend

floweringplotgradient
save_plot("figs/fig-phenotypes-ecotype_frequencies-annualtempearture-flowering-gradient.pdf",
          floweringplotgradient,
          base_height = 5, base_width = 6)
save_plot("figs/fig-phenotypes-ecotype_frequencies-annualtempearture-flowering-gradient.png",
          floweringplotgradient,
          base_height = 5, base_width = 6)


dormancyplotgradient<-
  mmm %>%
  dplyr::filter(longitude > -15) %>%
  # dplyr::filter(year == 3) %>%
  ggplot(.,
         aes(x=sites_bio1, y=r_dormancy, color=r_dormancy, size=year,shape=r_p_dormancy<0.05)
  )+
  # geom_point(size=5.2, color='white')+
  # geom_point(size=5.2, color='black')+
  geom_point()+
  scale_size_continuous(breaks = c(1,2,3),range = c(3,6))+
  scale_colour_gradient2("dormancy time evolution (r)",
                         low = "#BF812D",mid = "lightgrey",high = "black",
                         midpoint = 0
  )+
  geom_hline(yintercept = 0,lty="dotted")+
  scale_shape_manual("p<0.05",values = c("TRUE"=16,"FALSE"=1))+
  stat_smooth(method='glm',color='grey', se=F, formula = y~poly(x,2),
              data = mmm %>% filter(longitude > -15, r_p_dormancy < 0.05) # Filtering for p < 0.05
  )+
  # scale_x_log10()+
  ylab("Dormancy evolution (r)")+
  xlab("Annual temperature (C)")+
  theme_minimal()+
  guides(size = "none") # Removes the shape legend

dormancyplotgradient
save_plot("figs/fig-phenotypes-ecotype_frequencies-annualtempearture-dormancy-gradient.pdf",
          dormancyplotgradient,
          base_height = 5, base_width = 6)
save_plot("figs/fig-phenotypes-ecotype_frequencies-annualtempearture-dormancy-gradient.png",
          dormancyplotgradient,
          base_height = 5, base_width = 6)




phenotypecorrelationclimates<-
  plot_grid(floweringplotgradient,dormancyplotgradient,wueplotgradient ,ncol=1,align = "hv",labels = c("A","B",'C'))
phenotypecorrelationclimates

save_plot("figs/fig-phenotypes-ecotype_frequencies_correlation_traits_vs_climates_grid.png",
          phenotypecorrelationclimates,
          base_height = 12, base_width = 7)
save_plot("figs/fig-phenotypes-ecotype_frequencies_correlation_traits_vs_climates_grid.pdf",
          phenotypecorrelationclimates,
          base_height = 12, base_width = 7)



# #### Correlation plot vs latitude
#
# floweringlatitude<-
# ggplot(mmm)+
#   geom_point(aes(y=latitude, x=sites_bio1),color="white",
#              size=5.2)+
#   geom_point(aes(y=latitude, x=sites_bio1),color="black",
#              size=5.1)+
#   geom_point(aes(y=latitude, x=sites_bio1, color=r_ft10),
#               size=5)+
#   scale_colour_gradient2("FT10 (r)",
#                          # low = "#E69D4E",mid = "lightgrey",high = "#728DC4",
#                          low = "#76AB33",mid = "lightgrey",high = "#798AE9",
#                          # low = "#76AB33",mid = "grey",high = "#798AE9",
#                          # low = "#C97B3A",mid = "lightgrey",high = "#2E5EAB",
#                          midpoint = 0
#   )+
#   ylab("Latitude (N)")+
#   xlab("Annual temperature (C)")+
#   # theme_classic()
#   theme_minimal()
# floweringlatitude
#
# floweringlatitude
# save_plot("figs/fig-phenotypes-ecotype_frequencies-climate-flowering-and-latitude.pdf",floweringlatitude,
#           base_height = 4, base_width = 5)
# save_plot("figs/fig-phenotypes-ecotype_frequencies-climate-flowering-and-latitude.png",floweringlatitude,
#           base_height = 4, base_width = 5)
#
#
# floweringcorrelation<-
#   ggplot(mmm)+
#   geom_point(aes(y=r_ft10, x=sites_bio1),color="white",
#              size=5.2)+
#   geom_point(aes(y=r_ft10, x=sites_bio1),color="black",
#              size=5.1)+
#   geom_point(aes(y=r_ft10, x=sites_bio1, color=r_ft10),
#              size=5)+
#   scale_colour_gradient2("FT10 (r)",
#                          # low = "#E69D4E",mid = "lightgrey",high = "#728DC4",
#                          low = "#76AB33",mid = "lightgrey",high = "#798AE9",
#                          # low = "#76AB33",mid = "grey",high = "#798AE9",
#                          # low = "#C97B3A",mid = "lightgrey",high = "#2E5EAB",
#                          midpoint = 0
#   )+
#   ylab("Correlation ecotype change with flowering (r)")+
#   xlab("Annual temperature (C)")+
#   geom_hline(yintercept = 0,color='grey',lty='dotted')+
#   # theme_classic()
#   theme_minimal()
# floweringcorrelation
#
#
# save_plot("figs/fig-phenotypes-ecotype_frequencies-climate-flowering-correlation.pdf",floweringcorrelation,
#           base_height = 4, base_width = 5)
# save_plot("figs/fig-phenotypes-ecotype_frequencies-climate-flowering-correlation.png",floweringcorrelation,
#           base_height = 4, base_width = 5)
#
#
# floweringcorrelationsummmer<-
#   ggplot(mmm)+
#   geom_point(aes(y=r_ft10, x=sites_bio18),color="white",
#              size=5.2)+
#   geom_point(aes(y=r_ft10, x=sites_bio18),color="black",
#              size=5.1)+
#   geom_point(aes(y=r_ft10, x=sites_bio18, color=r_ft10),
#              size=5)+
#   scale_colour_gradient2("FT10 (r)",
#                          # low = "#E69D4E",mid = "lightgrey",high = "#728DC4",
#                          low = "#76AB33",mid = "lightgrey",high = "#798AE9",
#                          # low = "#76AB33",mid = "grey",high = "#798AE9",
#                          # low = "#C97B3A",mid = "lightgrey",high = "#2E5EAB",
#                          midpoint = 0
#   )+
#   ylab("Correlation accession change with flowering (r)")+
#   xlab("Summer precipitation (mm)")+
#   geom_hline(yintercept = 0,color='grey',lty='dotted')+
#   # theme_classic()
#   theme_minimal()
# floweringcorrelationsummmer
#
#
# save_plot("figs/fig-phenotypes-ecotype_frequencies-climate-summerprecipitation-flowering-correlation.pdf",
#           floweringcorrelationsummmer,
#           base_height = 4, base_width = 5)
# save_plot("figs/fig-phenotypes-ecotype_frequencies-climate-summerprecipitation-flowering-correlation.png",
#           floweringcorrelationsummmer,
#           base_height = 4, base_width = 5)

# # Flowerig time against summer precipitation
# mmm %>%
#   dplyr::filter(longitude > -15) %>%
#   ggplot(.,
#          aes(x=sites_bio18, y=r_ft16, color=r_c13)
#   )+
#   geom_point(size=5.2, color='white')+
#   geom_point(size=5.2, color='black')+
#   geom_point(size=5)+
#   scale_colour_gradient2("WUE (r)",
#                          low = "#E3AF70",mid = "lightgrey",high = "#7B97CA",
#                          midpoint = 0
#   )+
#   geom_hline(yintercept = 0,lty="dotted")+
#   stat_smooth(method='glm',color='grey', se=F, formula = y~poly(x,2))+
#   # scale_x_log10()+
#   ylab("FT (r)")+
#   xlab("Summer precipitation (mm)")+
#   theme_minimal()
#
# # Flowerig time against summer precipitation
# mmm %>%
#   dplyr::filter(longitude > -15) %>%
#   ggplot(.,
#          aes(x=sites_bio18, y=r_ft16, color=b_c13, shape=p_c13<0.05)
#   )+
#   geom_point(size=5.2, color='white')+
#   geom_point(size=5.2, color='black')+
#   geom_point(size=5)+
#   scale_colour_gradient2("WUE (r)",
#                          low = "#E3AF70",mid = "lightgrey",high = "#7B97CA",
#                          midpoint = 0
#   )+
#   geom_hline(yintercept = 0,lty="dotted")+
#   stat_smooth(method='glm',color='grey', se=F, formula = y~poly(x,2))+
#   # scale_x_log10()+
#   ylab("FT (r)")+
#   xlab("Summer precipitation (mm)")+
#   theme_minimal()


##### Gradients of selection plots

# wueselection<-
#   mmm %>%
#   dplyr::filter(longitude > -15) %>%
#   ggplot(.,
#          aes(x=sites_bio18, y=b_c13, color=b_c13, shape=p_c13<0.05)
#   )+
#   geom_point(size=5.2, color='white')+
#   geom_point(size=5.2, color='black')+
#   geom_point(size=5)+
#   scale_shape_manual(values = c("TRUE"=16,"FALSE"=1))+
#   scale_colour_gradient2("WUE (b)",
#                          low = "#E3AF70",mid = "lightgrey",high = "#7B97CA",
#                          midpoint = 0
#   )+
#   geom_hline(yintercept = 0,lty="dotted")+
#   stat_smooth(
#     aes(x=sites_bio18, y=b_c13),
#     method='glm',color='grey', se=F, formula = y~poly(x,2))+
#   # scale_x_log10(limits=c(1,1000))+
#   ylab("WUE marginal correlation (r)")+
#   xlab("Summer precipitation (mm)")+
#   theme_minimal()
# wueselection


#### Phenotype trajectory plot across cliamte ####
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

phenotype_trajectory_from_founder<-
  sitephenotypes %>%
    dplyr::filter(longitude > -15) %>%
    ggplot(.)+
    # geom_point(data=pheno,aes(y=FT16, x=Delta_13C),shape=3)+ # the ecotypes
    geom_segment(aes(yend=FT_mean, xend=WUE_mean, y=founderphenotypes$FT_mean, x=founderphenotypes$WUE_mean,color=sites_bio1),
                 arrow = arrow(angle = 30, length = unit(0.1, "inches"),
                    ends = "last", type = "open")
    )+
    geom_point(aes(y=FT_mean, x=WUE_mean,  color=sites_bio1), size=3)+ # the sites
    geom_point(data=founderphenotypes,aes(y=FT_mean, x=WUE_mean),color="black", size=3)+
    xlab("WUE (d13C)")+ylab("Flowering time (days)")+
    scale_color_gradientn("Temperature (C)",colors = rev(brewer.pal(5,"RdBu")))+
    theme_minimal()

phenotype_trajectory_from_founder

save_plot("figs/fig-phenotypes-trajectories-from-founder-FT16-WUE.pdf",
          phenotype_trajectory_from_founder,
          base_height = 4, base_width = 5)
save_plot("figs/fig-phenotypes-trajectories-from-founder-FT16-WUE.png",
          phenotype_trajectory_from_founder,
          base_height = 4, base_width = 5)

