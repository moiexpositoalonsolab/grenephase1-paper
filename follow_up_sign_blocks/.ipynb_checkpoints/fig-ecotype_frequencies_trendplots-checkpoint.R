################################################################################
### Goal
### Do a plot like in the E coli trend plots, Muller plots, but with ecotypes
### 

################################################################################
# 
# #### If in google drive
# from google.colab import drive
# drive.mount('/content/drive')
# %load_ext rpy2.ipython


################################################################################

# library(grene)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

################################################################################

# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE

setwd(myfolder)

################################################################################

    # "#4C2737",
    # "#3D323C",
    # "#55586D",
aracolors=c(
    "#8088A6",
    "#B4ABBE",
    "#A084BD",
    "#B376AF",
    "#C89EB8",
    "#B44472",
    "#BC686F",
    "#A36378",
    "#F4B5C4",
    "#F19BA0",
    "#A35B5F",
    "#AD7F7B",
    "#7F5541",
    "#855E46",
    "#E87858",
    "#F2A445",
    "#E89D4C",
    "#E0C5B3",
    "#A69868",
    "#AFAB61",
    "#B4C838",
    "#B4D77C"
)
aracolors<-rev(aracolors)

################################################################################
# Load datasets
# These were the old ones
# load("~/Google Drive/RESEARCH/MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/figures-raw/plot_ecotype_trend/ecotype_cline_bio1.RData")
# setwd("~/Google Drive/RESEARCH/MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/figures-raw/plot_ecotype_trend/")

# Frozen datasets
load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate.rda")
load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate.rda")
p0=read.table("data/founder_ecotype_frequency.txt",header = F) %>% 
  dplyr::rename(id=V1,freq=V2)
head(p0)


# My main dataset
ecofreq<-ecotype_allyears_frequencies_long_raw_climate

# make the first generation the same across sites
extractgeneration0<- ecofreq %>%
  dplyr::filter(year==1) %>% # just get the generation 1 data for everyone, which is the most complete
  dplyr::mutate(year=0) %>%
  dplyr::mutate(freq=startfreq) 

ecofreq<- rbind(ecofreq, extractgeneration0)

ecofreq$id<-as.factor(ecofreq$id)
ecofreq$year<-as.numeric(ecofreq$year)

################################################################################
# Test it worked
ex<-ecofreq %>% 
  dplyr::filter(site==4, rep==3)

fig_ecotype_mullerplot_site4_rep3<-
ggplot(ex, 
       aes(x = year, y = freq, fill = as.factor(id), group = as.factor(id))) + 
  geom_area()+
  scale_fill_manual(values=colorRampPalette(aracolors)(length(unique(ex$id))))+
  xlab("Years")+ylab("Ecotype frequency")+
  theme_minimal()+
  theme(legend.position = "none")



save_plot("figs/fig_ecotype_mullerplot_site4_rep3.pdf",fig_ecotype_mullerplot_site4_rep3,
          base_width = 6,base_height = 5 )
save_plot("figs/fig_ecotype_mullerplot_site4_rep3.png",fig_ecotype_mullerplot_site4_rep3,
          base_width = 6,base_height = 5 )

# Run for all sites and rep combinations
uniquesites<-unique(ecofreq$site)
for(mysite in uniquesites){
tmp<-ecofreq %>% 
  dplyr::filter(site==mysite)

uniquereps<-unique(tmp$rep)
for(myrep in uniquereps){
ex<-tmp %>% 
  dplyr::filter(rep==myrep)

fig_ecotype_mullerplot_site4_rep3<-
  ggplot(ex, 
         aes(x = year, y = freq, fill = as.factor(id), group = as.factor(id))) + 
  geom_area()+
  scale_x_continuous(breaks = c(0:4))+
  scale_fill_manual(values=colorRampPalette(aracolors)(length(unique(ex$id))))+
  xlab("Years")+ylab("Ecotype frequency")+
  theme_minimal()+
  theme(legend.position = "none")
fig_ecotype_mullerplot_site4_rep3

save_plot(paste0("figs/fig_ecotype_mullerplots/fig_ecotype_mullerplot_site",mysite,"_rep",myrep,".pdf"),fig_ecotype_mullerplot_site4_rep3,
          base_width = 6,base_height = 5 )
save_plot(paste0("figs/fig_ecotype_mullerplots/fig_ecotype_mullerplot_site",mysite,"_rep",myrep,".png"),fig_ecotype_mullerplot_site4_rep3,
          base_width = 6,base_height = 5 )

}
}

################# Trying the line plot
mysite=4
myrep=3

fig_ecotype_lineplot_site4_rep3<-
  
ggplot(ex, 
       aes(x = year, y = freq+0.001, fill = as.factor(id),color = as.factor(id), group = as.factor(id))) + 
  # geom_area()+
  geom_line()+
  scale_x_continuous(breaks = c(0:4))+
  scale_y_log10()+
  scale_fill_manual(values=colorRampPalette(aracolors)(length(unique(ex$id))))+
  scale_color_manual(values=colorRampPalette(aracolors)(length(unique(ex$id))))+
  xlab("Years")+ylab("Ecotype frequency")+
  theme_minimal()+
  theme(legend.position = "none")

save_plot(paste0("figs//fig_ecotype_lineplot_site",mysite,"_rep",myrep,".pdf"),
          fig_ecotype_lineplot_site4_rep3,
          base_width = 6,base_height = 5 )
save_plot(paste0("figs//fig_ecotype_lineplot_site",mysite,"_rep",myrep,".png"),
          fig_ecotype_lineplot_site4_rep3,
          base_width = 6,base_height = 5 )


################################################################################
### Now averages per site

ex<-ecofreq %>%
  # dplyr::filter(site==4) %>%
  group_by(id,site,year) %>%
  summarise(freq=mean(freq),
            year=head(year,1)
  ) %>%
  ungroup %>%
  group_by(year) %>%
  mutate(freqsum=sum(freq)) %>%
  ungroup %>%
  mutate(freqscale=freq/freqsum)

ggplot(ex, 
       aes(x = year, y = freq, fill = as.factor(id), group = as.factor(id))) + 
  geom_area()+
  scale_fill_manual(values=colorRampPalette(aracolors)(length(unique(ex$id))))+
  theme_minimal()+
  theme(legend.position = "none")



ex<-ecofreq %>%
  # dplyr::filter(site==4) %>%
  group_by(id,site,year,sites_bio1, ecotypes_bio1) %>%
  summarise(freq=mean(freq),
            year=head(year,1)
  ) %>%
  ungroup %>%
  group_by(year,site,sites_bio1) %>%
  mutate(freqsum=sum(freq)) %>%
  ungroup %>%
  mutate(freqscale=freq/freqsum)


# Save plot
mullerallsites<-
ggplot(ex, 
       aes(x = year, y = freq, fill = as.factor(id), group = as.factor(id))) + 
  geom_area()+
  scale_fill_manual(values=colorRampPalette(aracolors)(length(unique(ex$id))))+
  theme_minimal()+
  xlab("Years")+ylab("Ecotype frequency")+
  theme(legend.position = "none")+
  facet_grid(~site)
# facet_wrap(~site,ncol = 1, nrow=31)

save_plot(paste0("figs//fig_ecotype_mullerplot_allsites.pdf"),mullerallsites,
          base_width = 16,base_height = 5 )
save_plot(paste0("figs//fig_ecotype_mullerplot_allsites.png"),mullerallsites,
          base_width = 16,base_height = 5 )


##########
# oder the factor based on bio1
library(forcats)
sites_simple_names<-read.csv("grene/data/sites_simple_names.csv",header = T)
simplenames<-paste0(sites_simple_names$city,", ", sites_simple_names$country)
names(simplenames)<-sites_simple_names$site

# Recode the sites
ex$site <- recode(as.character(ex$site), !!!simplenames)

# Recoder sites and ecotypes based on bio1
ex$site <- fct_reorder(ex$site, ex$sites_bio1)
ex$id <- fct_reorder(ex$id, ex$ecotypes_bio1)

mullerallsitesordered<-
ggplot(ex, 
       aes(x = year, y = freq, fill = as.factor(id), group = as.factor(id))) + 
  geom_area()+
  scale_fill_manual(values=colorRampPalette(aracolors)(length(unique(ex$id))))+
  theme_minimal()+
  xlab("Years")+ylab("Ecotype frequency")+
  theme(legend.position = "none",
        strip.text = element_text(angle = 90, hjust = 0))+
  facet_grid(~site)
  # facet_wrap(~site,ncol = 1, nrow=31)
  
save_plot(paste0("figs//fig_ecotype_mullerplot_allsites-orderedbio1.pdf"),mullerallsitesordered,
          base_width = 20,base_height = 5 )
save_plot(paste0("figs//fig_ecotype_mullerplot_allsites-orderedbio.png"),mullerallsitesordered,
          base_width = 20,base_height = 5 )

################################################################################
# #### Visualize classic evolutionary trends
# 
# # Example with site 4
# extractgeneration0<- ecofreq %>%
#     # dplyr::filter(site=="4",rep==1 ) %>%
#     dplyr::filter(year==1) %>%
#     dplyr::mutate(year=0) %>%
#     dplyr::mutate(frequency=p0) %>%
#     dplyr::select(ecotypeid, frequency, generation,p0) %>%
#     dplyr::select(-p0)
# ex<-
#     ecofreq %>%
#     dplyr::filter(site=="4",plot==1 ) %>%
#     dplyr::select(ecotypeid, frequency, generation,p0) %>%
#     dplyr::select(-p0) %>%
#     rbind(.,
#             extractgeneration0
#           )
# ex %>% head
# ex %>% tail
# 
# # Plot from generation 0
# # ex$ecotypeid <- as.factor(ex$ecotypeid)  # Convert ecotypeid to factor
# 
# 
# 
# ggplot(ex, aes(x = generation, y = frequency, fill = ecotypeid)) +
#     geom_area()+
#     scale_fill_manual(values=colorRampPalette(aracolors)(length(unique(ex$ecotypeid))))+
#     theme(legend.position = "none")

################################################################################
#### Visualize classic evolutionary trends - SUM
# mysite=6
# mysite=52
# 
# # Example with site 4
# extractgeneration0<- ecofreq %>%
#     dplyr::filter(site==mysite,plot==1 ) %>%
#     dplyr::filter(generation==1) %>%
#     dplyr::mutate(generation=0) %>%
#     # dplyr::mutate(frequency=p0) %>%
#     dplyr::mutate(frequency=p0 * 12) %>%
#     dplyr::select(ecotypeid, frequency, generation,p0) %>%
#     dplyr::select(-p0)
# ex<-
#     ecofreq %>%
#     dplyr::filter(site==mysite) %>%
#     dplyr::select(ecotypeid, frequency, generation) %>%
#     group_by(ecotypeid, generation) %>%
#     dplyr::summarise(frequency=sum(frequency),
#                      generation=head(generation,1),
#                      ) %>%
#     rbind(.,
#           extractgeneration0
#     )
# ex %>% head
# ex %>% tail
# 
# # Plot from generation 0
# ex$ecotypeid <- as.factor(ex$ecotypeid)  # Convert ecotypeid to factor
# 
# ggplot(ex, aes(x = generation, y = frequency, fill = ecotypeid)) +
#     geom_area()+
#     scale_fill_manual(values=colorRampPalette(aracolors)(length(unique(ex$ecotypeid))))+
#     theme(legend.position = "none")
# 
# 
# #### Visualize classic evolutionary trends - MEAN
# mysite=52
# mysite=6
# mysite=43
# mysite=5
# mysite=4
# 
# bacteriaplot<-function(mysite=4){
#     # Example with site 4
#     extractgeneration0<- ecofreq %>%
#         dplyr::filter(site==mysite,plot==1 ) %>%
#         dplyr::filter(generation==1) %>%
#         dplyr::mutate(generation=0) %>%
#         # dplyr::mutate(frequency=p0) %>%
#         dplyr::mutate(frequency=p0 ) %>%
#         dplyr::select(ecotypeid, frequency, generation,p0) %>%
#         dplyr::select(-p0)
#     ex<-
#         ecofreq %>%
#         dplyr::filter(site==mysite) %>%
#         dplyr::select(ecotypeid, frequency, generation) %>%
#         group_by(ecotypeid, generation) %>%
#         dplyr::summarise(frequency=mean(frequency),
#                          generation=head(generation,1),
#         ) %>%
#         rbind(.,
#               extractgeneration0
#         )
#     ex %>% head
#     ex %>% tail
# 
#     # Plot from generation 0
#     ex$ecotypeid <- as.factor(ex$ecotypeid)  # Convert ecotypeid to factor
# 
#     ggplot(ex, aes(x = generation, y = frequency, fill = ecotypeid)) +
#         geom_area()+
#         scale_fill_manual(values=
#                              (colorRampPalette(rev(aracolors))(length(unique(ex$ecotypeid))))
# 
#                           )+
#         # scale_fill_manual(values=colorRampPalette(brewer.pal(5,"Spectral"))(length(unique(ex$ecotypeid))))+
#         theme(legend.position = "none",
#               axis.title.x = element_blank(),
#               axis.title.x.bottom = element_blank(),
#               axis.line.x  = element_blank(),
#               axis.text.x =  element_blank()
#               )
# }
# 
# 
# sites=unique(ecofreq$site)
# allplots<-
#     lapply((sites),
#            function(i) bacteriaplot(i)
#            )
# 
# plot_grid(plotlist = allplots,ncol=1,
#           hjust = -0.5,
#           vjust = 1
#           )
# 



# 
# ################################################################################
# 
# # try wrap
# 
# mysite=c(4,5)
# mysit=sites
# extractgeneration0<- ecofreq %>%
#     dplyr::filter(plot==1 ) %>%
#     # dplyr::filter(site %in% mysite) %>%
#     dplyr::filter(generation==1) %>%
#     dplyr::mutate(generation=0) %>%
#     dplyr::mutate(frequency=p0 ) %>%
#     dplyr::select(ecotypeid, frequency, generation, site)
# ex<-
#     ecofreq %>%
#     # dplyr::filter(site %in% mysite) %>%
#     dplyr::select(ecotypeid, frequency, generation, site) %>%
#     group_by(ecotypeid, generation, site) %>%
#     dplyr::summarise(frequency=mean(frequency),
#                      generation=head(generation,1),
#     ) %>%
#     rbind(.,
#           extractgeneration0
#     )
# ex %>% head
# ex %>% tail
# 
# # Plot from generation 0
# ex$ecotypeid <- as.factor(ex$ecotypeid)  # Convert ecotypeid to factor
# 
# # Plot
# allplotsbacteria<-
# ggplot(ex, aes(x = generation, y = frequency, fill = ecotypeid)) +
#     geom_area()+
#     scale_fill_manual(values=
#                           (colorRampPalette(rev(aracolors))(length(unique(ex$ecotypeid))))
# 
#     )+
#     # scale_fill_manual(values=colorRampPalette(brewer.pal(5,"Spectral"))(length(unique(ex$ecotypeid))))+
#     facet_wrap(~site, ncol= 1 )+
#     theme(legend.position = "none",
#           axis.title.x = element_blank(),
#           axis.title.x.bottom = element_blank(),
#           axis.line.x  = element_blank(),
#           axis.text.x =  element_blank()
#     )
# 
# 
# save_plot("fig-ecotype_frequencies_trendplots.pdf",allplotsbacteria,
#           base_width = 4,base_height = 20 )
# 
# 
# 
# ################################################################################
# 
# 
# # try wrap
# 
# mysite=c(4,55)
# mysit=sites
# extractgeneration0<- ecofreq %>%
#     dplyr::filter(plot==1 ) %>%
#     dplyr::filter(site %in% mysite) %>%
#     dplyr::filter(generation==1) %>%
#     dplyr::mutate(generation=0) %>%
#     dplyr::mutate(frequency=p0 ) %>%
#     dplyr::select(ecotypeid, frequency, generation, site)
# ex<-
#     ecofreq %>%
#     dplyr::filter(site %in% mysite) %>%
#     dplyr::select(ecotypeid, frequency, generation, site) %>%
#     group_by(ecotypeid, generation, site) %>%
#     dplyr::summarise(frequency=mean(frequency),
#                      generation=head(generation,1),
#     ) %>%
#     rbind(.,
#           extractgeneration0
#     )
# ex %>% head
# ex %>% tail
# 
# # Plot from generation 0
# ex$ecotypeid <- as.factor(ex$ecotypeid)  # Convert ecotypeid to factor
# 
# # Plot
# allplotsbacteria<-
#     ggplot(ex, aes(x = generation, y = frequency, fill = ecotypeid)) +
#     geom_area()+
#     scale_fill_manual(values=
#                           (colorRampPalette(rev(aracolors))(length(unique(ex$ecotypeid))))
# 
#     )+
#     # scale_fill_manual(values=colorRampPalette(brewer.pal(5,"Spectral"))(length(unique(ex$ecotypeid))))+
#     facet_wrap(~site, ncol= 1 )+
#     theme(legend.position = "none",
#           axis.title.x = element_blank(),
#           axis.title.x.bottom = element_blank(),
#           axis.line.x  = element_blank(),
#           axis.text.x =  element_blank()
#     )
# 
# allplotsbacteria
# save_plot("fig-ecotype_frequencies_trendplots-examplesites.pdf",allplotsbacteria,
#           base_width = 4,base_height = 5 )
# 

