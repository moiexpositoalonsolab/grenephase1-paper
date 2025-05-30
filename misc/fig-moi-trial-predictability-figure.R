################################################################################
### Goal
### Visualize predictability
### combine with survival in different ways

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



library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(latex2exp)
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)
source("function-colors-grenenet.R")

# prefix <- "FUTURE_PREDICTION_"
# PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
# PATH <-  "~/Shareddrives/MOI-LAB/PROJECTS/grenenet//GRENE-net_PHASE1/grenephase1-analyses/"
# setwd(PATH)

#####*********************************************************************######
#### Load climate info #####

## MOD 1
# climate <- read.delim("data-intermediate/Leave_one_out_prediction/climate_distance/Genomic_offset_site_climate_distance.txt")
climate <- read.delim("data-intermediate/Leave_one_out_prediction/climate_distance/Genomic_offset_climate_distance_loo_prediction_summary.txt")

# MOD 2
# Go of p1/p0
go <- read.delim("data-intermediate/Leave_one_out_prediction/binomial_regression/genomic_offset_firstgen_binom_reg_p1_p0_results.csv",sep=",") %>%
# Go of log p/p0
# go <- read.delim("data-intermediate/Leave_one_out_prediction/binomial_regression/genomic_offset_firstgen_binom_reg_log_p1_p0_results.csv",sep=",") %>%
# Go of delta p
# go <- read.delim("data-intermediate/Leave_one_out_prediction/binomial_regression/genomic_offset_firstgen_binom_reg_results.csv",sep=",") %>%
  mutate(source = "genomic offset",
         r2 = r_squared,
         pearson_r = pearsonr,
         sp_r = sp_correlation) %>%
  dplyr::select(site,plot,pearson_r,sp_r,r2,source)
replicates <- read.delim("data-intermediate/Leave_one_out_prediction/replicates/Genomic_offset_plot_replication_prediction_summary.txt")

# MOD 3 and 4
stabilizing_era5 <- read.delim("data-intermediate/Leave_one_out_prediction/stabilizing_selection_era5/Genomic_offset_stabilizing_selection_loo_era5_prediction_summary.txt")
stabilizing_prs <- read.delim("data-intermediate/Leave_one_out_prediction/stabilizing_selection_prs/Genomic_offset_stabilizing_selection_loo_p1_p0_prs_prediction_summary.txt")


# Baseline parallelism
parallelism <- read.delim("data-intermediate/generation_1_parallelism.txt")

# Climate data
load("grene/data/locations_data.rda")
locations_data %>% head
load("grene/data/worldclim_sitesdata.rda")
worldclim_sitesdata<-
  worldclim_sitesdata %>%
  dplyr::select(site,starts_with("bio")) %>%
  merge(.,locations_data,by="site")

#####*********************************************************************######
##### Combine datasets with climate #####
df<-rbind(
  # climate,
  go,
  replicates,
  stabilizing_era5,
  stabilizing_prs
  ) %>%
    data.frame
# Add an order
ordersource<-
  c(
    "plot_mean"=5,
    "stabilizing_prs"=4,
    "stabilizing_era5"=3,
    "genomic offset"=2,
    "climate_distance"=1
    )
df$sort<-ordersource[df$source]
# Add the worldclim data
df <-
  merge(df, worldclim_sitesdata,by="site")

#####*********************************************************************######
##### Plot general predictability #####
predictability_plot_general<-
ggplot(df)+
  # dots spearman
  geom_jitter(
     aes(x=sp_r,
         y=sort,
         color=sort
         ),
     height=0.1,width=0,
     alpha=0.5, size=0.25)+
  geom_boxplot(
    aes(x=sp_r,
        y=sort+0.25,
        group=sort,
        color=sort
    ),
      alpha=0.5,
    width=0.1,
    outliers = F)+
  # dots r2
  geom_jitter(
    aes(x=r2+1,
        y=sort,
        color=sort
        ),
    height=0.1,width = 0,
    alpha=0.5, size=0.25)+
  geom_boxplot(fun = mean,
               aes(x=r2+1,
                   y=sort+0.25,
                   group=sort,
                   color=sort
               ),
               alpha=0.5,
               width=0.1,
               outliers = F)+
  scale_color_gradientn(colours = rev(eggcolors))+
  xlim(c(-0.55, 2))+
  scale_x_continuous(breaks = round(seq(-0.6,2, by=0.2)*100)/100,
                     labels = c(round(seq(-0.6,0.8, by=0.2)*100)/100, round(seq(0,1, by=0.2)*100)/100)
                     )+
  theme_minimal()+
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())+
  theme(panel.grid.major.y = element_blank())+
  theme(panel.grid.minor.y = element_blank())+
  ylab("Spearman's r.            Variance explained (r2)")+
  ylab("")
predictability_plot_general

ggsave("figs/fig-Moitrial-predictability-LOO.pdf",
          predictability_plot_general,
          height = 2, width = 6, units="in")
ggsave("figs/fig-Moitrial-predictability-LOO.png",
          predictability_plot_general,
       height = 2, width = 6, units="in")

##### Plot each of the sites

plot_predictability_bysite<-
df %>%
  # dplyr::filter(source=="stabilizing_era5") %>%
ggplot(.)+
  ### Fit regressions
  geom_point(
    aes(
        y=sp_r,
        x=latitude,
        color=sort
    ),
    height=0.1,width=0,
    alpha=0.5, size=0.25)+
  geom_point(
    aes(
      y=r2+1,
      x=latitude,
      color=sort
    ),
    height=0.1,width=0,
    alpha=0.5, size=0.25)+
  # geom_boxplot(
  #              aes(
  #                  x=sp_r,
  #                  y=bio1,
  #                  group=paste(site, source),
  #                  color=sort
  #              ),
  #              alpha=0.5,
  #              width=0.3,
  #              outliers = F)+
  ### Fit regressions
  stat_smooth(
      aes(
        y=sp_r,
        x=(latitude),
        group=sort,
        color=sort
      ), method="glm", formula = y~poly(x,2))+
  ### Fit regressions
  stat_smooth(
    aes(
      y=r2+1,
      x=(latitude),
      group=sort,
      color=sort
    ), method="glm", formula = y~poly(x,2))+
  ### Themes
  coord_flip()+
  theme_minimal()+
  theme(legend.position = "none")+
  # theme(axis.text.y = element_blank())+
  # theme(panel.grid.major.y = element_blank())+
  # theme(panel.grid.minor.y = element_blank())+
  scale_color_gradientn(colours = rev(eggcolors))+
  ylim(c(-0.55, 2))+
  scale_y_continuous(breaks = round(seq(-0.6,2, by=0.2)*100)/100,
                     labels = c(round(seq(-0.6,0.8, by=0.2)*100)/100, round(seq(0,1, by=0.2)*100)/100)
  )+
  ylab("Spearman's r.            Variance explained (r2)")+
  xlab("Latitude")
plot_predictability_bysite

ggsave("figs/fig-Moitrial-predictability-LOO-bysite-latitude.pdf",
       plot_predictability_bysite,
       height = 5, width = 6, units="in")
ggsave("figs/fig-Moitrial-predictability-LOO-bysite-latitude.png",
       plot_predictability_bysite,
       height = 5, width = 6, units="in")


#####*********************************************************************######
#### Predictability above replicates
dfm<-
  df %>%
  dplyr::mutate(sp_r = ifelse(sp_r>0, sp_r, 0.001)) %>%
  dplyr::mutate(r2 = ifelse(r2>0, r2, 0.001)) %>%
  dplyr::group_by(site, plot) %>%
  dplyr::mutate(sp_r_reference=sp_r[source == "plot_mean"]) %>%
  dplyr::mutate(r2_reference=r2[source == "plot_mean"]) %>%
  dplyr::mutate(sprelative =  (sp_r-sp_r_reference)  ,
                spratio =  (sp_r/sp_r_reference)  ,
                r2ratio =  (r2/r2_reference)
                ) %>%
  dplyr::ungroup()

dfm$spratio
dfm$r2ratio %>% hist
dfm$r2ratio %>% summary
table(dfm$r2ratio<0)
dfm$r2ratio[dfm$r2ratio<0]

ggplot(dfm)+
# dots spearman
geom_jitter(
  aes(x=sprelative,
      y=sort,
      color=sort
  ),
  height=0.1,width=0,
  alpha=0.5, size=0.25)+
  theme_minimal()



gettest<-function(x=x){
  tmp<-t.test(x)
  tmp$estimate
}
gettestp<-function(x=x){
  tmp<-t.test(x)
  tmp$p.value
}
difwithtop<-
dfm %>%
  dplyr::filter(source=="stabilizing_era5") %>%
  group_by(site) %>%
  summarise(tdif=gettest(sprelative),
            tp=gettestp(sprelative))


(difwithtop$tp < 0.05/30) %>% table
summary(difwithtop$tdif)


#####*********************************************************************######
### Survival data #####
sur=read.csv("data/survival.csv")
head(sur)

# # Correct site 4, they did surviva to 5
# sur[sur$site==4,]$X5_survival<-1

# Make survival year 5 of NA
sur$X5_survival

## Survival predicted
popsize=read.csv("pop_size_estimation/pop_size_estimations.csv")
popsizeraw=popsize

# The predictability
stabilizing_prs %>% head

# Merge
predsur<-merge(sur, stabilizing_prs,by=c("site","plot")) %>%
  dplyr::filter(X3_survival >-1)

predsur %>%
  dplyr::filter(X3_survival >-1) %>%
ggplot(.)+
  geom_point(aes(y=X3_survival, x=sp_r))

# Merge and average 3 year
predsur<-merge(sur, stabilizing_prs,by=c("site","plot")) %>%
  dplyr::filter(X3_survival >-1) %>%
  group_by(site) %>%
  summarise(
              meansur= mean(X3_survival[X3_survival %in% c(0,1)],na.rm=T),
              meansur5= mean(X5_survival[X5_survival %in% c(0,1)],na.rm=T),
             meanprs=mean(sp_r,na.rm=T),
             meanr2=mean(r2,na.rm=T) )
predsur %>%
  ggplot(.)+
  geom_point(aes(y=meansur5, x=meanprs))+
  stat_smooth(aes(y=meansur5, x=meanprs),method="glm")+
  xlab("Mean stabilizing prs LOO rho (y1)")+
  ylab("Survival proportion of plots (y5)")

predsur %>%
ggplot(.)+
  geom_point(aes(y=meansur, x=meanprs))+
  stat_smooth(aes(y=meansur, x=meanprs),method="glm")+
  xlab("Mean stabilizing prs LOO rho (y1)")+
  ylab("Survival proportion of plots (y3)")

sink("tables/test-correlation-prsstabilizing-survival3years.txt")
cat("
Test of correlation between the stabilizing selection PRS
and survival at 3rd year.
Both are averaged across plots
"
)
cor.test(predsur$meanprs, predsur$meansur)
sink()

t.test(predsur$meanr2>0.1, predsur$meansur)


# try the GO
predsur<-merge(sur, stabilizing_prs,by=c("site","plot")) %>%
  dplyr::filter(X3_survival >-1) %>%
  group_by(site) %>%
  summarise(
    meansur= mean(X3_survival,na.rm=T),
    meansur5= mean(X3_survival,na.rm=T),
             meanprs=mean(pearson_r,na.rm=T) )
predsur %>%
  ggplot(.)+
  geom_point(aes(y=meansur, x=meanprs))
cor.test(predsur$meanprs, predsur$meansur)


# Merge and average 1 year
predsur<-merge(sur, stabilizing_prs,by=c("site","plot")) %>%
  dplyr::filter(X1_survival >-1) %>%
  group_by(site) %>%
  summarise( meansur= mean(X1_survival,na.rm=T),
             meanprs=mean(sp_r,na.rm=T) )
predsur %>%
  ggplot(.)+
  geom_point(aes(y=meansur, x=meanprs))

# Merge and average 2 year
predsur<-merge(sur, stabilizing_prs,by=c("site","plot")) %>%
  dplyr::filter(X2_survival >-1) %>%
  group_by(site) %>%
  summarise( meansur= mean(X2_survival,na.rm=T),
             meanprs=mean(sp_r,na.rm=T) )
predsur %>%
  ggplot(.)+
  geom_point(aes(y=meansur, x=meanprs))

cor.test(predsur$meanprs, predsur$meansur)


# Merge and average 4 year
predsur<-merge(sur, stabilizing_prs,by=c("site","plot")) %>%
  dplyr::filter(X4_survival >-1) %>%
  group_by(site) %>%
  summarise( meansur= mean(X4_survival,na.rm=T),
             meanprs=mean(sp_r,na.rm=T) )
predsur %>%
  ggplot(.)+
  geom_point(aes(y=meansur, x=meanprs))
cor.test(predsur$meanprs, predsur$meansur)

# Merge and average 5 year
predsur<-merge(sur, stabilizing_prs,by=c("site","plot")) %>%
  dplyr::filter(X5_survival >-1) %>%
  group_by(site) %>%
  summarise( meansur= mean(X5_survival,na.rm=T),
             meanprs=mean(sp_r,na.rm=T) )
predsur %>%
  ggplot(.)+
  geom_point(aes(y=meansur, x=meanprs))
cor.test(predsur$meanprs, predsur$meansur)



# Check year 3 and 5
predsur<-merge(sur, stabilizing_prs,by=c("site","plot")) %>%
  group_by(site) %>%
  summarise(
    meansur5= mean(X5_survival[X5_survival>-1],na.rm=T),
    meansur4= mean(X4_survival[X4_survival>-1],na.rm=T),
    meansur3= mean(X3_survival[X3_survival>-1],na.rm=T),
    meansur2= mean(X2_survival[X2_survival>-1],na.rm=T),
    meansur2= mean(X2_survival[X2_survival>-1],na.rm=T),
    meansur1= mean(X2_survival[X1_survival>-1],na.rm=T),
    meanr2=mean(r2,na.rm=T) ,
    meanprs=mean(sp_r,na.rm=T) )

t.test(predsur$meanr2>0.1, predsur$meansur3)
t.test(predsur$meanr2>0.1, predsur$meansur5)

ggplot(predsur)+
  geom_point(aes(y=meansur3, x=meanr2))
ggplot(predsur)+
  geom_point(aes(y=meansur5, x=meanr2))

ggplot(predsur)+
  geom_point(aes(y=meansur5, x=meanprs))+
  stat_smooth(aes(y=meansur5, x=meanprs),method="glm")+
  xlab("Mean stabilizing prs LOO rho (y1)")+
  ylab("Survival proportion of plots (y5)")

# 3 vs 5
ggplot(predsur)+
  geom_segment(aes(
                   y=meansur3, x=meanprs,
                   yend=meansur5, xend=meanprs
  ), color="black", arrow = arrow(length = unit(0.1, "inches")))+
  geom_point(aes(y=meansur3, x=meanprs), color="red")+
  geom_point(aes(y=meansur5, x=meanprs), color="blue")

# 2 vs 3
ggplot(predsur)+
  geom_segment(aes(
    y=meansur2, x=meanprs,
    yend=meansur3, xend=meanprs
  ), color="black", arrow = arrow(length = unit(0.1, "inches")))+
  geom_point(aes(y=meansur2, x=meanprs), color="red")+
  geom_point(aes(y=meansur3, x=meanprs), color="blue")


# 2 vs 5
ggplot(predsur)+
  geom_segment(aes(
    y=meansur2, x=meanprs,
    yend=meansur5, xend=meanprs
  ), color="black", arrow = arrow(length = unit(0.1, "inches")))+
  geom_point(aes(y=meansur2, x=meanprs), color="red")+
  geom_point(aes(y=meansur5, x=meanprs), color="blue")

# 2 vs 4
ggplot(predsur)+
  geom_segment(aes(
    y=meansur2, x=meanprs,
    yend=meansur4, xend=meanprs
  ), color="black", arrow = arrow(length = unit(0.1, "inches")))+
  geom_point(aes(y=meansur2, x=meanprs), color="red")+
  geom_point(aes(y=meansur4, x=meanprs), color="blue")

# 4 vs 5
ggplot(predsur)+
  geom_segment(aes(
    y=meansur4, x=meanprs,
    yend=meansur5, xend=meanprs
  ), color="black", arrow = arrow(length = unit(0.1, "inches")))+
  geom_point(aes(y=meansur4, x=meanprs), color="red")+
  geom_point(aes(y=meansur5, x=meanprs), color="blue")


#### Subset survival data ####

##  data for logistic regressions
gonew<-read.csv("data-intermediate/go_sign_snps_from_founder_pop.csv")

mersur<-mersurclimate<-merge(sur, stabilizing_prs,by=c("site","plot")) %>% # OPTION ONE IGNORE LOST REPLICATES
  dplyr::mutate(
    X1_survival=ifelse(X1_survival>-1,X1_survival,NA),
    X2_survival=ifelse(X2_survival>-1,X2_survival,NA),
    X3_survival=ifelse(X3_survival>-1,X3_survival,NA),
    X4_survival=ifelse(X4_survival>-1,X4_survival,NA),
    X5_survival=ifelse(X5_survival>-1,X5_survival,NA)
                ) %>%
  merge(.,worldclim_sitesdata,by="site") %>%
  merge(.,gonew,by="site")
gonew<-read.csv("data-intermediate/go_sign_snps_from_founder_pop.csv")
mersurgo<-merge(sur, go,by=c("site","plot")) %>% # OPTION ONE IGNORE LOST REPLICATES.  <<<<<<<<<<, IMPORTNAT MAKES EVERYTHING GO PREDICTIONS
  dplyr::mutate(
    X1_survival=ifelse(X1_survival>-1,X1_survival,NA),
    X2_survival=ifelse(X2_survival>-1,X2_survival,NA),
    X3_survival=ifelse(X3_survival>-1,X3_survival,NA),
    X4_survival=ifelse(X4_survival>-1,X4_survival,NA),
    X5_survival=ifelse(X5_survival>-1,X5_survival,NA)
  )%>%
  merge(.,worldclim_sitesdata,by="site") %>%
  merge(.,gonew,by="site")
merrep<-merge(sur, replicates,by=c("site","plot")) %>% # OPTION ONE IGNORE LOST REPLICATES
  dplyr::mutate(
    X1_survival=ifelse(X1_survival>-1,X1_survival,NA),
    X2_survival=ifelse(X2_survival>-1,X2_survival,NA),
    X3_survival=ifelse(X3_survival>-1,X3_survival,NA),
    X4_survival=ifelse(X4_survival>-1,X4_survival,NA),
    X5_survival=ifelse(X5_survival>-1,X5_survival,NA)
  )%>%
  merge(.,worldclim_sitesdata,by="site")

#####*********************************************************************######
#### FIGS POPSIZE ####
plota<-
  ggplot(mersurclimate)+
  geom_point(aes(y=X1_flowerstotal+X2_flowerstotal+X3_flowerstotal,
                 # x=sp_r),
                 x=r2),
             color="black",size=2)+
  geom_point(aes(y=X1_flowerstotal+X2_flowerstotal+X3_flowerstotal,
                 # x=sp_r),
                 x=r2,
                 color=bio1), size=1.5)+
  # scale_y_log10()+
  # scale_x_log10()+
  scale_color_gradientn("",colours = rev(redblue))+
  xlab("Predictability year 1 (stabilizing r2 LOO)")+
  ylab("Population size summed year 1-3 (# adults)")+
  theme_minimal()+
  theme(legend.position = "none")
plota

plotb<-
  ggplot(mersurclimate)+
  geom_point(aes(y=X1_flowerstotal+X2_flowerstotal+X3_flowerstotal,
                 x=sp_r), color="black",size=2)+
  geom_point(aes(y=X1_flowerstotal+X2_flowerstotal+X3_flowerstotal,
                 x=sp_r,
                 color=bio1), size=1.5)+
  scale_y_log10()+
  scale_x_log10()+
  scale_color_gradientn("",colours = rev(redblue))+
  xlab("Predictability year 1 (stabilizing r2 LOO)")+
  ylab("Population size summed year 1-3 (# adults)")+
  theme_minimal()+
  theme(legend.position = "none")
plotb
# legend
temp_plot <- ggplot(mersurclimate) + # Create one plot with a visible legend to extract
  geom_point(aes(x = r2, y = X1_flowerstotal + X2_flowerstotal + X3_flowerstotal,
                 color = bio1)) +
  scale_color_gradientn("Temp. (C)", colours = rev(redblue)) +
  theme_minimal() +
  theme(legend.position = "right")  # show legend
legend_grob <- get_legend(temp_plot) # Extract legend as a grob
plota <- ggdraw() + # Overlay the legend grob in front of the plot
  draw_plot(plota) +
  draw_plot(legend_grob, x = 0.7, y = 0.6, width = 0.25, height = 0.25)  # adjust position here

# grid
predictabilitysize<-
  plot_grid(plota,plotb, labels= "AUTO")
predictabilitysize

save_plot(filename = "figs/fig-predictability-stabilizing-year-1-vs-popsize-allyears-bio1.pdf",
          predictabilitysize,
          base_height = 4,base_width = 8)
save_plot(filename = "figs/fig-predictability-stabilizing-year-1-vs-popsize-allyears-bio1.png",
          predictabilitysize,
          base_height = 4,base_width = 8)


save_plot(filename = "figs/fig-predictability-stabilizing-year-1-vs-popsize-allyears-bio1-a.pdf",
          plota,
          base_height = 2,base_width = 4)


#### FIGS LANDSCAPE ####
# Trials
df %>%
  dplyr::filter(source=="plot_mean") %>% # check out plot mean (repeatability)
  ggplot(.)+
  geom_point(
    aes(
      y=r2,
      # y=sp_r,
      x=bio1),
    color="black",
    alpha=1, size=3.5)+
  geom_point(
    aes(
      y=r2,
      # y=sp_r,
      x=bio1,
      color=bio1
    ),
    alpha=1, size=3)+
  stat_smooth(aes(
    y=r2,
    # y=sp_r,
    x=bio1),
    se=F, color="grey",
    method = "glm",formula=y~poly(x,2))+
  stat_smooth(aes(
    y=r2,
    # y=sp_r,
    x=bio1),
    se=F, color="grey",
    method = "glm",formula=y~poly(x,1))+
  scale_color_gradientn("",colours = rev(redblue))+
  theme_minimal()+
  coord_flip()

# Export
predictability<-
  df %>%
  # dplyr::filter(source=="genomic offset") %>%
  dplyr::filter(source=="stabilizing_prs") %>%
  ggplot(.)+
  geom_point(
    aes(
      y=r2,
      # y=sp_r,
      x=bio1),
    color="black",
    alpha=1, size=3.5)+
  geom_point(
    aes(
      y=r2,
      # y=sp_r,
      x=bio1,
      color=bio1
    ),
    alpha=1, size=3)+
  stat_smooth(aes(
    y=r2,
    # y=sp_r,
    x=bio1),
    se=F, color="grey",
    method = "glm",formula=y~poly(x,2))+
  stat_smooth(aes(
    y=r2,
    # y=sp_r,
    x=bio1),
    se=F, color="grey",
    method = "glm",formula=y~poly(x,1))+
  scale_color_gradientn("",colours = rev(redblue))+
  theme_minimal()+
  coord_flip()
predictability
predictability<-
  df %>%
  # dplyr::filter(source=="genomic offset") %>%
  dplyr::filter(source=="stabilizing_prs") %>%
  ggplot(.)+
  geom_point(
    aes(
      # y=r2,
      y=sp_r,
      x=bio1),
    color="black",
    alpha=1, size=3.5)+
  geom_point(
    aes(
      # y=r2,
      y=sp_r,
      x=bio1,
      color=bio1
    ),
    alpha=1, size=3)+
  stat_smooth(aes(
    # y=r2,
    y=sp_r,
    x=bio1),
    se=F, color="grey",
    method = "glm",formula=y~poly(x,2))+
  scale_color_gradientn("",colours = rev(redblue))+
  theme_minimal()+
  coord_flip()
predictability
save_plot(filename = "figs/fig-predictability-stabilizing_prs-year-1-bio1.pdf",
          predictability,
          base_height = 5,base_width = 6
)

predictability<-
  df %>%
  # dplyr::filter(source=="genomic offset") %>%
  dplyr::filter(source=="stabilizing_prs") %>%
  ggplot(.)+
  geom_point(
    aes(
      y=r2,
      # y=sp_r,
      x=bio1),
    color="black",
    alpha=1, size=3.5)+
  geom_point(
    aes(
      y=r2,
      # y=sp_r,
      x=bio1,
      color=bio1
    ),
    alpha=1, size=3)+
  stat_smooth(aes(
    y=r2,
    # y=sp_r,
    x=bio1),
    se=F, color="grey",
    method = "glm",formula=y~poly(x,2))+
  scale_color_gradientn("",colours = rev(redblue))+
  theme_minimal()+
  coord_flip()
predictability
save_plot(filename = "figs/fig-predictability-stabilizing_prs-year-1-bio1-r2.pdf",
          predictability,
          base_height = 5,base_width = 6
)

predictability<-
  df %>%
  # dplyr::filter(source=="genomic offset") %>%
  dplyr::filter(source=="genomic offset") %>%
  ggplot(.)+
  geom_point(
    aes(
      # y=r2,
      y=sp_r,
      x=bio1),
    color="black",
    alpha=1, size=3.5)+
  geom_point(
    aes(
      # y=r2,
      y=sp_r,
      x=bio1,
      color=bio1
    ),
    alpha=1, size=3)+
  stat_smooth(aes(
    # y=r2,
    y=sp_r,
    x=bio1),
    se=F, color="grey",
    method = "glm",formula=y~poly(x,2))+
  scale_color_gradientn("",colours = rev(redblue))+
  theme_minimal()+
  coord_flip()
predictability
predictability<-
  df %>%
  # dplyr::filter(source=="genomic offset") %>%
  dplyr::filter(source=="genomic offset") %>%
  ggplot(.)+
  geom_point(
    aes(
      y=sp_r,
      # y=r2,
      x=bio1),
    color="black",
    alpha=1, size=3.5)+
  geom_point(
    aes(
      y=sp_r,
      # y=r2,
      x=bio1,
      color=bio1
    ),
    alpha=1, size=3)+
  stat_smooth(aes(
    y=sp_r,
    # y=r2,
    x=bio1),
    se=F, color="grey",
    method = "glm",formula=y~poly(x,2))+
  scale_color_gradientn("",colours = rev(redblue))+
  theme_minimal()+
  coord_flip()
predictability
save_plot(filename = "figs/fig-predictability-genomicoffset-year-1-bio1.pdf",
          predictability,
          base_height = 5,base_width = 6
)
predictability<-
  df %>%
  # dplyr::filter(source=="genomic offset") %>%
  dplyr::filter(source=="genomic offset") %>%
  ggplot(.)+
  geom_point(
    aes(
      y=r2,
      x=bio1),
    color="black",
    alpha=1, size=3.5)+
  geom_point(
    aes(
      y=r2,
      x=bio1,
      color=bio1
    ),
    alpha=1, size=3)+
  stat_smooth(aes(
    y=r2,
    x=bio1),
    se=F, color="grey",
    method = "glm",formula=y~poly(x,2))+
  scale_color_gradientn("",colours = rev(redblue))+
  theme_minimal()+
  coord_flip()
predictability
save_plot(filename = "figs/fig-predictability-genomicoffset-year-1-bio1-r2.pdf",
          predictability,
          base_height = 5,base_width = 6
)


df %>%
  # dplyr::filter(source=="genomic offset") %>%
  dplyr::filter(source=="genomic offset") %>%
  ggplot(.)+
  geom_point(
    aes(
      # y=r2,
      y=sp_r,
      x=bio1),
    color="black",
    alpha=1, size=3.5)+
  geom_point(
    aes(
      # y=r2,
      y=sp_r,
      x=bio1,
      color=bio1
    ),
    alpha=1, size=3)+
  stat_smooth(aes(
    # y=r2,
    y=sp_r,
    x=bio1),
    se=F, color="grey",
    method = "glm",formula=y~poly(x,2))+
  scale_color_gradientn("",colours = rev(redblue))+
  theme_minimal()+
  coord_flip()

df %>%
  # dplyr::filter(source=="genomic offset") %>%
  dplyr::filter(source=="genomic offset") %>%
  ggplot(.)+
  geom_point(
    aes(
      # y=r2,
      y=r2,
      x=bio1),
    color="black",
    alpha=1, size=3.5)+
  geom_point(
    aes(
      # y=r2,
      y=r2,
      x=bio1,
      color=bio1
    ),
    alpha=1, size=3)+
  stat_smooth(aes(
    # y=r2,
    y=r2,
    x=bio1),
    se=F, color="grey",
    method = "glm",formula=y~poly(x,2))+
  scale_color_gradientn("",colours = rev(redblue))+
  theme_minimal()+
  coord_flip()

#####*********************************************************************######
##### LOGISTIC  REG ######

sink("tables/logistic-regressions-predictability-popsurvival.txt")
# model y2
modely2rho<-model <- glm(X2_survival ~ sp_r, data =mersur , family = binomial)
modely2r2<-model <- glm(X2_survival ~ r2, data =mersur , family = binomial)
summary(model)

mersur$X2_log_r2<-predict(modely2r2, newdata = mersur, type = "response")
mersur$X2_log_rho<-predict(modely2rho, newdata = mersur, type = "response")

# model y3
modely3rho<-model <- glm(X3_survival ~ sp_r, data =mersur , family = binomial)
modely3r2<-model <- glm(X3_survival ~ r2, data =mersur , family = binomial)
summary(model)

mersur$X3_log_r2<-predict(modely3r2, newdata = mersur, type = "response")
mersur$X3_log_rho<-predict(modely3rho, newdata = mersur, type = "response")

# model y4
modely4rho<-model <- glm(X4_survival ~ sp_r, data =mersur , family = binomial)
modely4r2<-model <- glm(X4_survival ~ r2, data =mersur , family = binomial)
summary(model)

mersur$X4_log_r2<-predict(modely4r2, newdata = mersur, type = "response")
mersur$X4_log_rho<-predict(modely4rho, newdata = mersur, type = "response")

# model y5
modely5rho<-model <- glm(X5_survival ~ sp_r, data =mersur , family = binomial)
modely5r2<-model <- glm(X5_survival ~ r2, data =mersur , family = binomial)
summary(model)

mersur$X5_log_r2<-predict(modely5r2, newdata = mersur, type = "response")
mersur$X5_log_rho<-predict(modely5rho, newdata = mersur, type = "response")

sink()


sink("tables/logistic-regressions-predictability-popsurvival-genomicoffset.txt")
# model y2
modely2rho<-model <- glm(X2_survival ~ sp_r, data =mersurgo , family = binomial)
modely2r2<-model <- glm(X2_survival ~ r2, data =mersurgo , family = binomial)
summary(modely2rho)
summary(modely2r2)

mersurgo$X2_log_r2<-predict(modely2r2, newdata = mersurgo, type = "response")
mersurgo$X2_log_rho<-predict(modely2r2, newdata = mersurgo, type = "response")

# model y3
modely3rho<-model <- glm(X3_survival ~ sp_r, data =mersurgo , family = binomial)
modely3r2<-model <- glm(X3_survival ~ r2, data =mersurgo , family = binomial)
summary(modely3rho)
summary(modely3r2)

mersurgo$X3_log_r2<-predict(modely3r2, newdata = mersurgo, type = "response")
mersurgo$X3_log_rho<-predict(modely3rho, newdata = mersurgo, type = "response")

# model y4
modely4rho<-model <- glm(X4_survival ~ sp_r, data =mersurgo , family = binomial)
modely4r2<-model <- glm(X4_survival ~ r2, data =mersurgo , family = binomial)
summary(modely4rho)
summary(modely4r2)

mersurgo$X4_log_r2<-predict(modely4r2, newdata = mersurgo, type = "response")
mersurgo$X4_log_rho<-predict(modely4rho, newdata = mersurgo, type = "response")

# model y5
modely5rho<-model <- glm(X5_survival ~ sp_r, data =mersurgo , family = binomial)
modely5r2<-model <- glm(X5_survival ~ r2, data =mersurgo , family = binomial)
summary(modely5rho)
summary(modely5r2)


mersurgo$X5_log_r2<-predict(modely5r2, newdata = mersurgo, type = "response")
mersurgo$X5_log_rho<-predict(modely5rho, newdata = mersurgo, type = "response")

sink()

########### LOGISTIC REG MULTIVARIATE ####
# Check logistic with plot replicate
glm(X2_survival ~ r2, data =merrep , family = binomial) %>% summary
glm(X3_survival ~ r2, data =merrep , family = binomial) %>% summary
glm(X4_survival ~ r2, data =merrep , family = binomial) %>% summary
glm(X5_survival ~ r2, data =merrep , family = binomial) %>% summary


# Check logistic with just temperature
glm(X3_survival ~ 1, data =mersurclimate , family = binomial) %>% summary
glm(X3_survival ~ bio1, data =mersurclimate , family = binomial) %>% summary
glm(X3_survival ~ sp_r, data =mersurclimate , family = binomial) %>% summary
glm(X3_survival ~ r2, data =mersurclimate , family = binomial) %>% summary
glm(X3_survival ~ r2 + bio1, data =mersurclimate , family = binomial) %>% summary
glm(X3_survival ~ r2 + bio1 , data =mersurclimate , family = binomial) %>% summary


glm(X5_survival ~ bio1, data =mersurclimate , family = binomial) %>% summary
glm(X5_survival ~ r2 + bio1, data =mersurclimate , family = binomial) %>% summary
glm(X5_survival ~ r2 * bio1, data =mersurclimate , family = binomial) %>% summary


glm(X3_survival ~ r2 + bio1 , data =dplyr::filter(mersurclimate,bio1>5) , family = binomial) %>% summary


glm(X5_survival ~ r2 * go, data =mersurclimate , family = binomial) %>% summary


glm(X2_survival ~ r2 * bio1, data =mersurclimate , family = binomial) %>% summary
glm(X3_survival ~ r2 * bio1, data =mersurclimate , family = binomial) %>% summary
glm(X4_survival ~ r2 * bio1, data =mersurclimate , family = binomial) %>% summary
glm(X5_survival ~ r2 * bio1, data =mersurclimate , family = binomial) %>% summary


glm(X2_survival ~ r2 + bio1, data =mersurgo , family = binomial) %>% summary
glm(X3_survival ~ r2 + bio1, data =mersurgo , family = binomial) %>% summary
glm(X4_survival ~ r2 + bio1, data =mersurgo , family = binomial) %>% summary
glm(X5_survival ~ r2 * bio1, data =mersurgo , family = binomial) %>% summary

#### Figs
pdf("figs/fig-predictability-r2-year-1-correlation-survival-populations-across-years-r2.pdf")

ggplot(mersur)+
  geom_jitter(aes(y=X2_survival,x=r2),
              height = 0.025)+
  geom_line(aes(y=X2_log_r2,x=r2))+
  xlab("Predictability year 1 (stabilizing spr LOO)")+
  ylab("Survival at year 2 (survival/death)")+
  theme_minimal()

ggplot(mersur)+
  geom_jitter(aes(y=X3_survival,x=r2),
              height = 0.025)+
  geom_line(aes(y=X3_log_r2,x=r2))+
  xlab("Predictability year 1 (stabilizing spr LOO)")+
  ylab("Survival at year 3 (survival/death)")+
  theme_minimal()

ggplot(mersur)+
  geom_jitter(aes(y=X4_survival,x=r2),
              height = 0.025)+
  geom_line(aes(y=X4_log_r2,x=r2))+
  xlab("Predictability year 1 (stabilizing spr LOO)")+
  ylab("Survival at year 4 (survival/death)")+
  theme_minimal()


ggplot(mersur)+
  geom_jitter(aes(y=X5_survival,x=r2),
              height = 0.025)+
  geom_line(aes(y=X5_log_r2,x=r2))+
  xlab("Predictability year 1 (stabilizing spr LOO)")+
  ylab("Survival at year 4 (survival/death)")+
  theme_minimal()

dev.off()


pdf("figs/fig-predictability-year-1-correlation-survival-populations-across-years-spearmanr.pdf")
# ggplot(mersur)+
#   geom_point(aes(y=X1_flowerstotal,x=sp_r))+
#   geom_line(aes(y=X1_log,x=sp_r))+
#   xlab("Predictability year 1 (stabilizing spr LOO)")+
#   ylab("Survival at year 2 (flowers total)")+
#   theme_minimal()


ggplot(mersur)+
  geom_jitter(aes(y=X2_survival,x=sp_r),
              height = 0.025)+
  geom_line(aes(y=X2_log_rho,x=sp_r))+
  xlab("Predictability year 1 (stabilizing spr LOO)")+
  ylab("Survival at year 2 (survival/death)")+
  theme_minimal()

ggplot(mersur)+
  geom_jitter(aes(y=X3_survival,x=sp_r),
             height = 0.025)+
  geom_line(aes(y=X3_log_rho,x=sp_r))+
  xlab("Predictability year 1 (stabilizing spr LOO)")+
  ylab("Survival at year 3 (survival/death)")+
  theme_minimal()

ggplot(mersur)+
  geom_jitter(aes(y=X4_survival,x=sp_r),
             height = 0.025)+
  geom_line(aes(y=X4_log_rho,x=sp_r))+
  xlab("Predictability year 1 (stabilizing spr LOO)")+
  ylab("Survival at year 4 (survival/death)")+
  theme_minimal()

ggplot(mersur)+
  geom_jitter(aes(y=X5_survival,x=sp_r),
             height = 0.025)+
  geom_line(aes(y=X5_log_rho,x=sp_r))+
  xlab("Predictability year 1 (stabilizing spr LOO)")+
  ylab("Survival at year 5 (survival/death)")+
  theme_minimal()
dev.off()


pdf("figs/fig-predictability-r2-year-1-correlation-survival-populations-across-years-r2-genomicoffset.pdf")

ggplot(mersurgo)+
  geom_jitter(aes(y=X2_survival,x=r2),
              height = 0.025)+
  geom_line(aes(y=X2_log_r2,x=r2))+
  xlab("Predictability year 1 (stabilizing spr LOO)")+
  ylab("Survival at year 2 (survival/death)")+
  theme_minimal()

ggplot(mersurgo)+
  geom_jitter(aes(y=X3_survival,x=r2),
              height = 0.025)+
  geom_line(aes(y=X3_log_r2,x=r2))+
  xlab("Predictability year 1 (stabilizing spr LOO)")+
  ylab("Survival at year 3 (survival/death)")+
  theme_minimal()

ggplot(mersurgo)+
  geom_jitter(aes(y=X4_survival,x=r2),
              height = 0.025)+
  geom_line(aes(y=X4_log_r2,x=r2))+
  xlab("Predictability year 1 (stabilizing spr LOO)")+
  ylab("Survival at year 4 (survival/death)")+
  theme_minimal()


ggplot(mersurgo)+
  geom_jitter(aes(y=X5_survival,x=r2),
              height = 0.025)+
  geom_line(aes(y=X5_log_r2,x=r2))+
  xlab("Predictability year 1 (stabilizing spr LOO)")+
  ylab("Survival at year 4 (survival/death)")+
  theme_minimal()

dev.off()


#  ##### LOGISTIC QUADRATICS ######
#
# # Model  quadratics
# mersur$r2_sq <- mersur$r2^2
# mersur$rho_sq <- mersur$sp_r^2
#
# glm_quadratic <- glm(X2_survival ~ r2 + r2_sq, family = binomial, data = mersur)
# summary(glm_quadratic)
# mersur$X2_log_quad<-predict(glm_quadratic, newdata = mersur, type = "response")
#
#
# glm_quadratic <- glm(X3_survival ~ r2 + r2_sq, family = binomial, data = mersur)
# summary(glm_quadratic)
# mersur$X3_log_quad<-predict(glm_quadratic, newdata = mersur, type = "response")
#
# glm_quadratic <- glm(X4_survival ~ r2 + r2_sq, family = binomial, data = mersur)
# summary(glm_quadratic)
# mersur$X4_log_quad<-predict(glm_quadratic, newdata = mersur, type = "response")
#
# glm_quadratic <- glm(X5_survival ~ r2 + r2_sq, family = binomial, data = mersur)
# summary(glm_quadratic)
# mersur$X5_log_quad<-predict(glm_quadratic, newdata = mersur, type = "response")

# # does this happen in the sites that have 5 years in gen 2 and 3?
# siteswithcontinousdata<-summarize(group_by(mersur,site), continuoussites=any(X5_survival==1,na.rm = T))
# mersurcontinuous<-
#   mersur %>%
#   dplyr::filter(site %in% siteswithcontinousdata$site[siteswithcontinousdata[,2]==TRUE])
#
# ggplot(mersurcontinuous)+
#   geom_point(aes(y=X2_survival,x=r2))+
#   stat_smooth(aes(y=X2_survival,x=r2),method="glm",formula=y~poly(x,2),se=F)+
#     xlab("Predictability year 1 (stabilizing spr LOO)")+
#   ylab("Survival at year 2 (survival/death)")+
#   theme_minimal()+
#   labs(title="year 2 quadratic")
# ggplot(mersurcontinuous)+
#   geom_point(aes(y=X3_survival,x=r2))+
#   stat_smooth(aes(y=X3_survival,x=r2),method="glm",formula=y~poly(x,2),se=F)+
#   xlab("Predictability year 1 (stabilizing spr LOO)")+
#   ylab("Survival at year 3 (survival/death)")+
#   theme_minimal()+
#   labs(title="year 3 quadratic")
# ggplot(mersurcontinuous)+
#   geom_point(aes(y=X5_survival,x=r2))+
#   stat_smooth(aes(y=X5_survival,x=r2),method="glm",formula=y~poly(x,2))+
#   # geom_line(aes(y=X2_log_quad,x=r2))+
#   xlab("Predictability year 1 (stabilizing spr LOO)")+
#   ylab("Survival at year 5 (survival/death)")+
#   theme_minimal()+
#   labs(title="year 5 quadratic")
#
# ggplot(mersurcontinuous)+
#   geom_point(aes(y=X5_survival,x=sp_r))+
#   stat_smooth(aes(y=X5_survival,x=sp_r),method="glm",formula=y~poly(x,2))+
#   # geom_line(aes(y=X2_log_quad,x=r2))+
#   xlab("Predictability year 1 (stabilizing spr LOO)")+
#   ylab("Survival at year 5 (survival/death)")+
#   theme_minimal()+
#   labs(title="year 5 quadratic")
#
# ggplot(mersurcontinuous)+
#   geom_jitter(aes(y=log(X5_flowerstotal+1),x=r2), height = 0.5, width = 0.05)+
#   stat_smooth(aes(y=log(X5_flowerstotal+1),x=r2),method="glm",formula=y~poly(x,2))+
#   # geom_line(aes(y=X2_log_quad,x=r2))+
#   xlab("Predictability year 1 (stabilizing spr LOO)")+
#   ylab("Survival at year 5 (survival/death)")+
#   theme_minimal()+
#   labs(title="year 5 quadratic")
#
#
# #### Figs quadratic
# pdf("figs/fig-predictability-r2-year-1-correlation-survival-populations-across-years-r2-quadratic.pdf")
#
# ggplot(mersur)+
#   geom_point(aes(y=X2_survival,x=r2))+
#   geom_line(aes(y=X2_log_quad,x=r2))+
#   xlab("Predictability year 1 (stabilizing spr LOO)")+
#   ylab("Survival at year 2 (survival/death)")+
#   theme_minimal()+
#   labs(title="year 2 quadratic")
#
# ggplot(mersur)+
#   geom_point(aes(y=X3_survival,x=r2))+
#   geom_line(aes(y=X3_log_quad,x=r2))+
#   xlab("Predictability year 1 (stabilizing spr LOO)")+
#   ylab("Survival at year 3 (survival/death)")+
#   theme_minimal()+
#   labs(title="year 3 quadratic")
#
# ggplot(mersur)+
#   geom_point(aes(y=X4_survival,x=r2))+
#   geom_line(aes(y=X4_log_quad,x=r2))+
#   xlab("Predictability year 1 (stabilizing spr LOO)")+
#   ylab("Survival at year 4 (survival/death)")+
#   theme_minimal()+
#   labs(title="year 4 quadratic")
#
# ggplot(mersur)+
#   geom_point(aes(y=X5_survival,x=r2))+
#   geom_line(aes(y=X5_log_quad,x=r2))+
#   xlab("Predictability year 1 (stabilizing spr LOO)")+
#   ylab("Survival at year 5 (survival/death)")+
#   theme_minimal()+
#   labs(title="year 5 quadratic")
#
# dev.off()




########### LOGISTIC PLOTS ####

# LOGISTIC PLOTS
# # Function to create a logistic curve layer
add_logistic_curve <- function(intercept, slope, xrange = c(0, 1)) {
  # Create a data frame with r2 values and predicted probabilities
  curve_data <- data.frame(
    r2 = seq(xrange[1], xrange[2], length.out = 200)
  )
  curve_data$probability <- 1 / (1 + exp(-(intercept + slope * curve_data$r2)))

  # Return the geom_line layer
  geom_line(data = curve_data, aes(x = r2, y = probability), size = 1.2)
}


myplot <-
  ggplot(mersurclimate) +
  geom_jitter(aes(y = X2_survival,
                  # x=r2,
                  x=sp_r,
                  color=bio1),
              height=0.05)+
  scale_color_gradientn("", colours = rev(redblue)) +
  add_logistic_curve(modely2rho$coefficients[1],
                     modely2rho$coefficients[2] ,
                     # xrange = range(mersurclimate$r2)
                     xrange = range(mersurclimate$sp_r)
                     ) +  # <- now adds properly
  # boxplot
  stat_boxplot(aes(x=sp_r, group=X2_survival, y=X2_survival), alpha=0)+
  xlab("Predictability year 1 (stabilizing r2 LOO)") +
  ylab("Survival at year 2 (survival/death)") +
  theme_minimal() +
  theme(legend.position = "none")
myplot
save_plot(myplot,
        filename=  "figs/fig-predictability-year-1-correlation-survival-year-2-populations-across-years-bio1.pdf",
    base_height=2, base_width=3)

myplot <-
ggplot(mersurclimate)+
  geom_jitter(aes(y=X3_survival,
                  x=sp_r,
                  color=bio1),
              height=0.05)+
  scale_color_gradientn("",colours = rev(redblue))+
  add_logistic_curve(modely3rho$coefficients[1],
                     modely3rho$coefficients[2] ,
                     xrange = range(mersurclimate$sp_r)
  ) +  # <- now adds properly
  # boxplot
  stat_boxplot(aes(x=sp_r, group=X3_survival, y=X3_survival), alpha=0)+
  xlab("Predictability year 1 (stabilizing r2 LOO)")+
  ylab("Survival at year 3 (survival/death)")+
  theme_minimal()+
  theme(legend.position = "none")
myplot
save_plot(myplot,
         filename= "figs/fig-predictability-year-1-correlation-survival-year-3-populations-across-years-bio1.pdf",
          base_height=2, base_width=3)


myplot <-
  ggplot(mersurclimate) +
  geom_jitter(aes(y = X2_survival,
                  x=r2,
                  color=bio1),
              height=0.05)+
  scale_color_gradientn("", colours = rev(redblue)) +
  add_logistic_curve(modely2r2$coefficients[1],
                     modely2r2$coefficients[2] ,
                     xrange = range(mersurclimate$r2)
  ) +
  # boxplot
  stat_boxplot(aes(
                  x=r2,
                   group=X2_survival, y=X2_survival), alpha=0)+
  xlab("Predictability year 1 (stabilizing r2 LOO)") +
  ylab("Survival at year 2 (survival/death)") +
  theme_minimal() +
  theme(legend.position = "none")
myplot
save_plot(myplot,
          filename=  "figs/fig-predictability-year-1-correlation-r2-survival-year-2-populations-across-years-bio1.pdf",
          base_height=2, base_width=3)

myplot <-
  ggplot(mersurclimate)+
  geom_jitter(aes(y=X3_survival,
                  x=r2,
                  color=bio1),
              height=0.05)+
  scale_color_gradientn("",colours = rev(redblue))+
  add_logistic_curve(modely3r2$coefficients[1],
                     modely3r2$coefficients[2] ,
                     xrange = range(mersurclimate$r2)
  ) +  # <- now adds properly
  # boxplot
  stat_boxplot(aes(
    x=sp_r,
                   x=r2,
                   group=X3_survival, y=X3_survival), alpha=0)+
  xlab("Predictability year 1 (stabilizing r2 LOO)")+
  ylab("Survival at year 3 (survival/death)")+
  theme_minimal()+
  theme(legend.position = "none")
myplot
save_plot(myplot,
          filename= "figs/fig-predictability-year-1-correlation-r2-survival-year-3-populations-across-years-bio1.pdf",
          base_height=2, base_width=3)




#####*********************************************************************######
##### Heatmap visualize logistic reg ####

dat<-mersurgo
dat<-mersur
dat<-merrep

# dat$bio1<-mersurgo$go

# Grid for 2D prediction
grid <- expand.grid(
  r2 = seq(min(dat$r2), max(dat$r2), length.out = 100),
  bio1 = seq(min(dat$bio1), max(dat$bio1), length.out = 100)
)


# YEAR 2
grid$prob <- predict(glm(X2_survival ~ r2 * bio1, data = dat, family = binomial),
                     newdata = grid, type = "response")
# Heatmap
p2<-
  ggplot(grid, aes(x = r2, y = bio1, fill = prob)) +
  geom_tile() +
  # scale_fill_gradient(low = "blue", high = "red") +
  scale_fill_gradientn(colors=redblue) +
  labs(fill = "Survival",title="year 2") +
  theme_minimal()

# YEAR 3
grid$prob <- predict(glm(X3_survival ~ r2 * bio1, data = dat, family = binomial),
                     newdata = grid, type = "response")
# Heatmap
p3<-
  ggplot(grid, aes(x = r2, y = bio1, fill = prob)) +
  geom_tile() +
  # scale_fill_gradient(low = "blue", high = "red") +
  scale_fill_gradientn(colors=redblue) +
  labs(fill = "Survival",title="year 3") +
  theme_minimal()

# YEAR 4

grid$prob <- predict(glm(X4_survival ~ r2 * bio1, data = dat, family = binomial),
                     newdata = grid, type = "response")
# Heatmap
p4<-
  ggplot(grid, aes(x = r2, y = bio1, fill = prob)) +
  geom_tile() +
  # scale_fill_gradient(low = "blue", high = "red") +
  scale_fill_gradientn(colors=redblue) +
  labs(fill = "Survival",title="year 4") +
  theme_minimal()



# YEAR 5
grid$prob <- predict(glm(X5_survival ~ r2 * bio1, data = dat, family = binomial),
                     newdata = grid, type = "response")
# Heatmap
p5<-
  ggplot(grid, aes(x = r2, y = bio1, fill = prob)) +
  geom_tile() +
  # scale_fill_gradient(low = "blue", high = "red") +
  scale_fill_gradientn(colors=redblue) +
  labs(fill = "Survival",title="year 5") +
  theme_minimal()

heatmapsur<-
  plot_grid(p2,p3,p4,p5, ncol=4)
heatmapsur

p3b<-
  ggplot(grid, aes(x = r2, y = bio1, fill = factor(prob>0.5)) )+
  geom_tile() +
  labs(fill = "Survival",title="year 3") +
  theme_minimal()
p3b


#
# save_plot(heatmapsur,
#           filename= "figs/fig-predictability-year-1-heatmap-allyears-survival.pdf",
#           base_height=3, base_width=8)

#####*********************************************************************######
##### LOGISTIC Combined years  #####

datlong <-
  # mersur %>% # stabilizng selection <<<< REPLICATES~!!!
  # mersurgo %>% # choose stabilizing or go <<<< REPLICATES~!!!
  merrep %>% # replicates <<<< REPLICATES~!!!
  pivot_longer(
    cols = matches("^X\\d+_survival$"),       # Match columns like X1_survival, X2_survival, etc.
    names_to = "year",                        # Create new column 'year'
    names_pattern = "X(\\d+)_survival",       # Extract the number as year
    values_to = "survival"                    # Values go into 'survival' column
  ) %>%
  mutate(year = as.integer(year))             # Convert year to integer


glm(formula=survival ~ r2 + bio1 + year + (1|site), data = datlong, family = binomial) %>%   summary
glmer(survival ~ r2 * bio1 + year+ (1 | site), data = datlong, family = binomial) %>%summary
glmer(survival ~ r2 * bio1 + (1 | year) + (1 | site), data = datlong, family = binomial) %>%summary


alllogistics<-
ggplot(datlong)+
  geom_jitter(aes(y=survival,x=r2, color=bio1),height = 0.1, width = 0)+
  stat_smooth(aes(y=survival,x=r2),method="glm",method.args = list(family = "binomial"),
              se=F, color="black")+
  scale_color_gradientn(colors=rev(redblue)) +
  facet_wrap(~paste0("year ",year) + ~ifelse(bio1>12,"hot","cold"), ncol=2)+
  xlab(TeX("Predictability (repeatability LOO $r^2$)"))+
  scale_y_continuous(breaks = c(0,1),labels = c(0,1),minor_breaks = c(0,1))+
  ylab("Survival")+
  theme_minimal()
alllogistics

save_plot(filename = "figs/fig-predictability-r2-year-1-correlation-survival-populations-allyears-replicability.pdf",
          alllogistics,
          base_height = 8,base_width = 6)
save_plot(filename = "figs/fig-predictability-r2-year-1-correlation-survival-populations-allyears-large-replicability.pdf",
          alllogistics,
          base_height = 10,base_width = 12)


datlong <-
  mersur %>% # stabilizng selection <<<< REPLICATES~!!!
  # mersurgo %>% # choose stabilizing or go <<<< REPLICATES~!!!
  # merrep %>% # replicates <<<< REPLICATES~!!!
  pivot_longer(
    cols = matches("^X\\d+_survival$"),       # Match columns like X1_survival, X2_survival, etc.
    names_to = "year",                        # Create new column 'year'
    names_pattern = "X(\\d+)_survival",       # Extract the number as year
    values_to = "survival"                    # Values go into 'survival' column
  ) %>%
  mutate(year = as.integer(year))             # Convert year to integer


glm(formula=survival ~ r2 + bio1 + year + (1|site), data = datlong, family = binomial) %>%   summary
glmer(survival ~ r2 * bio1 + year+ (1 | site), data = datlong, family = binomial) %>%summary
glmer(survival ~ r2 * bio1 + (1 | year) + (1 | site), data = datlong, family = binomial) %>%summary


alllogistics<-
  ggplot(datlong)+
  geom_jitter(aes(y=survival,x=r2, color=bio1),height = 0.1, width = 0)+
  stat_smooth(aes(y=survival,x=r2),method="glm",method.args = list(family = "binomial"),
              se=F, color="black")+
  scale_color_gradientn(colors=rev(redblue)) +
  facet_wrap(~paste0("year ",year) + ~ifelse(bio1>12,"hot","cold"), ncol=2)+
  xlab(TeX("Predictability (stab. LOO $r^2$)"))+
  scale_y_continuous(breaks = c(0,1),labels = c(0,1),minor_breaks = c(0,1))+
  ylab("Survival")+
  theme_minimal()
alllogistics

save_plot(filename = "figs/fig-predictability-r2-year-1-correlation-survival-populations-allyears.pdf",
          alllogistics,
          base_height = 8,base_width = 6)
save_plot(filename = "figs/fig-predictability-r2-year-1-correlation-survival-populations-allyears-large.pdf",
          alllogistics,
          base_height = 10,base_width = 12)


#### Log LOGISTIC


glmsummary<-function(survival, r2, bio1, site, output="r2", col=4){
  d<-data.frame(survival, r2, bio1)
  coefficients(summary(glm(survival ~ r2 * bio1 + (1|site),  family = binomial, data=d)))[output,col]
}

sink("tables/logistic-regressions-predictability-popsurvival-stabilizing-r2-bio1-peryear-warmonly.txt")
# sink("tables/logistic-regressions-predictability-popsurvival-stabilizing-r2-bio1-peryear.txt")
datlong %>%
  # dplyr::filter(bio1>11) %>%
  group_by(year) %>%
  summarise(
    r2beta=glmsummary(survival, r2, bio1, site, output="r2",1),
    bio1beta=glmsummary(survival, r2, bio1, site, output="bio1",1),
    bio1r2beta=glmsummary(survival, r2, bio1, site, output="r2:bio1",1),

    r2p=glmsummary(survival, r2, bio1, site, output="r2"),
    bio1p=glmsummary(survival, r2, bio1, site, output="bio1"),
    bio1r2p=glmsummary(survival, r2, bio1, site, output="r2:bio1")
  )
sink()

datlong %>%
  dplyr::filter(bio1>10) %>% # only hot!
  group_by(year) %>%
  summarise(
    r2beta=glmsummary(survival, r2, bio1, site, output="r2",1),
    bio1beta=glmsummary(survival, r2, bio1, site, output="bio1",1),
    bio1r2beta=glmsummary(survival, r2, bio1, site, output="r2:bio1",1),

    r2p=glmsummary(survival, r2, bio1, site, output="r2"),
    bio1p=glmsummary(survival, r2, bio1, site, output="bio1"),
    bio1r2p=glmsummary(survival, r2, bio1, site, output="r2:bio1")
  )

merrep %>%
  # dplyr::filter(bio1 > 10) %>%
  glmer(X3_survival ~ r2 * bio1 + (1 | site), data = ., family = binomial) %>%
  summary()

#####*********************************************************************######
##### PREDICTABILITY Va (w) #####

load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize_flowers.rda")
e = ecotype_allyears_frequencies_long_raw_climate_popsize_flowers

# Check that it works
esub<-e %>%
  dplyr::filter(year==1)
library(lme4)
mod<-
  lmer(
    formula= log((freq+0.00001)/startfreq) ~ 1 + (1|id),
    data = esub
  )
id<-as.data.frame(VarCorr(mod),comp="Variance")[1,"vcov"]
res<-as.data.frame(VarCorr(mod),comp="Variance")[2,"vcov"]
id/(id+res)

geth2freq<-function(freq,startfreq,id){
  mydata<-data.frame(
    freq,startfreq,id
  )
  mod<-
    lmer(
      # formula= log((freq+0.00001)/startfreq) ~ 1 + (1|id), #LOG
      formula= ((freq)/startfreq) ~ 1 + (1|id), #NOT LOG
      data = mydata
    )
  id<-as.data.frame(VarCorr(mod),comp="Variance")[1,"vcov"]
  res<-as.data.frame(VarCorr(mod),comp="Variance")[2,"vcov"]
  id/(id+res)
}
# Now per garden

myres<-c()
for(i in unique(esub$site)){
  esubsub<-
    esub %>% dplyr::filter(site==i)
  # esubsub$y<-log((esubsub$freq+0.001)/esubsub$startfreq)
  if(length(unique(esubsub$rep))>1){
    mod<-
      lmer(
        formula= log((freq+0.001)/startfreq) ~ 1 + (1|id),
        # formula= y ~ 1 + (1|id),
        data = esubsub
      )
    id<-as.data.frame(VarCorr(mod),comp="Variance")[1,"vcov"]
    res<-as.data.frame(VarCorr(mod),comp="Variance")[2,"vcov"]
    myres<-c(myres,id/(id+res))
  }else{
    myres<-c(myres,NA)
  }
}
myva<-data.frame(site=unique(esub$site),pseudoH2=myres)

vasur<-merge(myva,sur,by="site")
vasur<-merge(myva,mersurclimate,by="site")

vasurmean<-
  vasur %>%
  # group_by(site, pseudoH2) %>%
  group_by(site, pseudoH2,bio1) %>%
  summarise(
    predmean=mean(r2),
    survival1=mean(X1_survival>0),
    survival2=mean(X2_survival>0),
    survival3=mean(X3_survival>0),
    survival4=mean(X4_survival>0),
    survival5=mean(X5_survival>0),
    # pop=sum(X2_flowerstotal+X3_flowerstotal)
    pop=sum(X3_flowerstotal)
  )


ggplot(vasurmean)+
  geom_point(aes(x=pseudoH2,y=predmean,color=bio1),shape=16,size=3)+
  stat_smooth(aes(x=pseudoH2,y=predmean), method="glm", formula=y~poly(x,2), color="grey",se=F)+
  scale_color_gradientn("",colours = rev(redblue))+
  ylab("Predictability year 1 (median stabilizing r LOO)")+
  xlab("Pseudo-h2 of frequency of accessions")+
  theme_minimal()+
  theme(legend.position = "none")

ggplot(vasurmean)+
  geom_point(aes(x=pseudoH2,y=predmean,color=bio1, size=survival2),
             shape=16)+
  stat_smooth(aes(x=pseudoH2,y=predmean), method="glm", formula=y~poly(x,2), color="grey",se=F)+
  scale_color_gradientn("",colours = rev(redblue))+
  ylab("Predictability year 1 (median stabilizing r LOO)")+
  ylab("Pseudo-h2 of frequency of accessions")+
  theme_minimal()+
  theme(legend.position = "none")

ggplot(vasurmean)+
  geom_point(aes(x=pseudoH2,y=predmean,color=bio1, size=survival2),
             shape=16)+
  stat_smooth(aes(x=pseudoH2,y=predmean), method="glm", formula=y~poly(x,2), color="grey",se=F)+
  scale_color_gradientn("",colours = rev(redblue))+
  ylab("Predictability year 1 (median stabilizing r LOO)")+
  ylab("Pseudo-h2 of frequency of accessions")+
  theme_minimal()+
  theme(legend.position = "none")
ggplot(vasurmean)+
  geom_point(aes(x=pseudoH2,y=predmean,color=bio1, size=survival5),
             shape=16)+
  stat_smooth(aes(x=pseudoH2,y=predmean), method="glm", formula=y~poly(x,2), color="grey",se=F)+
  scale_color_gradientn("",colours = rev(redblue))+
  ylab("Predictability year 1 (median stabilizing r LOO)")+
  ylab("Pseudo-h2 of frequency of accessions")+
  theme_minimal()

ggplot(vasurmean)+
  geom_point(aes(x=pseudoH2,y=predmean,color=bio1, size=pop),
             shape=16)+
  stat_smooth(aes(x=pseudoH2,y=predmean), method="glm", formula=y~poly(x,2), color="grey",se=F)+
  scale_color_gradientn("",colours = rev(redblue))+
  ylab("Predictability year 1 (median stabilizing r LOO)")+
  ylab("Pseudo-h2 of frequency of accessions")+
  theme_minimal()+
  theme(legend.position = "none")

# regresion heritability and survival
pdf("figs/fig-pseudo-heritability-h2-averagesurvival.pdf",
          width = 5,height = 5)

ggplot(vasurmean)+
  geom_point(aes(x=pseudoH2,y=survival3,color=bio1), size=4)+
  stat_smooth(aes(x=pseudoH2,y=survival3,color=bio1),method="glm" ,method.args = list(family = "binomial"))+
  scale_color_gradientn("",colours = rev(redblue))+
  ylab("Year 3 survival")+
  xlab("Pseudo-h2 of frequency of accessions")+
  theme_minimal()+
  theme(legend.position = "none")
ggplot(vasurmean)+
  geom_point(aes(x=pseudoH2,y=survival4,color=bio1), size=4)+
  stat_smooth(aes(x=pseudoH2,y=survival4,color=bio1),method="glm" ,method.args = list(family = "binomial"))+
  scale_color_gradientn("",colours = rev(redblue))+
  ylab("Year 2 survival")+
  xlab("Pseudo-h2 of frequency of accessions")+
  theme_minimal()+
  theme(legend.position = "none")
ggplot(vasurmean)+
  geom_point(aes(x=pseudoH2,y=survival5,color=bio1), size=4)+
  stat_smooth(aes(x=pseudoH2,y=survival5,color=bio1),method="glm" ,method.args = list(family = "binomial"))+
  scale_color_gradientn("",colours = rev(redblue))+
  ylab("Year 2 survival")+
  xlab("Pseudo-h2 of frequency of accessions")+
  theme_minimal()+
  theme(legend.position = "none")
dev.off()

#####*********************************************************************######
#####* Repeatability vs predictability ####
#####*
df_wide <- df %>%
  pivot_wider(
    id_cols = c(site, plot),
    names_from = source,
    values_from = c(sp_r)
  )
ggplot(df_wide)+
  geom_point(aes(x=climate_distance,y=stabilizing_prs, color=bio1))+
  scale_color_gradientn("",colours = rev(redblue))+
  geom_abline()
ggplot(df_wide)+
  geom_point(aes(x=plot_mean,y=stabilizing_prs, color=bio1))+
  scale_color_gradientn("",colours = rev(redblue))+
  geom_abline()



#####*********************************************************************######
#####* Genomic Offset ####

go<-read.csv("data-intermediate//genomic_offset_per_site_firt_gen_prob_squared_l1out_all_proba.csv",header=T)

head(go)

go$X5 %>% hist
go$X5 %>% summary
go$X4 %>% summary
apply(go, 2, mean) %>% hist
apply(go, 2, var) %>% hist
apply(go, 2, sum) %>% hist

metric=apply(go, 2, min)
metric=apply(go, 2, var) /apply(go, 2, mean)
metric=apply(go, 2, var)
metric=apply(go, 2, sum)
metric=apply(go, 2, mean)
metricvar<-apply(go, 2, var)
metricmean<-apply(go, 2, mean)
plot(metricvar,metricmean)



# make go metrics a table
govar<-data.frame(site=colnames(go),
                  govar = metricvar,
                  gomean = metricmean
                  ) %>%
  dplyr::mutate(site=gsub(site,pattern="X",replace= "",fixed = T))

# merge with survival and climate
govarsur<-merge(govar,sur, by="site") %>%
  merge(worldclim_sitesdata,by="site")

# check govar metrics with climate

meanandvariancego<-
qplot(data=govarsur, y=govar,x=gomean, color=bio1)+
  scale_color_gradientn("",colours = rev(redblue))+
  xlab("GO mean")+
  ylab("GO variance")+
  theme_minimal()
meanandvariancego
save_plot(filename = "figs/fig-genomicoffset-meanandvar-bio1.pdf",
          meanandvariancego,
          base_height = 2.5,base_width = 4)


qplot(data=govarsur, y=govar/gomean,x=gomean, color=bio1, size=X3_survival)+
  scale_color_gradientn("",colours = rev(redblue))+
  xlab("GO mean")+
  ylab("GO variance")+
  scale_y_log10()+
  scale_x_log10()+
  theme_minimal()
qplot(data=govarsur, y=govar,x=gomean, color=bio1, size=X4_survival)+
  scale_color_gradientn("",colours = rev(redblue))+
  xlab("GO mean")+
  ylab("GO variance")+
  scale_y_log10()+
  scale_x_log10()+
  theme_minimal()
qplot(data=govarsur, y=govar,x=gomean, color=bio1, size=X5_survival)+
  scale_color_gradientn("",colours = rev(redblue))+
  xlab("GO mean")+
  ylab("GO variance")+
  scale_y_log10()+
  scale_x_log10()+
  theme_minimal()

ggplot(data=govarsur)+
  # year 1
  geom_jitter(aes(y=govar,x=gomean, color=bio1, size=X1_survival), shape=1)+
  # year 2
  geom_jitter(aes(y=govar,x=gomean, color=bio1, size=X2_survival), shape=2)+
  # year 3
  geom_jitter(aes(y=govar,x=gomean, color=bio1, size=X3_survival), shape=16)+
  # year 4
  geom_jitter(aes(y=govar,x=gomean, color=bio1, size=X4_survival), shape=17)+
  # year 5
  geom_jitter(aes(y=govar,x=gomean, color=bio1, size=X5_survival), shape=18)+
  scale_color_gradientn("",colours = rev(redblue))+
  xlab("GO mean")+
  ylab("GO variance")+
  scale_y_log10()+
  scale_x_log10()+
  theme_minimal()+
  theme(legend.pos="none")

qplot(data=govarsur, y=govar,x=gomean, color=bio1, shape=ifelse(X5_survival==1,1,0) )+
  scale_color_gradientn("",colours = rev(redblue))+
  xlab("GO mean")+
  ylab("GO variance")+
  theme_minimal()

ggplot(data=govarsur)+
  geom_jitter(aes(y=govar,x=gomean, color=bio1, size=X2_flowerstotal) )+
  scale_color_gradientn("",colours = rev(redblue))+
  xlab("GO mean")+
  ylab("GO variance")+
  theme_minimal()


# fix -1
govarsur<-govarsur%>%
  mutate(
  X1_survival=ifelse(X1_survival>-1,X1_survival,NA),
  X2_survival=ifelse(X2_survival>-1,X2_survival,NA),
  X3_survival=ifelse(X3_survival>-1,X3_survival,NA),
  X4_survival=ifelse(X4_survival>-1,X4_survival,NA),
  X5_survival=ifelse(X5_survival>-1,X5_survival,NA)
  )
glm(X2_survival ~ govar , family = binomial, data = govarsur) %>% summary
glm(X3_survival ~ govar , family = binomial, data = govarsur) %>% summary
glm(X4_survival ~ govar , family = binomial, data = govarsur) %>% summary
glm(X5_survival ~ govar , family = binomial, data = govarsur) %>% summary
glm(X5_survival ~ gomean  * bio1 , family = binomial, data = govarsur) %>% summary

glm(X5_survival ~ bio1 , family = binomial, data = govarsur) %>% summary
glm(X5_survival ~ gomean*govar , family = binomial, data = govarsur) %>% summary

ggplot(govarsur)+
  geom_jitter(aes(y=X2_survival, x=rank(govar)),height = 0.1)
ggplot(govarsur)+
  geom_jitter(aes(y=X3_survival, x=rank(govar)),height = 0.1)
ggplot(govarsur)+
  geom_jitter(aes(y=X4_survival, x=rank(govar)),height = 0.1)
ggplot(govarsur)+
  geom_jitter(aes(y=X5_survival, x=rank(govar)),height = 0.1)

ggplot(govarsur)+
  geom_point(aes(y=bio1, x=govar, color=bio1))+
  scale_color_gradientn("",colours = rev(redblue))+
  theme_minimal()

# Get average survival of replicates
goex<-
  govarsur %>% # old
  group_by(site) %>%
  summarise(
    meansur1=mean(X2_survival),
    meansur2=mean(X2_survival),
    meansur3=mean(X3_survival),
    meansur4=mean(X4_survival),
    meansur5=mean(X5_survival),
    # meango=mean(govar),
    meango=mean(govar),
    meangom=mean(gomean),
    meanbio=mean(bio1)
  )
goex$proxysur<-
  goex[,c("meansur2","meansur3","meansur4","meansur5")] %>% apply(.,1,function(x)mean(x,na.rm=T))


cor.test(goex$proxysur, goex$meango, method='s')
cor.test(goex$proxysur, goex$meangom, method='s')
cor.test(goex$meansur5, goex$meanbio, method='s')

qplot(data=goex, y=meango,x=meangom, color=proxysur, size=3)+
  scale_color_gradientn("",colours = brewer.pal(9,"Greys"))+
  xlab("GO mean")+
  ylab("GO variance")+
  theme_minimal()

qplot(data=goex, y=meango,x=meangom, color=meansur4, size=3)+
  scale_color_gradientn("",colours = (redblue))+
  xlab("GO mean")+
  ylab("GO variance")+
  theme_minimal()

qplot(data=goex, y=meango,x=meangom, color=meansur2+meansur3+meansur4+meansur5, size=3)+
  scale_color_gradientn("",colours = (redblue))+
  xlab("GO mean")+
  ylab("GO variance")+
  theme_minimal()

#### Go total population
gonew<-read.csv("data-intermediate/go_sign_snps_from_founder_pop.csv")

goex<-
# govarsur %>% # old
  sur %>%  # new and merge
  merge(gonew, by="site") %>%
  merge(worldclim_sitesdata,by="site") %>%
  group_by(site) %>%
  summarise(
            meansur2=mean(X2_survival),
            meansur3=mean(X3_survival),
            meansur4=mean(X4_survival),
            meansur5=mean(X5_survival),
            # meango=mean(govar),
            meango=mean(go),
            meanbio=mean(bio1)
            )

ggplot(goex)+
  # year 2
  # geom_point(aes(y=meansur2, x=meango))+
  stat_smooth(aes(y=meansur2, x=meango), method='glm',color="grey90", se=F)+
  stat_smooth(aes(y=meansur3, x=meango), method='glm',color="grey70", se=F)+
  stat_smooth(aes(y=meansur4, x=meango), method='glm',color="grey40", se=F)+
  stat_smooth(aes(y=meansur5, x=meango), method='glm',color="black", se=F)+
  theme_minimal()
ggplot(goex)+
  # year 2
  # geom_point(aes(y=meansur2, x=meango))+
  stat_smooth(aes(y=meansur2, x=meanbio), method='glm',color="grey90", se=F)+
  stat_smooth(aes(y=meansur3, x=meanbio), method='glm',color="grey70", se=F)+
  stat_smooth(aes(y=meansur4, x=meanbio), method='glm',color="grey40", se=F)+
  stat_smooth(aes(y=meansur5, x=meanbio), method='glm',color="black", se=F)+
  theme_minimal()
ggplot(goex)+
  geom_point(aes(y=meansur5, x=meango))+
  stat_smooth(aes(y=meansur5, x=meango), method='glm')+
  theme_minimal()
ggplot(goex)+
  geom_point(aes(y=meansur5, x=meanbio))+
  stat_smooth(aes(y=meansur5, x=meanbio), method='glm')+
  theme_minimal()

cor.test(goex$meanbio, goex$meansur2, method='s')
cor.test(goex$meango, goex$meansur2, method='s')
cor.test(goex$meanbio, goex$meansur3, method='s')
cor.test(goex$meango, goex$meansur3, method='s')
cor.test(goex$meanbio, goex$meansur4, method='s')
cor.test(goex$meango, goex$meansur4, method='s')
cor.test(goex$meanbio, goex$meansur5, method='s')
cor.test(goex$meango, goex$meansur5, method='s')

ggplot(goex)+
  geom_point(aes(meanbio, meango))


lm(data=goex, meansur5 ~ meango ) %>% summary
lm(data=goex, meansur5 ~ meanbio ) %>% summary
lm(data=goex, meansur5 ~ meango * meanbio) %>% summary


#####*********************************************************************######
##### Genomic offset and predictability ####

mm<-
  merge(mersurclimate,goex, by='site')
head(mm)

ggplot(mm)+
  geom_point(aes(y = meango/meangom,  x=r2))
ggplot(mm)+
  geom_point(aes(y = meango/meangom,  x=sp_r))

ggplot(mm)+
  geom_point(aes(y = meango,  x=sp_r))

ggplot(mm)+
  geom_point(aes(y = meangom,  x=r2))

ggplot(mm,
       aes(x = meango,  y=r2, group=site,color=bio1)
       )+
  scale_color_gradientn("",colors=rev(redblue))+
  stat_summary()+
  stat_smooth(,method="glm")

ggplot(mm,
       aes(x = meangom,  y=r2, group=site,color=bio1)
)+
  scale_color_gradientn("",colors=rev(redblue))+
  stat_summary()+
  stat_smooth(,method="glm")


ggplot(mm,
       aes(x = meangom,  y=sp_r, group=site,color=bio1)
)+
  scale_color_gradientn("",colors=rev(redblue))+
  stat_summary()+
  stat_smooth(,method="glm")

ggplot(mm,
       aes(x = meango,  y=sp_r, group=site,color=bio1)
)+
  scale_color_gradientn("",colors=rev(redblue))+
  stat_summary()+
  stat_smooth(,method="glm")



ggplot(mm,
       aes(x = meangom,  y=r2, group=site,color=meansur2)
)+stat_summary()
ggplot(mm,
       aes(x = meangom,  y=r2, group=site,color=meansur3)
)+stat_summary()
ggplot(mm,
       aes(x = meangom,  y=r2, group=site,color=meansur4)
)+stat_summary()
ggplot(mm,
       aes(x = meangom,  y=r2, group=site,color=meansur5)
)+stat_summary()

