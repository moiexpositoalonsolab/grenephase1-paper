################################################################################
### Goal
### Plot the Gaussian model lines
################################################################################

library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)

myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives/MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
## Raw data
file_edelta = 'data/delta_ecotype_freq.txt'
file_e0 = 'data/merged_ecotype_frequency.txt'
file_eid = 'data/greneNet_final_v1.0.fam'
file_climate_sites= "grene/data/worldclim_sitesdata.rda"
file_info_sites=  "grene/data/sites.clim.rda"
files_climate_ecotypes="grene/data/worldclim_ecotypesdata.rda"
files_info_ecotypes="grene/data/ecotypes_data.rda"

load(files_climate_ecotypes)
load(file_climate_sites)
load(files_info_ecotypes)

## Intermediate data (previous results Spring 2024)
Vsresults<-read.table("data-intermediate/stabilizing_selection_wmax_vs.txt",header = T)
Vsresults<-
  Vsresults %>% dplyr::rename(wmax=wmax_mean,Vs=vs_mean)

Vsresults_old<-Vsresults

# Merge
tmp= merge(ecotypes_data, by.x="ecotypeid", worldclim_ecotypesdata, by.y="ecotypeid")
dfeco= merge(Vsresults, by.x="ecotype", tmp, by.y="ecotypeid")

## New results Sep 2024
Vsresults=read.table("data-intermediate/LOCAL_ADAPTATION_ecotype_Wmax_Vs.txt",header=T)
head(Vsresults)
Vsresults<-Vsresults %>% dplyr::rename(wmax=Wmax,Vs=Vs)



################################################################################
# Plot a model example

lineprojection<-
  expand.grid(unique(dfeco$ecotype), seq(0,25,0.5)) %>% 
  rename(id=Var1, bio1transplant=Var2) %>% 
  merge(.,dfeco, by.x="id",by.y="ecotype") %>% 
  # mutate(Westimated=wmax-(1/Vs)*(bio1-bio1transplant)^2) # Gaussian model. incorrect wmax already in the real space
  mutate(Westimated=wmax * exp(( -(bio1-bio1transplant)^2) / Vs)) 
  # dfeco$vulnerability <- dfeco$wmax * exp(-(dfeco$temperaturechange) / dfeco$Vs)

redgreen<-c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) )
redgreen<-c( "#A50F15" , "#DE2D26" , "#FB6A4A" , "#FCAE91" , "#FEE5D9" , "#EDF8E9" , "#BAE4B3" , "#74C476" , "#31A354" , "#006D2C")
redgreen<-c( "#A50F15" , "#DE2D26" , "#FB6A4A" , "#FCAE91" , "#FEE5D9" , "#EDF8E9" , "#BAE4B3" , "#74C476" , "#74C476" , "#31A354", "#31A354" , "#006D2C","#006D2C")

# [1] "#F7FCF5" "#E5F5E0" "#C7E9C0" "#A1D99B" "#74C476" "#41AB5D" "#238B45" "#006D2C" "#00441B"
ggplot(lineprojection)+
  geom_line(aes(y=Westimated,x=bio1-bio1transplant,color=1/Vs, group=id))+
  scale_color_gradientn("1/Vs",colors=redgreen)+
  xlab("Temperature (C)")+
  ylab("W modeled")+
  theme_minimal()
modellines_for_talk<-
ggplot(lineprojection)+
  geom_line(aes(y=Westimated,x=bio1-bio1transplant,color=1/Vs, group=id))+
  scale_color_gradientn("1/Vs",colors=brewer.pal(9,"Greys")[-c(1:2) ])+
  xlab("Temperature (C)")+
  ylab("W modeled")+
  theme_minimal()
save_plot("figs/fig-Vs-wmax-modellines-greys.pdf",modellines_for_talk, base_height = 4,base_width = 5)
save_plot("figs/fig-Vs-wmax-modellines-greys.png",modellines_for_talk, base_height = 4,base_width = 5)

fig_stabilizing_modellines<-
  ggplot(lineprojection)+
    geom_line(aes(y=Westimated,x=bio1transplant,color=1/Vs, group=id))+
    scale_color_gradientn("1/Vs",colors=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ))+
    xlab("Temperature (C)")+
    ylab("W modeled")+
    theme_minimal()
fig_stabilizing_modellines

save_plot("figs/fig-Vs-wmax-modellines.pdf",fig_stabilizing_modellines, base_height = 4,base_width = 5)
save_plot("figs/fig-Vs-wmax-modellines.png",fig_stabilizing_modellines, base_height = 4,base_width = 5)

fig_stabilizing_modellines_naturalscale<-
  ggplot(lineprojection)+
  geom_line(aes(y=exp(Westimated),x=bio1transplant,color=1/Vs, group=id))+
  scale_color_gradientn("1/Vs",colors=c( rev(brewer.pal(5,"Reds")[1:5]),(brewer.pal(5,"Greens")[1:5]) ))+
  xlab("Temperature (C)")+
  ylab("W modeled")+
  theme_minimal()
fig_stabilizing_modellines_naturalscale

save_plot("figs/fig-Vs-wmax-modellines-naturalscale.pdf",fig_stabilizing_modellines_naturalscale, base_height = 5,base_width = 6)
save_plot("figs/fig-Vs-wmax-modellines-naturalscale.png",fig_stabilizing_modellines_naturalscale, base_height = 5,base_width = 6)
