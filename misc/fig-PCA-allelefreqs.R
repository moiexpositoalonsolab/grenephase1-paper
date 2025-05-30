################################################################################
### Goal
### Get a similar plot as the Diffusion map for those non-believers on
### things that are not a PCA for decomposition of variable
### Using allele frequencies instead ecotypes
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
# Load files
# Load ecotypes and their ID
eco0<-read.table("data/founder_ecotype_frequency.txt", header=F) %>% 
  dplyr::rename(id=V1,e0=V2)
ecofreq<-read.table("data/merged_ecotype_frequency.txt", header=T) 
ecodelta<-read.table("data/delta_ecotype_freq.txt", header=T)

a0<-read.table("data/average_seedmix_p0_LDpruned.txt", header=F) %>% 
  dplyr::rename(Chr=V1,Pos=V2, freq=V3)
afreq<-read.table("data/merged_hapFIRE_allele_frequency_LDpruned.txt", header=T)

load("grene/data/samples_data.rda")
samples_data<-
  samples_data %>% 
  data.frame %>% 
  dplyr::mutate(year=as.numeric(year)-2017)

################################################################################
# transpose for analyses, and parse
ta<-t(afreq)

# recode info
tam<-
  ta %>% 
    data.frame() %>% 
    dplyr::mutate(sampleid=rownames(ta)) %>% 
    mutate(sampleid=gsub( "X", "", sampleid)) %>% 
    mutate(sample=sampleid) %>% 
    separate(sampleid, c("site","year","rep"), "_")
  
# Merge, maybe not necessary
# tam<-
#   tam %>% 
#   merge(.,by.x=c("site","year","rep"),
#        samples_data, by.y=c("site","year","plot")
#        )
#   


################################################################################
#### PCA analysis
# run pca
pcamod<-prcomp( ta)

# plot(pcamod)
pdat<-pcamod$x
dim(pdat)
pdat<-data.frame(pdat)

# calculate variance explained
percentageesplained<- pcamod$sdev^2 / sum(pcamod$sdev^2)
percentageesplained<-round(percentageesplained*100,1)
head(percentageesplained)
percentageesplainedlabels<- paste0("PC",1:length(percentageesplained)," (", percentageesplained,"%)")
percentageesplainedlabels %>% head

# Parse columns to include site information
load("grene/data/locations_data.rda")
load("grene/data/worldclim_sitesdata.rda")
sites_simple_names<-read.csv("grene/data/sites_simple_names.csv")


pdat$sampleid<-rownames(pdat)
pdat$sample= pdat$sampleid # backup

p<-
pdat %>% 
  # head %>% 
  mutate(sampleid=gsub( "X", "", sampleid)) %>% 
  mutate(sample=sampleid) %>% 
  separate(sampleid, c("site","year","rep"), "_") %>% 
  merge(
      dplyr::select(locations_data,longitude,latitude,altitude, site), 
      by='site'
      ) %>% 
  merge(
    .,
    worldclim_sitesdata, by="site") %>% 
  merge(
    .,
    sites_simple_names, by="site") %>% 
  mutate(name=paste(city, country))

# 

p1<-
ggplot(p) +
  geom_point(aes(x=PC1, y=PC2, color=as.numeric(year)))+
  scale_color_continuous("year",type = 'viridis')+
  # scale_color_gradientn("year",colors=brewer.pal(6,'Greys')[c(2,4,6)])+
  xlab(percentageesplainedlabels[1])+
  ylab(percentageesplainedlabels[2])+
  theme_minimal()
p1
p2<-
ggplot(p) +
  geom_point(aes(x=PC1, y=PC2, color=bio1))+
  scale_color_gradientn("Temperature (C)",colors = brewer.pal(5,"Reds"))+
  xlab(percentageesplainedlabels[1])+
  ylab(percentageesplainedlabels[2])+
  theme_minimal() 
p3<-
ggplot(p) +
  geom_point(aes(x=PC1, y=PC2, color=bio12))+
  xlab(percentageesplainedlabels[1])+
  ylab(percentageesplainedlabels[2])+
  scale_color_gradientn("Precipitation (mm)",colors = brewer.pal(5,"Blues"))+
  theme_minimal()
p4<-
ggplot(p) +
  geom_point(aes(x=PC1, y=PC2, color=factor(name)))+
  xlab(percentageesplainedlabels[1])+
  ylab(percentageesplainedlabels[2])+
  theme_minimal()


save_plot("figs/fig-PCA-allele-frequencies-color-by-year.pdf",
          p1,base_height = 5, base_width = 5)
save_plot("figs/fig-PCA-allele-frequencies-color-by-year.png",
          p1,base_height = 5, base_width = 5)
save_plot("figs/fig-PCA-allele-frequencies-color-by-bio1.pdf",
          p2,base_height = 5, base_width = 6)
save_plot("figs/fig-PCA-allele-frequencies-color-by-bio1.png",
          p2,base_height = 5, base_width = 6)
save_plot("figs/fig-PCA-allele-frequencies-color-by-bio12.pdf",
          p3,base_height = 5, base_width = 6)
save_plot("figs/fig-PCA-allele-frequencies-color-by-bio12.png",
          p3,base_height = 5, base_width = 6)
save_plot("figs/fig-PCA-allele-frequencies-color-by-site.pdf",
          p4,base_height = 5, base_width =  9)
save_plot("figs/fig-PCA-allele-frequencies-color-by-site.png",
          p4,base_height = 5, base_width = 9)

################################################################################
## Training onlast generation and project back

# Step 2: Run PCA on the subset
tam %>% 
  dplyr::filter(year==3) %>%  # only train in last generation
  dplyr::select(starts_with("X")) %>% 
prcomp(.) -> pcamod_subset

# Extract percentage explained
percentageesplained<- pcamod_subset$sdev^2 / sum(pcamod_subset$sdev^2)
percentageesplained<-round(percentageesplained*100,1)
head(percentageesplained)
percentageesplainedlabels<- paste0("PC",1:length(percentageesplained)," (", percentageesplained,"%)")
percentageesplainedlabels %>% head

# Step 3: Project the original afreq data onto the principal components from the subset PCA
# Extract the rotation (loadings) matrix from the subset PCA
rotation_matrix <- pcamod_subset$rotation

# Project the original afreq data onto the principal components
projected_data <- ta %*% rotation_matrix
projected_data_founder <- a0$freq %*% rotation_matrix

# Convert the projected data to a data frame for further use
projected_data_df <- data.frame(projected_data)


pp<-
  projected_data_df %>% 
  mutate(sampleid=rownames(projected_data_df)) %>% 
  mutate(sampleid=gsub( "X", "", sampleid)) %>% 
  mutate(sample=sampleid) %>% 
  separate(sampleid, c("site","year","rep"), "_") %>% 
  merge(
    dplyr::select(locations_data,longitude,latitude,altitude, site), 
    by='site'
  ) %>% 
  merge(
    .,
    worldclim_sitesdata, by="site") %>% 
  merge(
    .,
    sites_simple_names, by="site") %>% 
  mutate(name=paste(city, country))

# Do the plotting
redblue<-brewer.pal(9,"RdBu")
redblue<-c( "#B2182B" , "#D6604D" , "#F4A582" , "#FDDBC7" , "#F7F7F7" , "#D1E5F0" , "#92C5DE" , "#4393C3" , "#2166AC")
redblue<-c( "#B2182B", "#B2182B" , "#D6604D" ,"#D6604D", "#F4A582" , "#F7F7F7" , "#D1E5F0" , "#92C5DE" , "#4393C3" , "#2166AC")

# First everything to visualize
pp %>%  
  ggplot(.)+ # simple one with everythign
  geom_point(aes(x=PC1, y=PC2, color=bio1))+
  scale_color_gradientn("Temperature (C)",colors = brewer.pal(5,"Reds"))+
    xlab(percentageesplainedlabels[1])+
    ylab(percentageesplainedlabels[2])+
    xlim(c(-10,30))+ylim(c(-10,45))+
    theme_classic()+
    theme(axis.line.x = element_blank())+
    theme(axis.line.y = element_blank())

# Do the plots one generation at a time
p0bio1<-
  ggplot()+
  geom_point(data=projected_data_founder,aes(x=PC1, y=PC2), color="black", size=3) +
  xlab(percentageesplainedlabels[1])+
  ylab(percentageesplainedlabels[2])+
  xlim(c(-10,30))+ylim(c(-10,45))+
  theme_classic()+
  theme(axis.line.x = element_blank())+
  theme(axis.line.y = element_blank())
p0bio1

p1bio1<-
pp %>% 
  dplyr::filter(year==1) %>%
ggplot(.) +
  geom_segment(aes(xend=PC1, yend=PC2,x=projected_data_founder$PC1, y=projected_data_founder$PC2 , color=bio1))+
  geom_point(aes(x=PC1, y=PC2), color="white", size=2.2)+
  geom_point(aes(x=PC1, y=PC2), color="black", size=1.8)+
  geom_point(aes(x=PC1, y=PC2, color=bio1))+
  # scale_color_gradientn("Temperature (C)",colors = brewer.pal(5,"Reds"))+
  scale_color_gradientn("Temperature (C)",colors = rev(redblue))+
  geom_point(data=projected_data_founder,aes(x=PC1, y=PC2), color="black", size=3) +# add the founder
  xlab(percentageesplainedlabels[1])+
  ylab(percentageesplainedlabels[2])+
  xlim(c(-10,30))+ylim(c(-10,45))+
  theme_classic()+
  theme(axis.line.x = element_blank())+
  theme(axis.line.y = element_blank())
p1bio1

p2bio1<-
pp %>% 
  dplyr::filter(year==2) %>% 
  ggplot(.) +
  geom_segment(aes(xend=PC1, yend=PC2,x=projected_data_founder$PC1, y=projected_data_founder$PC2 , color=bio1))+
  geom_point(aes(x=PC1, y=PC2), color="white", size=2.2)+
  geom_point(aes(x=PC1, y=PC2), color="black", size=1.8)+
  geom_point(aes(x=PC1, y=PC2, color=bio1))+
  # scale_color_gradientn("Temperature (C)",colors = brewer.pal(5,"Reds"))+
  scale_color_gradientn("Temperature (C)",colors = rev(redblue))+
  geom_point(data=projected_data_founder,aes(x=PC1, y=PC2), color="black", size=3) +# add the founder
  xlab(percentageesplainedlabels[1])+
  ylab(percentageesplainedlabels[2])+
  xlim(c(-10,30))+ylim(c(-10,45))+
  theme_classic()+
  theme(axis.line.x = element_blank())+
  theme(axis.line.y = element_blank())
p2bio1

p3bio1<-
pp %>% 
  dplyr::filter(year==3) %>% 
  ggplot(.) +
  geom_segment(aes(xend=PC1, yend=PC2,x=projected_data_founder$PC1, y=projected_data_founder$PC2 , color=bio1))+
  geom_point(aes(x=PC1, y=PC2), color="white", size=2.2)+
  geom_point(aes(x=PC1, y=PC2), color="black", size=1.8)+
  geom_point(aes(x=PC1, y=PC2, color=bio1))+
  # scale_color_gradientn("Temperature (C)",colors = brewer.pal(5,"Reds"))+
  scale_color_gradientn("Temperature (C)",colors = rev(redblue))+
  geom_point(data=projected_data_founder,aes(x=PC1, y=PC2), color="black", size=3) +# add the founder
  xlab(percentageesplainedlabels[1])+
  ylab(percentageesplainedlabels[2])+
  xlim(c(-10,30))+ylim(c(-10,45))+
  theme_classic()+
  theme(axis.line.x = element_blank())+
  theme(axis.line.y = element_blank())
p3bio1

cowplot::save_plot("figs/fig-PCA-allele-frequencies-color-by-bio1-generation0.pdf",
          p0bio1,base_height = 3, base_width = 5)
cowplot::save_plot("figs/fig-PCA-allele-frequencies-color-by-bio1-generation0.png",
          p0bio1,base_height = 3, base_width = 5)
cowplot::save_plot("figs/fig-PCA-allele-frequencies-color-by-bio1-generation1.pdf",
          p1bio1,base_height = 3, base_width = 5)
cowplot::save_plot("figs/fig-PCA-allele-frequencies-color-by-bio1-generation1.png",
          p1bio1,base_height = 3, base_width = 5)
cowplot::save_plot("figs/fig-PCA-allele-frequencies-color-by-bio1-generation2.pdf",
          p2bio1,base_height = 3, base_width = 5)
cowplot::save_plot("figs/fig-PCA-allele-frequencies-color-by-bio1-generation2.png",
          p2bio1,base_height = 3, base_width = 5)
cowplot::save_plot("figs/fig-PCA-allele-frequencies-color-by-bio1-generation3.pdf",
          p3bio1,base_height = 3, base_width = 5)
cowplot::save_plot("figs/fig-PCA-allele-frequencies-color-by-bio1-generation3.png",
          p3bio1,base_height = 3, base_width = 5)



#### Same but bio18 or bio12

p1bio18<-
  pp %>% 
  dplyr::filter(year==1) %>%
  ggplot(.) +
  geom_point(aes(x=PC1, y=PC2), color="white", size=2.2)+
  geom_point(aes(x=PC1, y=PC2), color="black", size=1.8)+
  geom_point(aes(x=PC1, y=PC2, color=bio18))+
  scale_color_gradientn("Summer precipitation (mm)",colors = brewer.pal(5,"Blues"))+
  geom_point(data=projected_data_founder,aes(x=PC1, y=PC2), color="black", size=3) +# add the founder
  xlab(percentageesplainedlabels[1])+
  ylab(percentageesplainedlabels[2])+
  xlim(c(-10,30))+ylim(c(-10,45))+
  theme_classic()+
  theme(axis.line.x = element_blank())+
  theme(axis.line.y = element_blank())
p1bio18

p2bio18<-
  pp %>% 
  dplyr::filter(year==2) %>% 
  ggplot(.) +
  geom_point(aes(x=PC1, y=PC2), color="white", size=2.2)+
  geom_point(aes(x=PC1, y=PC2), color="black", size=1.8)+
  geom_point(aes(x=PC1, y=PC2, color=bio18))+
  scale_color_gradientn("Summer precipitation (mm)",colors = brewer.pal(5,"Blues"))+
  geom_point(data=projected_data_founder,aes(x=PC1, y=PC2), color="black", size=3) +# add the founder
  xlab(percentageesplainedlabels[1])+
  ylab(percentageesplainedlabels[2])+
  xlim(c(-10,30))+ylim(c(-10,45))+
  theme_classic()+
  theme(axis.line.x = element_blank())+
  theme(axis.line.y = element_blank())
p2bio18

p3bio18<-
  pp %>% 
  dplyr::filter(year==3) %>% 
  ggplot(.) +
  geom_point(aes(x=PC1, y=PC2), color="white", size=2.2)+
  geom_point(aes(x=PC1, y=PC2), color="black", size=1.8)+
  geom_point(aes(x=PC1, y=PC2, color=bio18))+
  scale_color_gradientn("Summer precipitation (mm)",colors = brewer.pal(5,"Blues"))+
  geom_point(data=projected_data_founder,aes(x=PC1, y=PC2), color="black", size=3) +# add the founder
  xlab(percentageesplainedlabels[1])+
  ylab(percentageesplainedlabels[2])+
  xlim(c(-10,30))+ylim(c(-10,45))+
  theme_classic()+
  theme(axis.line.x = element_blank())+
  theme(axis.line.y = element_blank())
p3bio18

save_plot("figs/fig-PCA-allele-frequencies-color-by-bio18-generation1.pdf",
          p1bio18,base_height = 3, base_width = 5)
save_plot("figs/fig-PCA-allele-frequencies-color-by-bio18-generation1.png",
          p1bio18,base_height = 3, base_width = 5)
save_plot("figs/fig-PCA-allele-frequencies-color-by-bio18-generation2.pdf",
          p2bio18,base_height = 3, base_width = 5)
save_plot("figs/fig-PCA-allele-frequencies-color-by-bio18-generation2.png",
          p2bio18,base_height = 3, base_width = 5)
save_plot("figs/fig-PCA-allele-frequencies-color-by-bio18-generation3.pdf",
          p3bio18,base_height = 3, base_width = 5)
save_plot("figs/fig-PCA-allele-frequencies-color-by-bio18-generation3.png",
          p3bio18,base_height = 3, base_width = 5)

