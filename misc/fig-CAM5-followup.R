################################################################################
### Goal
### Follow up CAM5 temporal trends and across sites
### 

################################################################################

# library(grene)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
theme_set(theme_minimal())

################################################################################

myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
# setwd("~/Shareddrives/MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/data-external//arabidopsistranscriptomes/")

################################################################################
##### Read cam5 #####
cam<-read_merged_allele_subset(
  subsetfile = "data-intermediate/merged_hapFIRE_allele_frequency-AT2G27030-CAM5.csv",
  headfile="data-intermediate/merged_hapFIRE_allele_frequency-header.csv",
  snpfile="data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-AT2G27030-CAM5.csv"
)

head(cam)
dim(cam)

# get LMs just for corroboration
camres<-
  cam %>% 
  dplyr::group_by(snp) %>% 
  dplyr::summarize(r=cor.test(freq-startfreq,bio1)$estimate,
                   p=cor.test(freq-startfreq,bio1)$p.value) %>% 
  merge(.,cam,by="snp")
# minimanhattan  
camres %>% 
  ggplot(.)+
  geom_point(aes(y=-log10(p), x=pos))+
  geom_vline( xintercept = 11531967)+
  geom_vline( xintercept = 11534358)+
  xlim(c(11531967-500,11534358+500))+
  theme_minimal()
# create top snps and neutral
camrestop<-
  camres %>% 
  dplyr::filter(-log10(p)>14) %>% 
  dplyr::filter(startfreq<0.5) 
camresmax<-
  camres %>% 
  dplyr::filter(-log10(p)  == max(-log10(p)))
camresneutral<-
  camres %>% 
  dplyr::mutate(roundfreq=round(freq*1000)) %>% 
  dplyr::filter(roundfreq %in% round(camrestop$freq*1000)) %>% 
  dplyr::filter(!(snp%in%camrestop$snp), -log10(p)<5) %>% 
  dplyr::filter(snp %in% sample(unique(snp), 10))


#plots
topsnpplot_cam5<-
  camresmax %>% 
  dplyr::filter(year==1) %>% 
  ggplot(.)+
  geom_point(aes((freq-startfreq),x=bio1),alpha=0.5,color="#75C376")+
  stat_summary(aes((freq-startfreq),x=bio1), size=1, alpha=1,color="#75C376")+
  stat_smooth(aes((freq-startfreq),x=bio1), method='glm', color="#75C376")+
  geom_hline(yintercept = 0,lty="dotted")+
  # scale_y_log10()+
  ylab("p1-p0")+
  xlab("Annual temperature (C)")+
  ylim(c(-0.06,+0.14))+
  theme_minimal()
topsnpplot_cam5
save_plot(plot=topsnpplot_cam5, file="figs/fig-allele_climate_lfmm_adaptive-AT2G27030-CAM5.pdf",base_height = 3.5,base_width = 3.5)
save_plot(plot=topsnpplot_cam5, file="figs/fig-allele_climate_lfmm_adaptive-AT2G27030-CAM5.png",base_height = 3.5,base_width = 3.5)
# #plots ofneutral SNps
neutralsnpplot_cam5<-
  camresneutral %>% 
  dplyr::filter(year==1) %>% 
  ggplot(.)+
  geom_point(aes((freq-startfreq),x=bio1),alpha=0.5,color="grey")+
  stat_summary(aes((freq-startfreq),x=bio1), size=1, alpha=1,color="grey")+
  stat_smooth(aes((freq-startfreq),x=bio1), method='glm', color="grey")+
  geom_hline(yintercept = 0,lty="dotted")+
  ylab("p1-p0")+
  xlab("Annual temperature (C)")+
  ylim(c(-0.06,+0.14))+
  theme_minimal()
neutralsnpplot_cam5
save_plot(plot=neutralsnpplot_cam5, file="figs/fig-allele_climate_lfmm_adaptive-AT2G27030-CAM5-neutral.pdf",base_height = 3.5,base_width = 3.5)
save_plot(plot=neutralsnpplot_cam5, file="figs/fig-allele_climate_lfmm_adaptive-AT2G27030-CAM5-neutral.png",base_height = 3.5,base_width = 3.5)

# Do the temporal analysis

cam5snpstemporal<-
  camresmax %>% 
  # dplyr::filter(bio1>15 | bio1<8 ) %>% # only look at extremes
  dplyr::mutate(warmvscold= bio1>9.5) %>% 
  # dplyr::filter(year>0) %>% 
  ggplot(.)+
  geom_jitter(aes(y=freq-startfreq, x=year), alpha=0.2,width = 0.025)+
  stat_summary(aes(y=freq-startfreq, x=year, group=site))+
  stat_smooth(aes(y=freq-startfreq, x=year, group=warmvscold, color=warmvscold), method='glm')+
  scale_color_manual(values = c( "#4E99C5", "#C4403C"))+
  geom_hline(lty='dotted',yintercept = 0)+
  # geom_line(aes(y=freq, x=year, group=site))+
  facet_wrap(~warmvscold, ncol=1)+  
  theme_minimal()+
  ylab("p1-p0")+
  xlab("Year")+
  coord_cartesian(y=c(-0.15,+0.25))+
  theme(legend.position = "none")
cam5snpstemporal

save_plot(plot=cam5snpstemporal, file="figs/fig-allele_climate_lfmm_adaptive-AT2G27030-CAM5-warm-vs-cold-trajectories.pdf",base_height = 7,base_width = 3.5)
save_plot(plot=cam5snpstemporal, file="figs/fig-allele_climate_lfmm_adaptive-AT2G27030-CAM5-warm-vs-cold-trajectories.png",base_height = 7,base_width = 3.5)


camresmax %>% 
  dplyr::filter(bio1>15 | bio1<8 ) %>% # only look at extremes
  dplyr::mutate(warmvscold= bio1>10) %>% 
  group_by(warmvscold) %>% 
  summarise(b=coefficients(summary(lm((freq-startfreq) ~year)))[2,1],
            p=coefficients(summary(lm((freq-startfreq) ~year)))[2,4]
            )
camresmax %>% 
  dplyr::filter(bio1>15 | bio1<8 ) %>% 
  dplyr::select(bio1, site) %>% 
  arrange(site)

camresmax %>% 
  dplyr::filter(bio1>15 | bio1<8 ) %>% # only look at extremes
  dplyr::mutate(warmvscold= bio1>10) %>% 
  group_by(site) %>% 
  summarise(b=coefficients(summary(lm((freq-startfreq) ~year)))[2,1],
            p=coefficients(summary(lm((freq-startfreq) ~year)))[2,4]
  ) %>% 
  arrange(b)

  # pull(site) %>% unique

### Add the gene models below
# Cam5 gene models
tairraw<-read.table("data-external/TAIR10_GFF3_genes_transposons.gff")
tairsub<-
  tairraw %>% 
  dplyr::filter(grepl("AT2G27030",V9))
tairsub
exon_data <- tairsub %>%
  filter(V3 == "exon") %>%
  mutate(start = as.numeric(V4),
         end = as.numeric(V5),
         isoform = V9)  # Extract exon information
camres %>% 
  ggplot(.)+
  geom_point(aes(y=-log10(p), x=pos))+
  geom_vline( xintercept = 11531967)+
  geom_vline( xintercept = 11534358)+
  xlim(c(11531967-500,11534358+500))+
  geom_rect(data=exon_data,  # adding the genes
            aes(xmin = start, xmax = end,
                ymin = -as.numeric(as.factor(isoform)) - 0.4, ymax = -as.numeric(as.factor(isoform)) -0.1
            ), fill = "grey") +
  theme_minimal()


# 2_11533937
