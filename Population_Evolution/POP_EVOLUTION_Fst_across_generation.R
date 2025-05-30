rm(list=ls())

library(ggplot2)
library(patchwork)
library(dplyr)

## RAPID EVOLUTION

## This script plot the Fst comparison across generation and sites

prefix <- "POP_EVOLUTION_"
PATH<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE Moi
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

fst_generation <- read.delim("data-intermediate/fst_founder_vs_samples.csv",sep=',')
fst_pop <- read.delim("data-intermediate/fst_matrix.csv",sep=",")

fst_generation %>%
  group_by(Generation) %>%
  summarise(mean_fst=median(FST))


ggplot(fst_generation,aes(x=Generation,y=FST))+
  geom_boxplot(outlier.color = "red",outlier.shape = 8,outlier.size = 3)+
  geom_jitter(shape=16, position=position_jitter(0.3),fill="grey70",cex=5,alpha=0.2)+
  ylab("Fst between founder \n and evolved populations")+
  theme_classic(base_size = 18)

ggsave(file=paste0("figs/",prefix,"Fst_across_generation.pdf"),device = "pdf",width = 8,height = 6)

sink("tables/test_fst_average_across_generations.txt")
fst_generation %>% 
  group_by(Generation) %>% 
  summarise(mean=median(FST),
            IQRlow=quantile(FST,0.025),
            IQRup=quantile(FST,0.975)
            )

fst_generation %>% 
  group_by(Generation) %>% 
  summarise(mean=median(FST),
            IQRlow=quantile(FST,0.25),
            IQRup=quantile(FST,0.75)
  )

sink()


# Check whether Fst also decreases if we only focus on popualtions that have
fstcombine<-
  fst_pop %>% 
  dplyr::select(sample,founder) %>% 
  dplyr::mutate(sampleid=sample) %>% 
  dplyr::rename(FST=founder) %>% 
  tidyr::separate(sampleid,into = c("site","year","rep"))
head(fstcombine)

sink("tables/test_fst_average_across_generations-CORRECTED_ONLY_SITES_WITH_3_GEN.txt")

# Without removing lost sites
fstcombine %>% 
  group_by(year) %>% 
  summarise(
            mean=mean(FST),
            median=median(FST),
            IQRlow=quantile(FST,0.025),
            IQRup=quantile(FST,0.975)
  )

# Removing the sites that do not have geneartion 3
fstcombine %>%
  dplyr::group_by(site) %>% 
  dplyr::mutate(maxyear=max(year)) %>% 
  ungroup %>% 
  dplyr::filter(maxyear==3) %>% 
  dplyr::group_by(year) %>% 
  summarise(
    mean=mean(FST),
    median=median(FST),
    IQRlow=quantile(FST,0.025),
    IQRup=quantile(FST,0.975)
  )  


tmp<-
  fstcombine %>%
  dplyr::group_by(site) %>% 
  dplyr::mutate(maxyear=max(year)) %>% 
  ungroup %>% 
  dplyr::filter(maxyear==3)

t.test( tmp$FST[tmp$year==3])
t.test( tmp$FST[tmp$year==2])


sink()
