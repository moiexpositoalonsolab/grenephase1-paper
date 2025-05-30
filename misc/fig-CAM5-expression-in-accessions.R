################################################################################
### Goal
### Do study the differences in expression of CAM5 across accessions
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
setwd(myfolder)

################################################################################
# Read backbone fam 
# fam<-read.table("/2029g/2029g.fam",header=F)
# Load normalized expression data
load("data-external/arabidopsistranscriptomes/1001t/TG_data_20180606.Rdata")
load("data-external/arabidopsistranscriptomes/1001t/gene_infoV2.Rdata")

# remove batch MU, which was the original low quality
ind2use<- which(TG.meta$batch_comb != "MU") 
# subset log2 expression for all good individuals
d<-TG.genes$d_log2_batch[,ind2use] 
# Generate an expression matrix (individuals=rows, genes=cols)
d<-data.frame(TG.meta$index[ind2use], t(d))
colnames(d)<-c("id",gene_infoV2$Name)

##### Ectract
d2<-TG.trans$d_log2[,ind2use] 
# d2<-TG.trans$raw[,ind2use]
# d<-TG.genes$d_log2_batch[,ind2use] 
# Generate an expression matrix (individuals=rows, genes=cols)
d2<-data.frame(TG.meta$index[ind2use], t(d2))
colnames(d2)<-c("id",TG.trans$names)

################################################################################

# Try to find CAM5 in the expression data
genesubset<-d %>% dplyr::select(id,AT2G27030, AT2G22540)
transcriptsubset<-d2 %>% data.frame() %>%  dplyr::select(id,AT2G27030.1,AT2G27030.3)
head(transcriptsubset)


# Worldclim
worldclim_ecotypesdata<-read.csv("data-external/1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv")
head(worldclim_ecotypesdata)

# Merge
mer<-merge(genesubset,worldclim_ecotypesdata, by.x='id',by.y="ecotypeid")
mer2<-merge(transcriptsubset,worldclim_ecotypesdata, by.x='id',by.y="ecotypeid")

################################################################################
# Plot cam5 and temp
ggplot(mer)+
  geom_point(aes(y=AT2G27030, x=bio1))+
  stat_smooth(aes(y=AT2G27030, x=bio1),method='glm')
# # Plot svp and temp
# ggplot(mer)+
#   geom_point(aes(y=AT2G22540, x=bio1))+
#   stat_smooth(aes(y=AT2G22540, x=bio1),method='glm')

# cam5 splice forms
transcript1<-
ggplot(mer2)+
  geom_point(aes(y=AT2G27030.1, x=bio1))+
  stat_smooth(aes(y=AT2G27030.1, x=bio1),method='glm', color= "#2166AC")+
  labs(x="Annual temperature of accession (C)")
transcript1

transcript2<-
ggplot(mer2)+
  geom_point(aes(y=AT2G27030.3, x=bio1))+
  stat_smooth(aes(y=AT2G27030.3, x=bio1),method='glm', color="#A50F15")+
  labs(x="Annual temperature of accession (C)")
transcript2

#ratio
ggplot(mer2)+
  geom_point(aes(y=AT2G27030.3-AT2G27030.1, x=bio1))+
  stat_smooth(aes(y=AT2G27030.3-AT2G27030.1, x=bio1),method='glm', color="grey")+
  labs(x="Annual temperature (accession origin)")


transcriptcorrelations<- plot_grid(transcript1,transcript2)
save_plot(paste0("figs/fig_CAM5_expression.pdf"),transcriptcorrelations,
          base_width = 8,base_height = 5 )
save_plot(paste0("figs/fig_CAM5_expression.png"),transcriptcorrelations,
          base_width = 8,base_height = 5 )


sink("tables/test-CAM5-transcript-temperature-correlation.txt")

cor.test(mer$AT2G27030, mer$bio1, method='s')
cor.test(mer$AT2G27030, mer$bio1, method='p')

cor.test(mer2$AT2G27030.1, mer$bio1, method='p')
cor.test(mer2$AT2G27030.3, mer$bio1, method='p')

cor.test(mer2$AT2G27030.1, mer$bio1, method='s')
cor.test(mer2$AT2G27030.3, mer$bio1, method='s')

sink()


# # Check whether accessions with these expressions increase or decrease in grene-net
# # myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
# # setwd(myfolder)
# load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize.rda")
# e= ecotype_allyears_frequencies_long_raw_climate_popsize
# head(e)
# emer2<-merge(e,transcriptsubset,by="id")
# ressplice<-
#   emer2 %>% 
#   group_by(site,sites_bio1) %>% 
#   dplyr::filter(year==1) %>% 
#   dplyr::filter(sites_bio1< 18) %>% 
#   # dplyr::filter(site!=26, sites!=60, sites!=10) %>% # removing israel and davis
#   summarize(r3=cor(freq-startfreq, AT2G27030.3, method='p'),
#             r3p=cor.test(freq-startfreq, AT2G27030.3, method='p')$p.value,
#             r1=cor(freq-startfreq, AT2G27030.1, method='p'),
#             r1p=cor.test(freq-startfreq, AT2G27030.1, method='p')$p.value,
#   )
# cor.test(ressplice$r1,ressplice$sites_bio1)
# cor.test(ressplice$r3,ressplice$sites_bio1)
# cor.test(ressplice$r1,ressplice$sites_bio1, method='s')
# cor.test(ressplice$r3,ressplice$sites_bio1, method='s')
# 
# ggplot(ressplice)+
#   geom_point(aes(y=r3, x=sites_bio1, shape=r3p<0.01))+
#   geom_hline(yintercept = 0)
# ggplot(ressplice)+
#   geom_point(aes(y=r1, x=sites_bio1, shape=r1p<0.05))+
  geom_hline(yintercept = 0)

