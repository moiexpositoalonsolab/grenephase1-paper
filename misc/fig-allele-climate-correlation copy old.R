################################################################################
### Goal
### Estimate allele climate associations simple GLM

################################################################################
#### If transformed to collab can use the stuff below
# from google.colab import drive
# drive.mount('/content/drive')
# %load_ext rpy2.ipython
# %%R

################################################################################
#### Libraries
library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(patchwork)

################################################################################
#### Location if run locally
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
# myfolder<-"~/grenephase1-analyses"

myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################

# Specify the file path
file_edelta = '/data/delta_ecotype_freq.txt'
file_e0 = '/data/merged_ecotype_frequency.txt'
file_eid = '/data/greneNet_final_v1.0.fam'
file_climate_sites="/grene/data/worldclim_sitesdata.rda"
file_sites_info="/grene/data/sites.clim.rda"

file_pdelta = '/data/merged_hapFIRE_allele_frequency_LDpruned.txt'
file_p0 = '/data/average_seedmix_p0_LDpruned.txt'
# 
# file_path_p0 <- paste0(myfolder, file_p0)
# file_path_pdelta <- paste0(myfolder, file_pdelta)
# file_path_edelta <- paste0(myfolder, file_edelta)
# file_path_e0 <- paste0(myfolder, file_e0)
# file_path_eid <- paste0(myfolder, file_eid)
# file_path_climate_sites=paste0(myfolder, file_climate_sites)
# file_path_sites_info=paste0(myfolder, file_sites_info)

################################################################################

# Load frequencies 
pdelta<-read.table("data/merged_hapFIRE_delta_p_LDpruned.txt", header=T)
p0<-read.table("data/average_seedmix_p0_LDpruned.txt", header=F)
p0$snp<-paste0(p0$V1,"_",p0$V2)
pdelta$snp<-p0$snp
p0<-p0 %>% rename(p0=V3)

# check dataset
head(pdelta)
dim(pdelta)
head(p0)
dim(p0)

################################################################################
# Make long to parse the column names
plong<-
  pdelta %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to=c("site","year","rep"),
    values_to = c("freq"),
    names_pattern = "X?(.*)_(.)_(.*)"
  ) %>%
  dplyr::select(site, year, rep, freq, snp)

df<- plong



# Extract the terminal year of frequency change
df<-
  df %>%
  group_by(site, rep, snp) %>%
  mutate(year=as.numeric(year)) %>%
  reframe(max_year = if(all(is.na(freq))) NA else max(year[!is.na(freq)], na.rm = TRUE),
          maxfreq = freq[year==max_year])

# Add starting frequency
df<- merge(
  df,
  select(p0,p0,snp)
)
# Create change in frequency
df<- df %>% 
  mutate(deltafreq= maxfreq - p0)

# # Spread the data to wide format
# df_wide <-
#   df %>%
#   mutate(site_rep=paste0(site,"_",rep)) %>% 
#   pivot_wider(names_from = id, values_from = maxfreq, id_cols = site_rep)
# head(df_wide)

################################################################################
#### Add climate

##############
# Load climates
load(file_path_climate_sites)
load(file_path_sites_info)

worldclim_sitesdata<-worldclim_sitesdata %>% data.frame %>% mutate_all(as.numeric)
sites.clim<-sites.clim %>% data.frame %>%  mutate_all(as.numeric)


# Merge
dfclim<-
  merge(df, worldclim_sitesdata,by="site")

################################################################################
#### Simple GLM

testdelta<-dfclim$deltafreq
testyear<-dfclim$max_year
testbio<-dfclim$bio1

glmfunctionb<-function(testdelta,testyear,testbio){
  lmmod<-lm(testdelta ~ testyear + testbio)
  p=coefficients(summary(lmmod))[3,4]
  b=coefficients(summary(lmmod))[3,1]
  return(b)
}
glmfunctionp<-function(testdelta,testyear,testbio){
  lmmod<-lm(testdelta ~ testyear + testbio)
  p=coefficients(summary(lmmod))[3,4]
  b=coefficients(summary(lmmod))[3,1]
  return(p)
}
resglm<-
  dfclim %>% 
  group_by(snp) %>% 
  summarise(b_bio1= glmfunctionb(deltafreq, max_year,bio1),
            p_bio1= glmfunctionp(deltafreq, max_year,bio1)
            )


# QQ plot
ggplot(resglm)+
  geom_point(aes(y=sort(-log10(p_bio1)),
                 x=sort(-log10(runif(nrow(resglm))))
                 ))+
  geom_abline(slope = 1,intercept = 0)+
  xlab("Expected p-values")+
  ylab("Observed p-values")

# Allele to visualize
resglm %>% 
  arrange(p_bio1) %>%
  # arrange(snp) %>% 
  head
# 1_10001716
# myallele="1_10001716"
myallele="1_786902"
myallele="3_6191720"
myallele="1_17199262"

dfsub<-
  dfclim %>% 
  filter(snp==myallele)

alleleplot<-
ggplot(dfsub)+
  geom_hline(yintercept = 0,lty="dotted")+
  geom_point(aes(y=deltafreq,x=bio1, color=site))+
  stat_smooth(aes(y=deltafreq,x=bio1),
              method = "glm",
              color="grey")+
  scale_color_discrete(guide="none")+
  ylab(paste0("Frequency change (",myallele,")"))+
  xlab("Annual temperature (C)")+
  theme_minimal()
alleleplot

save_plot(alleleplot,
          filename=paste0(myfolder,"/fig-allele-climate-correlation-bio1-",myallele,".pdf"),
          base_height = 4.5,base_width = 5)

# parse chromosome and site
resglm$chr <- 
   unlist(lapply(resglm$snp, function(i) strsplit(i,split="_")[[1]][1])) %>% as.numeric
resglm$bp <- 
  unlist(lapply(resglm$snp, function(i) strsplit(i,split="_")[[1]][2])) %>% as.numeric

# manhattan
manhattanplot<-
ggmanhattan(data = resglm ,chr = "chr",bp = "bp",P = "p_bio1", SNP = "snp",
            sparsify=FALSE, 
            build="NA", # need to remove the default hg19 otherwise map incorrect
            expand.x = c(0.03, 0.1), expand.y = c(0.01,0.03),
            significance=0.05/nrow(resglm)
            )

manhattanplot

save_plot(manhattanplot,
          filename=paste0(myfolder,"/fig-allele-glm-manhattan-bio1.pdf"),
          base_height = 4.5,base_width = 12)
