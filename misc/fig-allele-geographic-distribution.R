
#### If in google drive
# from google.colab import drive
# drive.mount('/content/drive')
# %load_ext rpy2.ipython

################################################################################
library(sp)
library(raster)
library(rbioclim)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(rgbif)
library(dplyr)
library(moiR)
library(tidyverse)


################################################################################
# Set map limits
xlim = c(-200, +200)
ylim = c(0, +90)
xlim = c(-200, +200)
ylim = c(0, +90)
xlim=c(-10.5,+ 53)
ylim=c(32,65)
Range=extent(c(xlim,ylim))

# Download present data w2
now<-rbioclim::getData(name="worldclim2",var='bio', res=2.5, path = "~/")
now<- now %>% crop(.,Range)

# Download future data Max Planck model
fut<-rbioclim::getData(name="CMIP5",var='bio', res=2.5,model='MP', year=50, rcp=85,path = "~/")
fut<- fut %>% crop(.,Range)
names(fut) <- names(now) # Fix names in future dataset
for(i in 1:11) fut[[i]]<-fut[[i]]/10 # Fix units of temp in future dataset (wordlclim2 temp is not Cx10)


################################################################################
# %%R

# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
myfolder<-"~/grenephase1-analyses"

# Specify the file path

library(moiR)

file_edelta = '/data/delta_ecotype_freq.txt'
file_e0 = '/data/merged_ecotype_frequency.txt'
file_climate_sites="/grene/data/worldclim_sitesdata.rda"
file_sites_info="/grene/data/sites.clim.rda"

file_path_pdelta = paste0(myfolder, '/data/merged_hapFIRE_allele_frequency_LDpruned.txt' )
file_path_p0 = paste0(myfolder,'/data/average_seedmix_p0_LDpruned.txt')

file_ecotypes_wodclim="/grene/data/worldclim_ecotypesdata.rda"
file_ecotypes_data="/grene/data/ecotypes_data.rda"


myallele="1_17199262" # example from allele-climate correlation
################################################################################

# Load frequencies 
pdelta<-read.table(file_path_pdelta, header=T)
p0<-read.table(file_path_p0, header=F)
p0$snp<-paste0(p0$V1,"_",p0$V2)
pdelta$snp<-p0$snp
p0<-p0 %>% rename(p0=V3)

# check dataset
head(pdelta)
dim(pdelta)
head(p0)
dim(p0)

################################################################################
# %%R

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



################################################################################
################################################################################
# Map SNP based on geographic locations of 1001 G
# Map SNP based on climate GWAs with GRENE_net

# Load climates
load(file_path_climate_sites)
load(file_path_sites_info)

worldclim_sitesdata<-worldclim_sitesdata %>% data.frame %>% mutate_all(as.numeric)
sites.clim<-sites.clim %>% data.frame %>%  mutate_all(as.numeric)


# Merge
dfclim<-
  merge(df, worldclim_sitesdata,by="site")

################################################################################
# GLM and GRENE_net 

allelepalete<-
c(
  rev(brewer.pal(9,"Reds")[2:9]),
  "grey95",
  brewer.pal(9,"Greens")[2:9]
)

# This was an allele found in the naive GLM
myallele="1_17199262"

dfsub<-
  dfclim %>% 
  filter(snp==myallele)

alleleplot<-
  ggplot(dfsub)+
  geom_hline(yintercept = 0,lty="dotted")+
  geom_point(aes(y=deltafreq,x=bio1, color=deltafreq))+
  stat_smooth(aes(y=deltafreq,x=bio1),
              method = "glm",
              color="grey")+
  scale_color_gradientn(colours = allelepalete)+
  # scale_color_discrete(guide="none")+
  ylab(paste0("Frequency change (",myallele,")"))+
  xlab("Annual temperature (C)")+
  theme_minimal()
alleleplot
save_plot(alleleplot,
          filename=paste0(myfolder,"/fig-allele-climate-correlation-bio1-",myallele,".pdf"),
          base_height = 4.5,base_width = 5)


lmod<-lm(data=dfsub, deltafreq ~ bio1)
summary(lmod)

nowclimate<-now[[1]]
names(nowclimate)<-"bio1" # variable needs to be the same name
sdmpred= predict(nowclimate,lmod)

pdf(paste0(myfolder,"/fig-allele-climate-map-bio1-",myallele,".pdf"),height = 5,width = 5)
plot(sdmpred, col=allelepalete)
dev.off()
# save_plot(alleleplot,
#           filename=,
#           base_height = 4.5,base_width = 5)

################################################################################
# Now Maxent of allele frequencies fromt the 1001 genomes to know where alleles if ound

# Based on individual genomes
library(mar)
library(data.table)
setwd(myfolder)

famfile = paste0(myfolder,'/data/greneNet_final_v1.0.fam')
bimfile = paste0(myfolder,'/data/greneNet_final_v1.0.bim')
bedfile = paste0(myfolder,'/data/greneNet_final_v1.0.bed')
bedfile = "/Users/moisesexpositoalonso/grenephase1-analyses/data/greneNet_final_v1.0.bed"

fam<-data.table::fread(famfile)
map<-data.table::fread(bimfile)
bed<-readbed(bedfile,
             N = nrow(fam),p = nrow(map),
             myrows = 1:nrow(fam),mycols = 1:nrow(map))
# dim(bed)

genomesub<-bed[,which(map$V2==myallele)]
# rm(bed)

# Load the ecotype data
load(paste0(myfolder,file_ecotypes_data))
coords<-ecotypes_data %>% dplyr::select(longitude ,latitude)

# Make the SDM
library(dismo)
library(rJava)
# mod<-maxent(x = now, p=coords)

# save(mod,file= paste0(myfolder,"/intermediate-data/maxentmod-", myallele,".rda"))
load(paste0(myfolder,"/intermediate-data/maxentmod-", myallele,".rda"))

allelenow=predict(mod, now) 

pdf(paste0(myfolder,"/figs/fig-allele-climate-map-bioclim-maxent231G",myallele,".pdf"),height = 5,width = 5)
# plot(sdmpred, col=allelepalete)
plot(allelenow, col=allelepalete)
dev.off()
