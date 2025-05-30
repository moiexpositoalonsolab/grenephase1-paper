# library(grene)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(tidyr)
theme_set(theme_cowplot())

myfolder<-"~/grenephase1-analyses/"
# setwd("~/Google Drive/RESEARCH/MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/figures-raw/plot_ecotype_trend/")
setwd(myfolder)


################################################################################
# Load alleles already parsed
load("data-intermediate/allele_allyears_frequencies_long_raw.rda")
a<-allele_allyears_frequencies_long_raw

# Load worldclim
load("grene/data/worldclim_sitesdata.rda")

################################################################################
# Get one site example
site4start<-
    a %>% 
    dplyr::filter(site==4) %>%
    dplyr::filter(year==1) %>% # select one year example
    mutate(freq=startfreq) %>%
    mutate(year=0) %>%
    dplyr::select(snp,site,year,rep,freq, startfreq)
head(site4start)  
dim(site4start)  

# Subset to example site
site4all <-
  a %>% 
  dplyr::filter(site==4) %>% 
  rbind(.,site4start)
site4rep12<-
  site4all %>% 
  dplyr::filter(rep==12)

head(site4all)
tail(site4all)

################################################################################
# Get consistency per SNP with GLM
getbeta<-function(freq,yr)
{
    y=asin(as.numeric(freq))
    x=as.numeric(yr)
    effectsize=coefficients(summary(lm(y~x)))[2,1] / coefficients(summary(lm(y~x)))[2,2]
}

site4consistency<-
    site4all %>%
    group_by(snp) %>%
    summarise(consistency=getbeta(freq,year)) %>%
    arrange(consistency)

hist(site4consistency$consistency)

################################################################################
# Make a plot

site4all %>%
  dplyr::filter(snp %in% sample(unique(site4all$snp),size = 100) ) %>% # random sampling
  dplyr::filter(rep==1 ) %>% # random sampling
  ggplot()+
  geom_line(aes(y=freq,x=year,group=snp))


# plot one site for example, 100 random, 100 consisten SNPs
tmp1<-
  site4rep12 %>%
    dplyr::filter(snp %in% head(site4consistency$snp, 500)) # low sampling
tmp2<-
  site4rep12 %>%
    dplyr::filter(snp %in% tail(site4consistency$snp, 500))  # high sampling
tmp3<-
  site4rep12 %>%
    dplyr::filter(snp %in% sample(site4consistency$snp, 500))  # random sampling

# not scaled
plot4<-
    ggplot()+
    # geom_line(data=tmp3,aes(y=freq,x=year,group=snp), col=moiR::transparent("black",alpha = 0.1))+
    geom_line(data=tmp1,aes(y=freq,x=year,group=snp), col=moiR::transparent("#CB2026",alpha = 0.1))+
    geom_line(data=tmp2,aes(y=freq,x=year,group=snp), col=moiR::transparent("#6AA84E",alpha = 0.1))+
  theme_classic()
plot4

# save_plot("figs/fig-allele_frequencies_trendplots-examplesites-glm.pdf",plot4,
#           base_width = 4,base_height = 3 )

save_plot("figs/fig-allele_frequencies_trendplots-examplesites-glm-site4replicate12.pdf",plot4,
          base_width = 4,base_height = 3 )




################################################################################
################################################################################
## FUNCTION

allele_trendplot<-function(mysite=4, myreplicate=12, saveplot=T, choosetop=T){
  
################################################################################
message("Subset data of site")
# Get one site example
site4start<-
  a %>% 
  dplyr::filter(site==mysite) %>%
  dplyr::filter(year==1) %>% # select one year example
  mutate(freq=startfreq) %>%
  mutate(year=0) %>%
  dplyr::select(snp,site,year,rep,freq, startfreq)

# Subset to example site
site4all <-
  a %>% 
  dplyr::filter(site==mysite) %>% 
  rbind(.,site4start)

site4rep12<-
  site4all %>% 
  dplyr::filter(rep==myreplicate)

# MAKE A LOOP IN CASE REPLICATE IS NOT GOOD
if(choosetop==T){
  
  tmpreps<-
  site4all %>% 
    group_by(rep) %>% 
    summarise(totaldata=length(freq))

  toprep<-arrange(tmpreps,-totaldata,descending=T)$rep[1]
  message("Using replicate ", toprep," which contains most data")
  myreplicate=toprep
  site4rep12<-
    site4all %>% 
    dplyr::filter(rep==myreplicate)
}


################################################################################
# Get consistency per SNP with GLM
message("Calculate GLMs")
getbeta<-function(freq,yr)
{
  y=asin(as.numeric(freq))
  x=as.numeric(yr)
  effectsize=coefficients(summary(lm(y~x)))[2,1] / coefficients(summary(lm(y~x)))[2,2]
}

site4consistency<-
  site4all %>%
  group_by(snp) %>%
  summarise(consistency=getbeta(freq,year)) %>%
  arrange(consistency)

################################################################################
# Make a plot with random and consistent snps
message("Plot trends")
# plot one site for example, 100 random, 100 consisten SNPs
tmp1<-
  site4rep12 %>%
  dplyr::filter(snp %in% head(site4consistency$snp, 100)) # low sampling
tmp2<-
  site4rep12 %>%
  dplyr::filter(snp %in% tail(site4consistency$snp, 100))  # high sampling
tmp3<-
  site4rep12 %>%
  dplyr::filter(snp %in% sample(site4consistency$snp, 100))  # random sampling

# not scaled
plot4<-
  ggplot()+
  # geom_line(data=tmp3,aes(y=freq,x=year,group=snp), col=moiR::transparent("black",alpha = 0.1))+
  geom_line(data=tmp1,aes(y=freq,x=year,group=snp), col=moiR::transparent("#CB2026",alpha = 0.1))+
  geom_line(data=tmp2,aes(y=freq,x=year,group=snp), col=moiR::transparent("#6AA84E",alpha = 0.1))+
  labs(y="Allele frequency", x="Generation")+
  theme_classic()
plot4

tmp1$class="positive"
tmp2$class="negative"
tmp3$class="neutral"

returndata<-
  rbind(tmp1,tmp2, tmp3)

if(saveplot){
save_plot(
          paste0("figs/fig-allele_frequencies_trendplots-examplesites-glm-site",mysite,"replicate",myreplicate,".pdf"),
          plot4,
          base_width = 4,base_height = 3 )
}
return(returndata)
}

################################################################################
# Test with site4 replicate 12
mysites<-unique(a$site)
customordertemperature<-
  c(
    24,27,48,1,42,46,25,55,49,52,33,57,
    53, 54, 9, 12, 11, 6, 5, 2, 45,28,4, 32, 43, 13, 10
  )
length(mysites)

allele_trendplot(mysite=4,myreplicate=12)
allele_trendplot(mysite=25,myreplicate=1)

for(i in mysites){
  message(i)
  allele_trendplot(mysite=i,myreplicate=1)
}

tmpplot<-
  allele_trendplot(mysite=27,myreplicate=2)

# Selected site 4 examples

at1<-allele_trendplot(mysite=4,myreplicate=1,saveplot=F,choosetop = F)
at2<-allele_trendplot(mysite=4,myreplicate=2,saveplot=F,choosetop = F)
at3<-allele_trendplot(mysite=4,myreplicate=3,saveplot=F,choosetop = F)
at4<-allele_trendplot(mysite=4,myreplicate=4,saveplot=F,choosetop = F)
at5<-allele_trendplot(mysite=4,myreplicate=5,saveplot=F,choosetop = F)
at6<-allele_trendplot(mysite=4,myreplicate=6,saveplot=F,choosetop = F)
at7<-allele_trendplot(mysite=4,myreplicate=7,saveplot=F,choosetop = F)
at8<-allele_trendplot(mysite=4,myreplicate=8,saveplot=F,choosetop = F)
at9<-allele_trendplot(mysite=4,myreplicate=9,saveplot=F,choosetop = F)
at10<-allele_trendplot(mysite=4,myreplicate=10,saveplot=F,choosetop = F)
at11<-allele_trendplot(mysite=4,myreplicate=11,saveplot=F,choosetop = F)
at12<-allele_trendplot(mysite=4,myreplicate=12,saveplot=F,choosetop = F)


#  myexamples<- 
#    list(
#      allele_trendplot(mysite=4,myreplicate=1,saveplot=F,choosetop = F),
#      allele_trendplot(mysite=4,myreplicate=2,saveplot=F,choosetop = F),
#      allele_trendplot(mysite=4,myreplicate=3,saveplot=F,choosetop = F),
#      allele_trendplot(mysite=4,myreplicate=4,saveplot=F,choosetop = F),
#      allele_trendplot(mysite=4,myreplicate=5,saveplot=F,choosetop = F),
#      allele_trendplot(mysite=4,myreplicate=6,saveplot=F,choosetop = F),
#      allele_trendplot(mysite=4,myreplicate=7,saveplot=F,choosetop = F),
#      allele_trendplot(mysite=4,myreplicate=8,saveplot=F,choosetop = F),
#      allele_trendplot(mysite=4,myreplicate=9,saveplot=F,choosetop = F),
#      allele_trendplot(mysite=4,myreplicate=10,saveplot=F,choosetop = F),
#      allele_trendplot(mysite=4,myreplicate=11,saveplot=F,choosetop = F),
#      # allele_trendplot(mysite=4,myreplicate=12,saveplot=F,choosetop = F),
#      # allele_trendplot(mysite=4,myreplicate=12,saveplot=F,choosetop = F)
#    )
myexamples<-list(
at1,at2,at3,at4,
at5,at6,at7,at8,
at9,at10,at11,at12
)
myexamples<-
  do.call(rbind,myexamples)

examplesreplicates<-
ggplot(myexamples)+
  # geom_line(data=tmp1,aes(y=freq,x=year,group=snp), col=moiR::transparent("#CB2026",alpha = 0.1))+
  # geom_line(data=tmp2,aes(y=freq,x=year,group=snp), col=moiR::transparent("#6AA84E",alpha = 0.1))+
  geom_line(aes(y=freq,x=year,group=snp, color=class), 
            alpha = 0.1)+
  scale_color_manual(values = c("positive"="#CB2026", "negative"="#6AA84E", "neutral"="grey"))+
  theme_minimal()+
  theme(legend.position = "none",
        strip.text = element_text(angle = 90, hjust = 0))+
  facet_wrap(~as.numeric(rep), ncol=3, nrow=4)+
  xlab("Years")+ylab("Allele frequency")
examplesreplicates

ggsave(paste0("figs/fig-allele_frequencies_trendplots-site4-severalreplicates.png"),
       examplesreplicates,
       width = 10,height = 10 )
ggsave(paste0("figs/fig-allele_frequencies_trendplots-site4-severalreplicates.pdf"),
       examplesreplicates,
       width = 10,height = 10 )




################################################################################
# Get the data for plots to do a facet grid

plotlist<-list()
for(i in mysites){
  message("Working on site ", i)
  message(i)
  plotlist[[i]]<- allele_trendplot(mysite=i,myreplicate=12,saveplot=F,choosetop = T)
}

# Collapse
alldata<-
  do.call(rbind,plotlist)
alldata %>% head
alldata %>% tail

# Merge alldata with bio1 to order sites
alldata<-
  alldata %>% 
  merge(.,dplyr::select(worldclim_sitesdata,site,bio1), by="site")

# Make the giant plot!

# sort factor for 
library(forcats)
sites_simple_names<-read.csv("grene/data/sites_simple_names.csv",header = T)
simplenames<-paste0(sites_simple_names$city,", ", sites_simple_names$country)
names(simplenames)<-sites_simple_names$site
# Recode the sites
alldata$site <- recode(as.character(alldata$site), !!!simplenames)
# Recoder sites and ecotypes based on bio1
alldata$site <- fct_reorder(alldata$site, alldata$bio)
# alldata$id <- fct_reorder(alldata$id, alldata$ecotypes_bio1)

# plot
alltrends<-
ggplot(alldata)+
  # geom_line(data=tmp1,aes(y=freq,x=year,group=snp), col=moiR::transparent("#CB2026",alpha = 0.1))+
  # geom_line(data=tmp2,aes(y=freq,x=year,group=snp), col=moiR::transparent("#6AA84E",alpha = 0.1))+
  geom_line(aes(y=freq,x=year,group=snp, color=class), 
            alpha = 0.1)+
  scale_color_manual(values = c("positive"="#CB2026", "negative"="#6AA84E", "neutral"="grey"))+
  theme_minimal()+
  theme(legend.position = "none",
        strip.text = element_text(angle = 90, hjust = 0))+
  facet_grid(~site)+
  xlab("Years")+ylab("Allele frequency")


ggsave(paste0("figs/fig-allele_frequencies_trendplots-allsites-topreplicate.png"),
       alltrends,
       width = 20,height = 5 )
ggsave(paste0("figs/fig-allele_frequencies_trendplots-allsites-topreplicate.pdf"),
       alltrends,
       width = 20,height = 5 )


