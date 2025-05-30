################################################################################
### Goal
### Estimate some basic metrics of ecotypes increasing in frequency across
### sites, means, variances, optimal places. Including poulation size

################################################################################
library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(cowplot)

################################################################################
#### Location
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
# myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)


################################################################################
# Load ecotypes

load(file="data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize_flowers.rda")
ls()

e<-enew <- ecotype_allyears_frequencies_long_raw_climate_popsize_flowers
enew[1:5,1:10]
colnames(e) %>% tail


#****************************************************************************#
repeatability_evolution<-function(df, possiblesites=unique(df$site), reps=10){

# Make an empty data frame to store the results
results <- data.frame()

# Do 10 random draws
for(i in 1:reps) {
  ## To not allow overlaps uncomment
  # mysites<- sample(possiblesites,2, replace = F)
  # site1<-mysites[1]
  # site2<-mysites[2]
  ## To allow overlaps
  site1<-sample(possiblesites,1, replace = F)
  site2<-sample(possiblesites,1, replace = F)

  # Subset data for 2 sites
  df4 <- df %>% filter(site == site1)
  df5 <- df %>% filter(site == site2)

  # Randomly pick one replicate from site 4
  random_rep4 <- sample(unique(df4$rep), 1)
  # Randomly pick one replicate from site 5
  random_rep5 <- sample(unique(df5$rep), 1)

  # Filter each subset down to the chosen replicate
  sub4 <- df4 %>% filter(rep == random_rep4) %>% dplyr::select(id,freq,startfreq)
  sub5 <- df5 %>% filter(rep == random_rep5) %>% dplyr::select(id,freq,startfreq)

  # Merge the two data frames by a common key, e.g., "year"
  # (If you need to match on a different column or combination of columns, change accordingly)
  merged_data <- merge(sub4, sub5, by = "id", suffixes = c("_4", "_5")) %>%
                dplyr::mutate(
                              logp_4=log(freq_4/startfreq_4),
                              logp_5=log(freq_5/startfreq_5)
                              ) %>%
                dplyr::filter(logp_4 != -Inf, logp_5 != -Inf)
  cor.test(merged_data$logp_4,merged_data$logp_5)$estimate

  # Compute correlation of freq_4 vs freq_5
  corr_val <- cor.test(merged_data$logp_4,merged_data$logp_5)$estimate
  corr_p <- cor.test(merged_data$logp_4,merged_data$logp_5)$p.value

  # Store the result
  results <- rbind(
    results,
    data.frame(
      iteration  = i,
      siteA= site1,
      siteB= site2,
      rapA  = random_rep4,
      repB  = random_rep5,
      correlation = corr_val,
      pval=corr_p
    )
  )
}

results

}
#****************************************************************************#

# Example Spanish sites
warmsites<-c(4,5,32)
warm<-e %>%
  dplyr::filter(site %in% warmsites) %>%
  dplyr::filter(year==1)
# Run repeatability 100 times with warm sites
res<-repeatability_evolution(warm,reps = 1000) %>%
      dplyr::mutate(within= siteA==siteB) %>%
      dplyr::mutate(within= ifelse(within,"within","between")) %>%
      dplyr::filter(correlation!=1)
head(res)


### Example cold sites '
coldsites<-c(46,49,54)
cold<-e %>%
  dplyr::filter(site %in% coldsites) %>%
  dplyr::filter(year==1)
# Run repeatability 100 times with warm sites
rescold<-repeatability_evolution(cold,possiblesites = coldsites,reps = 1000) %>%
  dplyr::mutate(within= siteA==siteB) %>%
  dplyr::mutate(within= ifelse(within,"within","between")) %>%
  dplyr::filter(correlation!=1)
head(rescold)


# PLOT
# c( "#B2182B" , "#D6604D" , "#F4A582" , "#FDDBC7" , "#F7F7F7" , "#D1E5F0" , "#92C5DE" , "#4393C3" , "#2166AC")


library(ggnewscale)

ggplot()+
  # the warm sites
  geom_density(data=res,aes(x=correlation, fill=factor(within)), alpha=0.5)+
  scale_fill_manual("",values = c("#B2182B", "#FCBBA1"))+
  # start a new scale
  new_scale_fill() +
  # the cold
  geom_density(data=rescold,aes(x=correlation, fill=factor(within)), alpha=0.5)+
  scale_fill_manual("",values = c( "#08306B",   "#9ECAE1"))+
  scale_x_continuous(limits = c(-0.5,1))+
  theme_minimal()

################################################################################
# Randomization

# Run repeatoing 100 times, with randomized data
fullyrandomized<-e %>%
  dplyr::filter(year==1) %>%
  dplyr::mutate(freq=sample(freq), startfreq=sample(startfreq))
# log(fullyrandomized$freq/fullyrandomized$startfreq) %>% hist

resrandom<-repeatability_evolution(df = fullyrandomized,
                                   possiblesites=unique(fullyrandomized$site),
                                   reps = 1000)%>%
  dplyr::mutate(within= siteA==siteB) %>%
  dplyr::mutate(within= ifelse(within,"within","between")) %>%
  dplyr::filter(correlation!=1)
head(resrandom)
tail(resrandom)

# # Run randomization and repeatability
nestedrandomized<-e %>%
  dplyr::filter(year==1) %>%
  group_by(site,rep) %>%
  dplyr::mutate(freq=sample(freq))
resnested<-repeatability_evolution(nestedrandomized,reps = 1000)%>%
  dplyr::mutate(within= siteA==siteB) %>%
  dplyr::mutate(within= ifelse(within,"within","between")) %>%
  dplyr::filter(correlation!=1)

# PLOT
ggplot()+
  geom_histogram(data=resrandom,aes(x=correlation, fill=factor(within)), alpha=0.5)+
  scale_fill_manual("",values = c("black", "grey"))+
  theme_minimal()

ggplot()+
  geom_histogram(data=resnested,aes(x=correlation, fill=factor(within)), alpha=0.5)+
  scale_fill_manual("",values = c("black", "grey"))+
  theme_minimal()


################################################################################
# PLot everything

nbins=50
myalpha=0.5

evo<-
ggplot()+
  # Neutral
  geom_density(data=resrandom,aes(x=correlation, fill=factor(within)), alpha=myalpha,bins = nbins, fill="black")+
  geom_density(data=resnested,aes(x=correlation, fill=factor(within)), alpha=myalpha,bins = nbins, fill="grey")+
  # scale_fill_manual("",values = c("black"))+
  new_scale_fill() +
  # start a new scale and plot cold
  geom_density(data=rescold,aes(x=correlation, fill=factor(within)), alpha=myalpha,bins = nbins)+
  scale_fill_manual("",values = c( "#08306B",   "#9ECAE1"))+
  new_scale_fill() +
  # the warm sites
  geom_density(data=res,aes(x=correlation, fill=factor(within)), alpha=myalpha, bins = nbins)+
  scale_fill_manual("",values = c("#B2182B", "#FCBBA1"))+
  # settings
  geom_vline(xintercept = 0,color="grey")+
  scale_x_continuous(limits = c(-0.5,1),breaks =seq(-0.4,1,0.2) )+
  xlab("Pearson's correlation (r)")+
  theme_minimal()
evo

ggsave("figs/fig-evo-parallelism-and-repeatability.pdf",
       evo,
       height = 5, width = 6, units="in")
ggsave("figs/fig-evo-parallelism-and-repeatability.png",
       evo,
       height = 5, width = 6, units="in")


sink("tables/parallelism-vs-repeatability.txt")
# print("warm neutral")

print("warm")
res %>%
  # head %>%
  group_by(within) %>%
  summarise(mean= t.test(correlation)$estimate,
            CI95up=t.test(correlation)$conf.int[2],
            CI95low=t.test(correlation)$conf.int[1]
            )

# print("cold neutral")

print("cold")
rescold %>%
  # head %>%
  group_by(within) %>%
  summarise(mean= t.test(correlation)$estimate,
            CI95up=t.test(correlation)$conf.int[2],
            CI95low=t.test(correlation)$conf.int[1]
  )

sink()
