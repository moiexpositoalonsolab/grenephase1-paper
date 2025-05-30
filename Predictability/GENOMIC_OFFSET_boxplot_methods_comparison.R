rm(list = ls())
setwd("~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/data-intermediate/Leave_one_out_prediction/")
library(dplyr)

replicates <- read.delim("replicates/Genomic_offset_plot_replication_prediction_summary.txt")
prs <- read.delim("stabilizing_selection_prs/Genomic_offset_stabilizing_selection_loo_p1_p0_prs_prediction_summary.txt")
#lfmm <- read.delim("~/Dropbox/Carnegie_Stanford/Projects/GrENE-net/genomic_offset/lfmm/Genomic_offset_site_lfmm.txt")
lfmm <- read.delim("binomial_regression/genomic_offset_firstgen_binom_reg_p1_p0_results.csv",sep=",") %>%
  mutate(source = "lfmm")
colnames(lfmm)[3] <- "sp_r"
colnames(lfmm)[5] <- "r2"

climate_distance <- read.delim("climate_distance/Genomic_offset_climate_distance_loo_prediction_summary.txt")



data <- rbind(replicates[,c(4,5,6)],climate_distance[,c(4,5,6)],lfmm[,c(3,5,16)],prs[,c(4,5,6)])
data$source <- factor(data$source,levels = c("lfmm","climate_distance","stabilizing_prs","plot_mean"))



ggplot(data,aes(x = source,y=sp_r))+
  geom_boxplot(fill="grey90")+
  theme_classic(base_size = 18)+
  theme(legend.position = "n")

ggplot(data,aes(x = source,y=r2))+
  geom_boxplot(fill="grey90")+
  theme_classic(base_size = 18)+
  theme(legend.position = "n")






era5 <- read.delim("stabilizing_selection_era5/Genomic_offset_stabilizing_selection_loo_era5_prediction_summary.txt")
data <- rbind(era5,replicates)
data$site <- as.factor(data$site)

replicate_summary <- replicates %>%
  group_by(site) %>%
  summarise(mean = mean(pearson_r))

selection_sites <- replicate_summary$site[replicate_summary$mean > 0.2]

replicates_selection <- replicates %>%
  filter(site %in% selection_sites)
prs_selection <- prs %>%
  filter(site %in%selection_sites )
lfmm_selection <- lfmm %>%
  filter(site %in%selection_sites )

climate_selection <- climate_distance %>%
  filter(site %in%selection_sites )



data_selection <- rbind(replicates_selection[,c(4,5,6)],climate_selection[,c(4,5,6)],lfmm_selection[c(3,4,5)],prs_selection[,c(4,5,6)])
data_selection$source <- factor(data_selection$source,levels = c("lfmm","climate_distance","stabilizing_prs","plot_mean"))
p1 <- ggplot(data_selection,aes(x = source,y=sp_r))+
  geom_boxplot(fill="grey90")+
  ylab("Spearman correlation (rho)")+
  xlab("")+
  theme_classic(base_size = 18)+
  theme(legend.position = "n",axis.text.x = element_text(angle = 60,vjust =1, hjust = 1))

p2 <- ggplot(data_selection,aes(x = source,y=r2))+
  geom_boxplot(fill="grey90")+
  ylab("prediction R2")+
  xlab("")+
  theme_classic(base_size = 18)+
  theme(legend.position = "n",axis.text.x = element_text(angle = 60, vjust =1, hjust = 1))

p1 + p2








era5_summary <- era5 %>%
  group_by(site) %>%
  summarise(mean_r2 = mean(r2),
            mean_pearson = mean(pearson_r))

mean(era5_summary$mean_pearson)
data$site <- factor(data$site,levels = era5_summary$site[order(era5_summary$mean_r2,decreasing = F)])
ggplot(data[data$source=="stabilizing_era5",],aes(x = site,y = sp_r))+
  geom_boxplot()+
  xlab("Experimental garden")+
  scale_color_manual(values = c("grey50","salmon"))+
  ylab("spearman correlation between observed and \n predicted accession frequency change")+
  #annotate(geom = "text",x = 20,y=0.5,label="Upper bound",size=10,color="grey50")+
  theme_classic(base_size = 18)


mean(era5_summary$mean)
mean(era5_summary$mean_pearson[replicate_summary$mean>0.2])



