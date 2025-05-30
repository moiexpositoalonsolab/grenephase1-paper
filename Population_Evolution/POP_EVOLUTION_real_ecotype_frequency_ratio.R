
rm(list=ls())

library(ggplot2)
library(dplyr)
library(patchwork)

## POPULATION EVOLUTION

## This script will analyze the real experimental ecotype frequency ecotype across 3 generations
## to show the changes in population composition and to indicate the strength of selection and drift

prefix <- "POP_EVOLUTION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

#####  load meta information, SNP and ecotype frequencies ######

meta <- read.delim("data/merged_sample_table.csv",sep=",",stringsAsFactors = F)
terminal_index <- read.delim("data/site_plot_terminal_index.txt")

## ecotype frequencies
ecotype_p0 <- read.delim("data/founder_ecotype_frequency.txt",header=F)
ecotype_p <- read.delim("data/merged_ecotype_frequency.txt",header=T,check.names = F)
ecotype_p[ecotype_p<0.001] <- 0

##### calculate the per generation delta frequency: compared to the previous generation

## define some global parameters

site_plot_index <- c(terminal_index$index1[terminal_index$generation1==1],
                     terminal_index$index2[terminal_index$generation1==1 & terminal_index$generation2==2],
                     terminal_index$index3[terminal_index$generation1==1 & terminal_index$generation2==2 & terminal_index$generation3==3])
experimental_sites <- unique(meta$site[site_plot_index])

experimental_meta <- meta[site_plot_index,]

num_ecotype <- nrow(ecotype_p0)
num_pop <- length(site_plot_index)
num_experiment_site <- length(experimental_sites)
founder_ecotype <- ecotype_p0$V1
founder_ecotype_frequency <- ecotype_p0$V2


## from generation 0 to generation 1
log_p1_p0 <- c()
gen1_index <- terminal_index$index1[terminal_index$generation1==1]
gen1_ecotype_frequency <- ecotype_p[,gen1_index]

for(i in 1:ncol(gen1_ecotype_frequency)){
  log_p1_p0 <- c(log_p1_p0,log(as.numeric(gen1_ecotype_frequency[i,]) / founder_ecotype_frequency[i]))
}
hist(log_p1_p0,breaks = 100)
mean(log_p1_p0)
## from generation 1 to generation 2
log_p2_p1 <- c()
gen2_index <- terminal_index$index2[terminal_index$generation1==1 & terminal_index$generation2==2 ]
gen1_index <- terminal_index$index1[terminal_index$generation1==1 & terminal_index$generation2==2 ]
gen2_ecotype_frequency <- ecotype_p[,gen2_index]
gen1_ecotype_frequency <- ecotype_p[,gen1_index]

for(i in 1:ncol(gen2_ecotype_frequency)){
  log_p2_p1 <- c(log_p2_p1,log(as.numeric(gen2_ecotype_frequency[i,]) / as.numeric(gen1_ecotype_frequency[i,])))
  #log_p2_p1 <- c(log_p2_p1,shannon_index(gen2_ecotype_frequency[,i]))
  
}
hist(log_p2_p1,breaks = 100)

## from generation 2 to generation 3
log_p3_p2 <- c()
gen3_index <- terminal_index$index3[terminal_index$generation1==1 & terminal_index$generation2==2 & terminal_index$generation3==3]
gen2_index <- terminal_index$index2[terminal_index$generation1==1 & terminal_index$generation2==2 & terminal_index$generation3==3]
gen3_ecotype_frequency <- ecotype_p[,gen3_index]
gen2_ecotype_frequency <- ecotype_p[,gen2_index]

log_p3_p2 <- c()
for(i in 1:ncol(gen3_ecotype_frequency)){
  log_p3_p2 <- c(log_p3_p2,log(as.numeric(gen3_ecotype_frequency[i,]) / as.numeric(gen2_ecotype_frequency[i,])))
  #log_p3_p2 <- c(log_p3_p2,shannon_index(gen3_ecotype_frequency[,i]))
  
}
hist(log_p3_p2,breaks = 100)

log_p_ratio <- c(log_p1_p0,log_p2_p1,log_p3_p2)


ecotype_frequency_ratio <- as.data.frame(matrix(nrow=length(log_p_ratio),ncol = 0)) %>%
  dplyr::mutate(log_p_ratio = log_p_ratio) %>%
  dplyr:::mutate(generation = c(rep(1,length(log_p1_p0)),rep(2,length(log_p2_p1)),rep(3,length(log_p3_p2))))
mean_ecotype_frequency_ratio <- ecotype_frequency_ratio %>%
  group_by(generation) %>%
  dplyr::summarise(mean_value = mean(log_p_ratio[!(is.na(log_p_ratio) | is.infinite(log_p_ratio))]))


p1 <- ggplot(ecotype_frequency_ratio,aes(x=log_p_ratio))+
  geom_histogram(alpha = 0.9,position = 'identity',bins = 30,fill="grey40")+
  #geom_density(alpha=0.8,fill="grey80")+
  facet_wrap(~generation)+
  geom_vline(xintercept = 0,linetype="dashed", color = "black")+
  geom_text(data = mean_ecotype_frequency_ratio, aes(x = 6, y = 10000, label = paste("Mean =", round(mean_value, 3))), vjust = -1) +
  ylab("Frequency")+
  theme_minimal()


p2 <- ggplot(ecotype_frequency_ratio,aes(x=log_p_ratio))+
  #geom_histogram(alpha = 0.9,position = 'identity',bins = 30,fill="grey40")+
  geom_density(alpha=0.8,fill="grey80")+
  facet_wrap(~generation)+
  geom_vline(xintercept = 0,linetype="dashed", color = "black")+
  geom_text(data = mean_ecotype_frequency_ratio, aes(x = 6, y = 0.2, label = paste("Mean =", round(mean_value, 3))), vjust = -1) +
  ylab("Density")+
  theme_minimal()

p1
p2
combined <- p1/p2
combined
ggsave(plot = combined,filename =  paste0("figs/",prefix,"real_grene_net_ecotype_frequency_ratio_across_years.pdf"),device = "pdf",width = 10,height = 8)

ggsave(plot = p2,filename =  paste0("figs/",prefix,"real_grene_net_ecotype_frequency_ratio_across_years_density.pdf"),device = "pdf",width = 8,height = 2)


## Now we can plot the per site per generation
## list N (experimental sites)
## list G (generation)


## generate N (site)
N <- c(rep(experimental_meta$site[experimental_meta$generation==1],231),
       rep(experimental_meta$site[experimental_meta$generation==2],231),
       rep(experimental_meta$site[experimental_meta$generation==3],231))

length(N)

ecotype_frequency_ratio <- ecotype_frequency_ratio %>%
  mutate(site=N)


RColorBrewer::brewer.pal(n = 9,name = "Greens")
cols <- c("generation 1"="C7E9C0","generation 2"="#41AB5D","generation 3"="#00441B")

ggplot(ecotype_frequency_ratio,aes(x=log_p_ratio))+
  geom_histogram(data=ecotype_frequency_ratio[ecotype_frequency_ratio$generation==1,],aes(x=log_p_ratio),alpha=0.6,fill="#C7E9C0",colour = F)+
  geom_histogram(data=ecotype_frequency_ratio[ecotype_frequency_ratio$generation==2,],aes(x=log_p_ratio),alpha=0.6,fill="#41AB5D",colour = F)+
  geom_histogram(data=ecotype_frequency_ratio[ecotype_frequency_ratio$generation==3,],aes(x=log_p_ratio),alpha=0.6,fill="#00441B",colour = F)+
  facet_wrap(~site)+
  geom_vline(xintercept = 0,linetype="dashed", color = "grey50")+
  ylab("Density")+
  theme_minimal()

ggsave(filename =  paste0("figs/",prefix,"real_grene_net_ecotype_frequency_ratio_across_sites_histogram.pdf"),device = "pdf",width = 10,height = 8)


## Now we are going to simulate the ecotype frequency change under random drift only given the known pop size
ecotype_p0

## generation 1

sample_size <- experimental_meta$total_flower_counts[experimental_meta$generation==1]

random_drift_df_1 <- as.data.frame(matrix(0,nrow=length(founder_ecotype),ncol=length(sample_size)+1))
random_drift_df_1[,1] <- founder_ecotype

for(i in 1:length(sample_size)){
  random_samples_ <- sample(1:231,sample_size[i],replace = T,prob = founder_ecotype_frequency)
  for(j in random_samples_){
    random_drift_df_1[j,i+1] <- random_drift_df_1[j,i+1] + 1
  }
}

random_drift_df_1_frequency <- apply(random_drift_df_1[,c(-1)],2,function(x) x/sum(x))

log_p1_p0_drift <- c()

for(i in 1:ncol(random_drift_df_1_frequency)){
  #log_p1_p0_drift <- c(log_p1_p0_drift,shannon_index(random_drift_df_1_frequency[,i]))
  log_p1_p0_drift <- c(log_p1_p0_drift,log(as.numeric(random_drift_df_1_frequency[i,]) / founder_ecotype_frequency[i]))
}
hist(log_p1_p0_drift,breaks = 100)
mean(log_p1_p0_drift)
