rm(list=ls())

library(ggplot2)
library(dplyr)

## POPULATION EVOLUTION

## This script will investigate which ecotypes increased or decreased their frequencies and is that related to climate

prefix <- "POP_EVOLUTION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

## first we need to establish the null frequency change expectation under drift for each ecotype

ecotype_p0 <- read.delim("data/founder_ecotype_frequency.txt",header=F)

meta <- read.delim("data/merged_sample_table.csv",sep=",")
terminal_index <- read.delim("data/site_plot_terminal_index.txt")

ecotype_climate <- read.delim("data-intermediate/CLIMATE_GWAS_bslmm_prs_ecotype_climate.txt")
#ecotype_climate <- ecotype_climate %>%
#  filter(ecotype %in% ecotype_p0$V1) %>%
#  arrange(match(ecotype,ecotype_p0$V1))


generation_1_parallelism <- read.delim("data-intermediate/generation_1_parallelism.txt")


gen1_index <- terminal_index$index1[terminal_index$generation1==1]

N <- meta$total_flower_counts[gen1_index]

ecotype_freq_sim <- matrix(0,nrow=231,ncol=500)
for(i in 1:500){
  sampled_size <- sample(N,1)
  sampled_ecotype <- sort(sample(1:231,sampled_size,replace = T,prob = ecotype_p0$V2))
  for(j in sampled_ecotype){
    ecotype_freq_sim[j,i] = ecotype_freq_sim[j,i] + 1
  }
  ecotype_freq_sim[,i] <- ecotype_freq_sim[,i] / sampled_size
  print(i)
}
expected_ecotype_freq_drift <- apply(ecotype_freq_sim,1,mean)
plot(expected_ecotype_freq_drift,ecotype_p0$V2)
abline(a=0,b=1)
ecotype_p0$expected_ecotype_freq_drift <- expected_ecotype_freq_drift

## now we can compare ecotypes in each site to the null expectation

ecotype_p <- read.delim("data/merged_ecotype_frequency.txt",check.names = F)

gen1_meta <- meta[gen1_index,] %>%
  filter(site != 33)
gen1_ecotype_p <- ecotype_p[,gen1_index]
gen1_ecotype_p <- gen1_ecotype_p[,-182]
table(colnames(gen1_ecotype_p) == gen1_meta$sample_name)

unique_sites <- unique(gen1_meta$site)

length(unique_sites)

site_climate <- read.delim("data-external/bioclimvars_sites_era5_year_2018.csv",sep=",") %>%
  filter(site %in% unique_sites) %>%
  arrange(match(site,unique_sites))

# ecotype_t_matrix <- matrix(nrow=231,ncol=length(unique_sites))
#
#
# for(i in 1:length(unique_sites)){
#   site_ecotype_p <- gen1_ecotype_p[,gen1_meta$site==unique_sites[i]]
#   for(j in 1:231){
#     m <- t.test(x=site_ecotype_p[j,],y=ecotype_freq_sim[j,])
#     ecotype_t_matrix[j,i] <- m$statistic
#   }
# }
#
# hist(ecotype_t_matrix)
#
#
# table(ecotype_t_matrix>1.96)
#
#
# df <- as.data.frame(matrix(nrow=231*30,ncol=5))
# colnames(df) <- c("ecotype","site","ecotype_bio1","site_bio1","t")
# df$ecotype <- rep(1:231,n=30)
# df$site <- rep(1:30,each=231)
# df$ecotype_bio1 <- rep(ecotype_climate$PRS_bio1,n=30)
# df$site_bio1 <- rep(site_climate$bio1,each=231)
# df$t <- as.numeric(ecotype_t_matrix)

#
# ggplot(df) +
#   geom_point(data=df[df$t > 4,],aes(x=site_bio1,y = ecotype_bio1),col="red")+
#   geom_abline(slope = 1)+
#   xlim(c(5,20))+
#   ylim(c(-3,20))+
#   theme_classic()
#

max_ecotype <- apply(gen1_ecotype_p,2,function(x) order(x,decreasing = T)[1:10])

max_ecotype_bio1_mean <- c()
for(i in 1:325){
  max_ecotype_bio1_mean <- c(max_ecotype_bio1_mean,mean(ecotype_climate$PRS_bio1[max_ecotype[,i]]))
}


gen1_meta_updated <- gen1_meta %>%
  left_join(site_climate[,1:2],by="site") %>%
  left_join(generation_1_parallelism[generation_1_parallelism$source=="ecotype",1:2],by="site")

plot(gen1_meta_updated$bio1,max_ecotype_bio1_mean,xlim=c(5,22),ylim=c(-3,20))
abline(a=0,b=1)
abline(lm(max_ecotype_bio1_mean~gen1_meta_updated$bio1),col="red")
gen1_meta_updated$mean[gen1_meta_updated$mean < 0] <- 0.01
summary(lm(max_ecotype_bio1_mean~gen1_meta_updated$bio1,weights = gen1_meta_updated$mean))
m <- lm(max_ecotype_bio1_mean~gen1_meta_updated$bio,weights = gen1_meta_updated$mean)
abline(m,col="blue")

df <- as.data.frame(matrix(nrow=325,ncol=4))
colnames(df) <- c("max_ecotype_bio1","site_bio1","parallelism","site")
df$max_ecotype_bio1<- max_ecotype_bio1_mean
df$site_bio1 <- gen1_meta_updated$bio1
df$parallelism <- gen1_meta_updated$mean
df$site <- gen1_meta_updated$site

summary(lm(data = df[df$parallelism>0.2,],formula =  max_ecotype_bio1~site_bio1 ,weights = parallelism))
summary(lm(data = df[df$parallelism<0.2,],formula = max_ecotype_bio1 ~site_bio1 ,weights = parallelism))

range(df$max_ecotype_bio1)
range(df$site_bio1)


df_summary <- df%>%
  group_by(site) %>%
  summarise(mean_ecotype = mean(max_ecotype_bio1),
            sd_ecotype = sd(max_ecotype_bio1),
            site = mean(site_bio1),
            r = mean(parallelism))


df %>%
  filter(parallelism>0.2) %>%
  ggplot(.)+
  geom_smooth(method = "lm",aes(x=site_bio1,y=max_ecotype_bio1,weight = parallelism))+
  geom_point(data=df_summary[df_summary$r>0.2,],aes(x=site,y=mean_ecotype),pch=16,cex=5)+
  geom_errorbar(data=df_summary[df_summary$r>0.2,],mapping=aes(x = site,ymin=mean_ecotype-sd_ecotype,ymax=mean_ecotype+sd_ecotype),width = 0.1)+
  annotate(geom = "text",label="y = 8.77 + 0.28x \n (R2: 0.240 P-value: 1.747e-9)",x=17,y=11,size=5)+
  xlab("Site Annual Mean Temp") +
  xlim(range(df_summary$site))+
  ylim(c(4,16))+
  ylab("Top Ecotype Annual Mean Temp")+
  theme_classic()

df %>%
  filter(parallelism<0.2) %>%
  ggplot(.)+
  geom_smooth(method = "lm",aes(x=site_bio1,y=max_ecotype_bio1,weight = parallelism))+
  geom_point(data=df_summary[df_summary$r<0.2,],aes(x=site,y=mean_ecotype),pch=16,cex=5,color="grey70")+
  geom_errorbar(data=df_summary[df_summary$r<0.2,],mapping=aes(x = site,ymin=mean_ecotype-sd_ecotype,ymax=mean_ecotype+sd_ecotype),width = 0.1,col="grey70")+
  annotate(geom = "text",label="y = 9.68 + 0.02x \n (R2: 0 P-value: 0.369)",x=17,y=11,size=5)+
  xlab("Site Annual Mean Temp") +
  xlim(range(df_summary$site))+
  ylim(c(4,16))+
  ylab("Top Ecotype Annual Mean Temp")+
  theme_classic()


df %>%
  ggplot(.)+
  ## sites under drift
  geom_smooth(data=df[df$parallelism<0.2,],method = "lm",aes(x=site_bio1,y=max_ecotype_bio1,weight = parallelism),col="grey80",se = F)+
  geom_point(data=df_summary[df_summary$r<0.2,],aes(x=site,y=mean_ecotype),pch=16,cex=5,color="grey80",)+
  geom_errorbar(data=df_summary[df_summary$r<0.2,],mapping=aes(x = site,ymin=mean_ecotype-sd_ecotype,ymax=mean_ecotype+sd_ecotype),width = 0.1,col="grey80")+
  annotate(geom = "text",label="y = 9.68 + 0.02x \n (R2: 0 P-value: 0.369)",x=17,y=7,size=6,color="grey80")+
  ## sites under selection
  geom_smooth(data=df[df$parallelism>0.2,],method = "lm",aes(x=site_bio1,y=max_ecotype_bio1,weight = parallelism),col="red",se=F)+
  geom_point(data=df_summary[df_summary$r>0.2,],aes(x=site,y=mean_ecotype),pch=16,cex=5,color="red")+
  geom_errorbar(data=df_summary[df_summary$r>0.2,],mapping=aes(x = site,ymin=mean_ecotype-sd_ecotype,ymax=mean_ecotype+sd_ecotype),width = 0.1,col="red")+
  annotate(geom = "text",label="y = 8.77 + 0.28x \n (R2: 0.240 P-value: 1.747e-9)",x=10,y=15,size=6,color="red")+
  geom_abline(slope = 1,lwd=1,lty=2)+
  xlab("Annual Mean Temp of Experimental Sites") +
  xlim(range(df_summary$site))+
  ylim(c(5,16))+
  ylab("Annual Mean Temp of \n Best-Performing Accessions")+
  theme_classic(base_size = 18)

ggsave(paste0("figs/",prefix,"site_best_performing_ecotype_bio1.pdf"),device = "pdf",width = 8,height = 6,units = "in")
ggsave(paste0("figs/",prefix,"site_best_performing_ecotype_bio1.png"),device = "png",width = 8,height = 6,units = "in",dpi = 600)

ggsave("~/Dropbox/Carnegie_Stanford/Conferences/BAPG2024/best_performing_ecotype_bio1.pdf",device = "pdf",width = 10,height = 7,units = "in")






## look at the pc space

ecotype_climate_pc <- prcomp(ecotype_climate[,2:20],center = T,scale = T)
site_climate_pc <- scale(site_climate[,2:20],center=ecotype_climate_pc$center,scale=ecotype_climate_pc$scale) %*% ecotype_climate_pc$rotation
site_climate$pc1 <- site_climate_pc[,1]

max_ecotype_pc_mean <- c()
for(i in 1:325){
  max_ecotype_pc_mean <- c(max_ecotype_pc_mean,mean(ecotype_climate_pc$x[max_ecotype[,i],1]))
}

gen1_meta_updated <- gen1_meta %>%
  left_join(site_climate[,c(1,21)],by="site") %>%
  left_join(generation_1_parallelism[generation_1_parallelism$source=="snp",1:2],by="site")

plot(gen1_meta_updated$pc1,max_ecotype_pc_mean,xlim=range(gen1_meta_updated$pc1),ylim=range(ecotype_climate_pc$x[,1]))
abline(a=0,b=1)
abline(lm(max_ecotype_pc_mean~gen1_meta_updated$pc1),col="red")
summary(lm(max_ecotype_pc_mean~gen1_meta_updated$pc1))


gen1_meta_updated$mean[ gen1_meta_updated$mean < 0] <- 0.001
summary(lm(max_ecotype_pc_mean~gen1_meta_updated$pc1,weights = gen1_meta_updated$mean))
m <- lm(max_ecotype_pc_mean~gen1_meta_updated$pc,weights = gen1_meta_updated$mean)
abline(m,col="blue")







plot(site_climate$bio1,apply(ecotype_t_matrix,2,function(x) sum(x<0.05)))
## load ecotype climate and site climate

ecotype_climate <- read.delim("data-intermediate/CLIMATE_GWAS_bslmm_prs_ecotype_climate.txt")
site_climate <- read.delim("data-external/bioclimvars_sites_era5_year_2018.csv",sep=",")
site_climate <- site_climate %>%
  filter(site %in% unique_sites) %>%
  arrange(match(site,unique_sites))

mean(site_climate$bio1)
mean(ecotype_climate$PRS_bio1)


ecotype_site_climate_distance <- matrix(nrow=231,ncol=30)

for(i in 1:length(unique_sites)){
  for(j in 1:231){
    ecotype_site_climate_distance[j,i] <- (site_climate$bio1[i] - ecotype_climate$PRS_bio1[j])
  }
}

plot(as.numeric(ecotype_site_climate_distance),as.numeric(ecotype_t_matrix))


plot(ecotype_climate$PRS_bio1,rep(y,231))
