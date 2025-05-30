rm(list=ls())

library(ggplot2)
library(dplyr)

## POPULATION EVOLUTION

## This script will plot the relationship between parallelism observed in the first generation and climate distantce
## and then plot the population size change over three generations

prefix <- "POP_EVOLUTION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)


generation_1_parallelism <- read.delim("data-intermediate/generation_1_parallelism.txt")
unique_sites <- unique(generation_1_parallelism$site)

## load the site climate data
site_climate <- read.delim("data-external/bioclimvars_sites_era5_year_2018.csv",sep=",")
site_climate <- site_climate %>%
  filter(site %in% unique_sites) %>%
  arrange(match(site,unique_sites))


## load the ecotype climate data
ecotype_climate <- read.delim("data-external/bioclimvars_ecotypes_era5.csv",sep = ",")
ecotype_p0 <- read.delim("data/founder_ecotype_frequency.txt",header=F)
ecotype_climate <- ecotype_climate %>%
  filter(ecotype %in% ecotype_p0$V1)

## load the population size 

pop_size <- read.delim("data-intermediate/pop_size_estimations.csv",sep=",")
pop_size_mean <- pop_size %>%
  filter(generation ==1 ) %>%
  filter(site %in% unique_sites) %>% 
  group_by(site) %>%
  dplyr::summarize(avg_flower = mean(flowerscollected_corrected))

pop_size_mean$avg_flower
site_climate$site[site_climate$site %in% pop_size_mean$site]


site_name <- read.delim("data/site_info.txt",check.names = F)

# ## run pca on ecotype and site climate
# ecotype_climate_pca <- prcomp(site_climate[,2:20],center = T,scale = T)
# summary(ecotype_climate_pca)
# 
# site_climate_scale <- scale(site_climate[,2:20],center=ecotype_climate_pca$center,scale = ecotype_climate_pca$scale)
# site_climate_pca <- site_climate_scale %*% ecotype_climate_pca$rotation
# 
# population_center <- apply(ecotype_climate_pca$x,2,mean)
# 
# 
# # climate distance between site and ecotype population center (PC1)
# distance <- (site_climate_pca[,1] - population_center[1])
# 
# y <- generation_1_parallelism$mean[generation_1_parallelism$source=="snp"] 
# summary(lm(y~poly(site_climate_pca[,1] ,2) ,weights = as.numeric(pop_size_mean$avg_flower)))
# 
# plot(distance,y,cex= as.numeric(pop_size_mean$avg_flower) / max(as.numeric(pop_size_mean$avg_flower)) * 5)
# 
# 
# 
# 
# plot(distance,y,cex=as.numeric(pop_size_mean$avg_flower)/max(as.numeric(pop_size_mean$avg_flower))*5)

# ## linear regression on bio14
# 
# population_center <- apply(ecotype_climate[,2:20],2,median)
# 
# # climate distance between site and ecotype population center (PC1)
# distance <- (site_climate$bio14 - population_center[14])^2
# y <- generation_1_parallelism$mean[generation_1_parallelism$source=="ecotype"] 
# summary(lm(y~distance  ,weights = as.numeric(pop_size_mean$avg_flower)))
# 
# 
# plot(distance,y,cex=as.numeric(pop_size_mean$avg_flower)/max(as.numeric(pop_size_mean$avg_flower))*5)
# 
# df <- as.data.frame(matrix(ncol=4,nrow=30))
# colnames(df) <- c("site","r","distance","pop_size")
# 
# df$site <- unique_sites
# df$r <- generation_1_parallelism$mean[generation_1_parallelism$source=="ecotype"] 
# df$distance <- distance
# df$pop_size <- as.numeric(pop_size_mean$avg_flower)
# 
# p1 <- ggplot(df,aes(x=distance,y=r))+
#   #geom_point(pch=16,cex = df$pop_size / median(df$pop_size ) * 2.5)+
#   geom_point(pch=16,cex = df$pop_size / 25,color="grey20")+
#   geom_smooth(method = "lm",mapping=aes(weight = pop_size))+
#   annotate("text",label="y=0.145 + 1.91x \n R2: 0.374 p-value 0.0001992",x = 1000,y=0.55,cex=5)+
#   xlab("Distance in Precipitation of Driest Month")+
#   ylab("Ecotype Parallelism")+
#   theme_classic()
# 
# legend_df <- df[df$site %in% c(48,2,53,23),]
# legend <- ggplot(legend_df)+
#   geom_point(x=1,y=c(1,2,3,4),cex=c(5,50,100,200)/25,color="grey20") +
#   xlim(c(0.9,1.1))+
#   ylim(c(0,5))+
#   geom_text(label = c(5,50,100,200),x = 1.02,y=c(1,2,3,4),cex=5)+
#   theme_classic()
# combined <- p1 | legend
# combined
# ggsave(plot = combined,filename =  paste0("figs/",prefix,"ecotype_parallelism_bio14_distance.pdf"),device = "pdf",width = 10,height = 5,units = "in")
# 


## now plot the pop size change over generations

## first load the meta data (only look at site_plots that are analyzed in the parallelism)

meta <- read.delim("data/merged_sample_table.csv",sep=",",stringsAsFactors = F)
terminal_index <- read.delim("data/site_plot_terminal_index.txt")

site_plot_index <- c(terminal_index$index1[terminal_index$generation1==1],
                     terminal_index$index2[terminal_index$generation1==1 & terminal_index$generation2==2],
                     terminal_index$index3[terminal_index$generation1==1 & terminal_index$generation2==2 & terminal_index$generation3==3])
experimental_meta <- meta[site_plot_index,]



pop_size_trimmed <- pop_size %>%
  filter(site %in% unique_sites) %>%
  filter(generation %in% c(1,2,3)) %>%
  mutate(sample_name = paste(site,generation,plot,sep="_")) %>%
  filter(sample_name %in% experimental_meta$sample_name)

#pop_size_trimmed$generation <- as.factor(pop_size_trimmed$generation)
ggplot(pop_size_trimmed,mapping=aes(x=factor(generation),y=flowerscollected_corrected))+
  geom_boxplot(width=0.3)+
  stat_summary(fun.y=mean, geom="point", shape=15, size=4,col="red")+
  xlab("Generation")+
  ylab("Estimated Population Size \n (Number of Flowers collected)")+
  geom_smooth(data = pop_size_trimmed,aes(x=generation,y=flowerscollected_corrected),formula = y~poly(x,2),method = "lm")+
  theme_minimal(base_size = 18)
ggsave(filename =  paste0("figs/",prefix,"estimated_pop_size_arcoss_generations.pdf"),device = "pdf",width = 4,height = 5,units = "in")

ggplot(pop_size_trimmed,mapping=aes(x=factor(generation),y=totalplantnumber_complete))+
  geom_boxplot(width=0.3)+
  stat_summary(fun.y=mean, geom="point", shape=15, size=2,col="red")+
  xlab("Generation")+
  ylab("Estimated Population Size")+
  geom_smooth(data = pop_size_trimmed,aes(x=generation,y=totalplantnumber_complete),formula = y~poly(x,2),method = "lm")+
  theme_minimal(base_size = 18,base_line_size = 0.2)
ggsave(filename =  paste0("figs/",prefix,"estimated_pop_size_arcoss_generations_totalplantnumbercomplete.pdf"),device = "pdf",width = 5,height = 4,units = "in")


experimental_meta$generation <- as.numeric(experimental_meta$generation)


## now plot the pop size change per site across all three generations
pop_size_summary <- as.data.frame(matrix(ncol = 5,nrow=30*3))
colnames(pop_size_summary) <- c("site","generation","mean","upper","lower")
pop_size_summary$generation <- as.factor(c(rep(1,30),rep(2,30),rep(3,30)))
pop_size_summary$site <- rep(unique(pop_size_trimmed$site),3)
pop_size_summary <- pop_size_summary %>%
  left_join(.,site_name[,c(1,2)],by="site")
order_site <- generation_1_parallelism$site[order(generation_1_parallelism$mean[generation_1_parallelism$source=="snp"],decreasing = T)]
pop_size_summary$site <- factor(pop_size_summary$site,levels = order_site)


for(i in 1:30){
  for(j in 1:3){
    index1 = which(pop_size_trimmed$site==unique(pop_size_trimmed$site)[i] & pop_size_trimmed$generation==j)
    if (length(index1) > 0){
      index2 = which(pop_size_summary$site == unique(pop_size_trimmed$site)[i] & pop_size_summary$generation==j)
      pop_size_summary[index2,3] <- mean(pop_size_trimmed$flowerscollected_corrected[index1])
      pop_size_summary[index2,4] <- mean(pop_size_trimmed$flowerscollected_corrected[index1]) + sd(pop_size_trimmed$flowerscollected_corrected[index1])
      pop_size_summary[index2,5] <- mean(pop_size_trimmed$flowerscollected_corrected[index1]) - sd(pop_size_trimmed$flowerscollected_corrected[index1])
      
    }
    
  }
}



ggplot(pop_size_summary,aes(x=site,y=mean,colour = generation))+
  coord_flip()+
  geom_point(position = position_dodge(.5),cex=3,pch=16)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1,position=position_dodge(.5),linewidth = 0.7)+
  scale_color_manual(values =c("grey20","grey50","grey80"))+
  ylab("Estimated Population Size \n (Number of Flowers collected)")+
  scale_x_discrete(breaks = pop_size_summary$site,labels = paste0("#",pop_size_summary$site," ",pop_size_summary$name))+
  xlab("GrENE-net garden")+
  theme_classic()
ggsave(filename =  paste0("figs/",prefix,"estimated_pop_size_per_site.pdf"),device = "pdf",width = 8,height = 10,units = "in")


site_32 <- pop_size_trimmed %>%
  filter(site==32)
site_32$generation <- as.numeric(site_32$generation)
summary(lm(site_32$flowerscollected_corrected~site_32$generation + I(site_32$generation^2)))

pop_size_trimmed$generation <- as.numeric(pop_size_trimmed$generation)
summary(lm(pop_size_trimmed$flowerscollected_corrected~pop_size_trimmed$generation + I(pop_size_trimmed$generation^2)))


ggplot(pop_size_trimmed,aes(x=generation,y = flowerscollected_corrected))+
  geom_point()+
  geom_smooth(formula = y~poly(x,2),method = "glm")+theme_classic()
