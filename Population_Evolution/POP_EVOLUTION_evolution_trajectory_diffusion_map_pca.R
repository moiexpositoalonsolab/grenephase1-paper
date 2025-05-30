rm(list=ls())

library(ggplot2)
library(dplyr)
library(patchwork)
library(GGally)
library(network)
library(sna)
library(igraph)
library(intergraph)
library(qgraph)
library(ggsci)
library(slingshot)
library(destiny)
library(mclust)
library(rgl)
library(scatterplot3d)
library(colorjam)
library(RColorBrewer)
library(plotly)
library(ggforce)

## POPULATION EVOLUTION

## This script will use the delta snp frequency to generate the diffusion map and PCA plots to
## describe the population change dynamics

prefix <- "POP_EVOLUTION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

#####  load meta information, SNP and ecotype frequencies ######

meta <- read.delim("data/merged_sample_table.csv",sep=",",stringsAsFactors = F)
terminal_index <- read.delim("data/site_plot_terminal_index.txt")

## SNP frequencies
ld_p0 <- read.delim("data/average_seedmix_p0_LDpruned.txt",header=F)
ld_p <- read.delim("data/merged_hapFIRE_allele_frequency_LDpruned.txt",header=T,check.names = F)

delta_p <- apply(ld_p,2,function(x) (x - ld_p0$V3))


## build PCA on the third generation and then project gen1 and gen2 onto the PC space

# delta_p_gen3 <- delta_p[,meta$generation==3]
# pca_gen3 <- prcomp(t(delta_p_gen3))
# summary(pca_gen3)
# pca_gen2 <- predict(pca_gen3,newdata = t(delta_p[,meta$generation==2]))
# pca_gen1 <- predict(pca_gen3,newdata = t(delta_p[,meta$generation==1]))
# pca_founder <- predict(pca_gen3,newdata= t(ld_p0$V3))
# 
# plot(pca_gen3$x[,1],pca_gen3$x[,2],xlim= c(-11.37498,28.95261),ylim=c(-12.68625,42.95254))
# plot(pca_gen2[,1],pca_gen2[,2],xlim= c(-11.37498,28.95261),ylim=c(-12.68625,42.95254))
# plot(pca_gen1[,1],pca_gen1[,2],xlim= c(-11.37498,28.95261),ylim=c(-12.68625,42.95254))

pca <- prcomp(t(delta_p))
pca_gen3 <- pca$x[meta$generation==3,]
pca_gen2 <- pca$x[meta$generation==2,]
pca_gen1 <- pca$x[meta$generation==1,]
pca_founder <- predict(pca,newdata= t(ld_p0$V3))


site_climate <- read.delim("data-external/bioclimvars_experimental_sites_era5.csv",sep=",")
site_climate <- site_climate %>%
  filter(site %in% unique(meta$site)) %>%
  arrange(match(site,unique(meta$site)))

## define the colors 
redblue<- rev(c( "#B2182B", "#B2182B" , "#D6604D" ,"#D6604D", "#F4A582" , "#F7F7F7" , "#D1E5F0" , "#92C5DE" , "#4393C3" , "#2166AC"))

quantile_custom <- function(x){
  q <- c()
  for(i in 1:length(x)){
    q <- c(q,round((x[i] - range(x)[1]) / (range(x)[2] - range(x)[1]) * 10))
  }
  q[q==0] <- 1
  return(q)
}


df_climate_quantile <- as.data.frame(matrix(nrow=nrow(site_climate),ncol=0)) %>%
  mutate(site = site_climate$site) %>%
  mutate(var1_quantile = quantile_custom(site_climate$bio1))


color <- c()
for(i in 1:nrow(df_climate_quantile)){
  color <- c(color,redblue[df_climate_quantile$var1_quantile[i]])
}
df_climate_quantile$plotcolor <- color


## make the dataframe

df3 <- tibble(
  site = meta$site[meta$generation==3],
  pc1 = pca_gen3[,1],
  pc2 = pca_gen3[,2]
) %>%
  left_join(.,site_climate,by="site") %>%
  left_join(.,df_climate_quantile,by="site")

df2 <- tibble(
  site = meta$site[meta$generation==2],
  pc1 = pca_gen2[,1],
  pc2 = pca_gen2[,2]
) %>%
  left_join(.,site_climate,by="site") %>%
  left_join(.,df_climate_quantile,by="site")


df1 <- tibble(
  site = meta$site[meta$generation==1],
  pc1 = pca_gen1[,1],
  pc2 = pca_gen1[,2]
) %>%
  left_join(.,site_climate,by="site") %>%
  left_join(.,df_climate_quantile,by="site")

cor.test(df3$pc1,df3$bio1)
cor.test(df3$pc2,df3$bio13)

range(pca$x[,1])
range(pca$x[,2])
summary(pca)

pca_gen3_plot <- df3 %>%
  ggplot(.) +
  geom_segment(aes(xend=pc1, yend=pc2,x=pca_founder[,1], y=pca_founder[,2],color=plotcolor))+
  scale_color_identity()+
  geom_point(aes(x=pc1, y=pc2), color="white", size=3.5)+
  geom_point(aes(x=pc1, y=pc2), color="black", size=3)+
  geom_point(aes(x=pc1, y=pc2, color=plotcolor),size=2.5)+
  geom_point(aes(x=pca_founder[,1], y=pca_founder[,2]), color="black", size=5)+ # add the founder
  xlab("PCA 8.7%")+
  ylab("PCA 5.1%")+
  xlim(c(-12,30))+ylim(c(-13,43))+
  theme_classic()+
  theme(axis.line.x = element_blank())+
  theme(axis.line.y = element_blank())
ggsave(plot = pca_gen3_plot,filename = "figs/POP_EVOLUTION_PCA_2D_Xing_generation3.pdf",device = "pdf",width = 5,height = 3)

pca_gen2_plot <- df2 %>%
  ggplot(.) +
  geom_segment(aes(xend=pc1, yend=pc2,x=pca_founder[,1], y=pca_founder[,2],color=plotcolor))+
  scale_color_identity()+
  geom_point(aes(x=pc1, y=pc2), color="white", size=3.5)+
  geom_point(aes(x=pc1, y=pc2), color="black", size=3)+
  geom_point(aes(x=pc1, y=pc2, color=plotcolor),size=2.5)+
  geom_point(aes(x=pca_founder[,1], y=pca_founder[,2]), color="black", size=5)+ # add the founder
  xlab("PCA 8.7%")+
  ylab("PCA 5.1%")+
  xlim(c(-12,30))+ylim(c(-13,43))+
  theme_classic()+
  theme(axis.line.x = element_blank())+
  theme(axis.line.y = element_blank())
ggsave(plot = pca_gen2_plot,filename = "figs/POP_EVOLUTION_PCA_2D_Xing_generation2.pdf",device = "pdf",width = 5,height = 3)

pca_gen1_plot  <- df1 %>%
  ggplot(.) +
  geom_segment(aes(xend=pc1, yend=pc2,x=pca_founder[,1], y=pca_founder[,2],color=plotcolor))+
  scale_color_identity()+
  geom_point(aes(x=pc1, y=pc2), color="white", size=3.5)+
  geom_point(aes(x=pc1, y=pc2), color="black", size=3)+
  geom_point(aes(x=pc1, y=pc2, color=plotcolor),size=2.5)+
  geom_point(aes(x=pca_founder[,1], y=pca_founder[,2]), color="black", size=5)+ # add the founder
  xlab("PCA 8.7%")+
  ylab("PCA 5.1%")+
  xlim(c(-12,30))+ylim(c(-13,43))+
  theme_classic()+
  theme(axis.line.x = element_blank())+
  theme(axis.line.y = element_blank())
ggsave(plot = pca_gen1_plot,filename = "figs/POP_EVOLUTION_PCA_2D_Xing_generation1.pdf",device = "pdf",width = 5,height = 3)



## build the diffusion map for all generations using delta snp frequency which is calculated by comparing the
## current generation to p0
set.seed(1)
dm <- DiffusionMap(t(delta_p),verbose = T,n_pcs = 20,density_norm = T)
plot(dm)
dm_matrix <- tibble(site = meta$site,
                    plot = meta$plot,
                    generation = meta$generation,
                    sample = meta$sample_name,
                    DC1 = dm$DC1,
                    DC2 = dm$DC2) %>%
  left_join(.,df_climate_quantile,by="site")
dm_gen3 <- dm_matrix[meta$generation==3,]
dm_gen2 <- dm_matrix[meta$generation==2,]
dm_gen1 <- dm_matrix[meta$generation==1,]


dm_gen3_plot <- dm_gen3 %>%
  ggplot(.) +
  #geom_segment(aes(xend=DC1, yend=DC2,x=pca_founder[,1], y=pca_founder[,2],color=plotcolor))+
  scale_color_identity()+
  geom_point(aes(x=DC1, y=DC2), color="white", size=3.5)+
  geom_point(aes(x=DC1, y=DC2), color="black", size=3)+
  geom_point(aes(x=DC1, y=DC2, color=plotcolor),size=2.5)+
  #geom_point(aes(x=pca_founder[,1], y=pca_founder[,2]), color="black", size=5)+ # add the founder
  #xlab("PCA 8.7%")+
  #ylab("PCA 5.1%")+
  #xlim(c(-12,30))+ylim(c(-13,43))+
  theme_classic()+
  theme(axis.line.x = element_blank())+
  theme(axis.line.y = element_blank())
dm_gen3_plot
dm_gen2_plot <- dm_gen2 %>%
  ggplot(.) +
  #geom_segment(aes(xend=DC1, yend=DC2,x=pca_founder[,1], y=pca_founder[,2],color=plotcolor))+
  scale_color_identity()+
  geom_point(aes(x=DC1, y=DC2), color="white", size=3.5)+
  geom_point(aes(x=DC1, y=DC2), color="black", size=3)+
  geom_point(aes(x=DC1, y=DC2, color=plotcolor),size=2.5)+
  #geom_point(aes(x=pca_founder[,1], y=pca_founder[,2]), color="black", size=5)+ # add the founder
  #xlab("PCA 8.7%")+
  #ylab("PCA 5.1%")+
  #xlim(c(-12,30))+ylim(c(-13,43))+
  theme_classic()+
  theme(axis.line.x = element_blank())+
  theme(axis.line.y = element_blank())
dm_gen2_plot
dm_gen1_plot <- dm_gen1 %>%
  ggplot(.) +
  #geom_segment(aes(xend=DC1, yend=DC2,x=pca_founder[,1], y=pca_founder[,2],color=plotcolor))+
  scale_color_identity()+
  geom_point(aes(x=DC1, y=DC2), color="white", size=3.5)+
  geom_point(aes(x=DC1, y=DC2), color="black", size=3)+
  geom_point(aes(x=DC1, y=DC2, color=plotcolor),size=2.5)+
  #geom_point(aes(x=pca_founder[,1], y=pca_founder[,2]), color="black", size=5)+ # add the founder
  #xlab("PCA 8.7%")+
  #ylab("PCA 5.1%")+
  #xlim(c(-12,30))+ylim(c(-13,43))+
  theme_classic()+
  theme(axis.line.x = element_blank())+
  theme(axis.line.y = element_blank())
dm_gen1_plot


ggsave(plot = dm_gen1_plot,filename = "figs/POP_EVOLUTION_DM_2D_Xing_generation1.pdf",device = "pdf",width = 5,height = 3)
ggsave(plot = dm_gen2_plot,filename = "figs/POP_EVOLUTION_DM_2D_Xing_generation2.pdf",device = "pdf",width = 5,height = 3)
ggsave(plot = dm_gen3_plot,filename = "figs/POP_EVOLUTION_DM_2D_Xing_generation3.pdf",device = "pdf",width = 5,height = 3)


combined <- dm_gen1_plot | dm_gen2_plot | dm_gen3_plot

combined
ggsave(plot = combined,filename = "figs/POP_EVOLUTION_DM_2D_Xing_combined.pdf",device = "pdf",width = 8,height = 3)

















scatterplot3d(dm$DC1[meta$generation==1],dm$DC2[meta$generation==1],dm$DC3[meta$generation==1], xlab="DC1",ylab="DC2",zlab="DC3",cex.symbols = 1,grid=TRUE, box=FALSE,
              angle = 45,xlim = range(dm$DC1),ylim = range(dm$DC2),zlim = range(dm$DC3))
scatterplot3d(dm$DC1[meta$generation==2],dm$DC2[meta$generation==2],dm$DC3[meta$generation==2], xlab="DC1",ylab="DC2",zlab="DC3",cex.symbols = 1,grid=TRUE, box=FALSE,
              angle = 45,xlim = range(dm$DC1),ylim = range(dm$DC2),zlim = range(dm$DC3))
scatterplot3d(dm$DC1[meta$generation==3],dm$DC2[meta$generation==3],dm$DC3[meta$generation==3], xlab="DC1",ylab="DC2",zlab="DC3",cex.symbols = 1,grid=TRUE, box=FALSE,
              angle = 45,xlim = range(dm$DC1),ylim = range(dm$DC2),zlim = range(dm$DC3))



pca <- prcomp(t(delta_p))
plot(pca$x[meta$generation==1,1],pca$x[meta$generation==1,2],xlim=range(pca$x[,1]),ylim=range(pca$x[,2]))
plot(pca$x[meta$generation==2,1],pca$x[meta$generation==2,2],xlim=range(pca$x[,1]),ylim=range(pca$x[,2]))
plot(pca$x[meta$generation==3,1],pca$x[meta$generation==3,2],xlim=range(pca$x[,1]),ylim=range(pca$x[,2]))


scatterplot3d(pca$x[meta$generation==1,1],pca$x[meta$generation==1,2],pca$x[meta$generation==1,3], xlab="PC1",ylab="PC2",zlab="PC3",cex.symbols = 1,grid=TRUE, box=FALSE,
              angle = 45,xlim = range(pca$x[,1]),ylim = range(pca$x[,2]),zlim = range(pca$x[,3]))
scatterplot3d(pca$x[meta$generation==2,1],pca$x[meta$generation==2,2],pca$x[meta$generation==2,3], xlab="PC1",ylab="PC2",zlab="PC3",cex.symbols = 1,grid=TRUE, box=FALSE,
              angle = 45,xlim = range(pca$x[,1]),ylim = range(pca$x[,2]),zlim = range(pca$x[,3]))
scatterplot3d(pca$x[meta$generation==3,1],pca$x[meta$generation==3,2],pca$x[meta$generation==3,3], xlab="PC1",ylab="PC2",zlab="PC3",cex.symbols = 1,grid=TRUE, box=FALSE,
              angle = 45,xlim = range(pca$x[,1]),ylim = range(pca$x[,2]),zlim = range(pca$x[,3]))


## now color code the points based on the climate of the experimental sites (bio1 and bio12)

site_climate <- read.delim("data-external/bioclimvars_experimental_sites_era5.csv",sep=",")
site_climate <- site_climate %>%
  filter(site %in% unique(meta$site)) %>%
  arrange(match(site,unique(meta$site)))


quantile_custom <- function(x){
  q <- c()
  for(i in 1:length(x)){
    q <- c(q,round((x[i] - range(x)[1]) / (range(x)[2] - range(x)[1]) * 10))
  }
  q[q==0] <- 1
  return(q)
}

red_blue_colors <- colorRampPalette(c("#004080","white","#e50000"))(10)
plot(1:10,col=red_blue_colors,pch=16)


df_climate_quantile <- as.data.frame(matrix(nrow=nrow(site_climate),ncol=0)) %>%
  mutate(site = site_climate$site) %>%
  mutate(var1_quantile = quantile_custom(site_climate$bio1))


color <- c()
for(i in 1:nrow(df_climate_quantile)){
  color <- c(color,red_blue_colors[df_climate_quantile$var1_quantile[i]])
}
df_climate_quantile$color <- color


dm_df <- as.data.frame(matrix(nrow=nrow(meta),ncol = 0)) %>%
  mutate(site = meta$site) %>%
  mutate(generation = meta$generation) %>%
  mutate(DC1 = dm$DC1) %>%
  mutate(DC2 = dm$DC2) %>%
  mutate(DC3 = dm$DC3) %>%
  left_join(df_climate_quantile[,c(1,3)],by = "site") %>%
  left_join(site_climate[,c(1,2)],by = "site")


levels(dm_df$generation) <- c("year 1", "year 2", "year 3")
ggplot(dm_df,aes(x=DC1,y=DC2))+
  geom_point(color=dm_df$color,shape=16,cex=10)+
  facet_wrap(~generation,labeller = as_labeller(c(
    "1" = "Generation 1",
    "2" = "Generation 2",
    "3" = "Generation 3"
  )))+
  theme_minimal(base_size = 25,base_line_size = 0.2)
range(dm_df$DC1)
range(dm_df$DC2)
dm_df %>%
  filter(generation == 1 ) %>% 
  ggplot(.,aes(x=DC1,y=DC2))+
  geom_point(color=dm_df$color[dm_df$generation==1],shape=16,cex=10)+
  xlim(c(-0.05,0.3))+
  ylim(c(-0.4,0.02))+
  theme_minimal(base_size = 25,base_line_size = 0.2)




pdf(file = paste0("figs/",prefix,"diffusion_map_generation1.pdf"),useDingbats = F,width = 6,height = 5)
scatterplot3d(dm_df$DC1[dm_df$generation==1],dm_df$DC2[dm_df$generation==1],dm_df$DC3[dm_df$generation==1], 
              xlab="DC1",ylab="DC2",zlab="DC3",cex.symbols = 3,grid=F, box=FALSE,tick.marks = F,scale.y = T,
              angle = 60,xlim = range(dm_df$DC1),ylim = range(dm_df$DC2),zlim = range(dm_df$DC3),
              color = dm_df$color[dm_df$generation==1],pch=16)

dev.off()
pdf(file = paste0("figs/",prefix,"diffusion_map_generation2.pdf"),useDingbats = F,width = 6,height = 5)
scatterplot3d(dm_df$DC1[dm_df$generation==2],dm_df$DC2[dm_df$generation==2],dm_df$DC3[dm_df$generation==2], 
              xlab="DC1",ylab="DC2",zlab="DC3",cex.symbols = 3,grid=F, box=FALSE,tick.marks = F,scale.y = T,
              angle = 60,xlim = range(dm_df$DC1),ylim = range(dm_df$DC2),zlim = range(dm_df$DC3),
              color = dm_df$color[dm_df$generation==2],pch=16)
dev.off()
pdf(file = paste0("figs/",prefix,"diffusion_map_generation3.pdf"),useDingbats = F,width = 6,height = 5)
scatterplot3d(dm_df$DC1[dm_df$generation==3],dm_df$DC2[dm_df$generation==3],dm_df$DC3[dm_df$generation==3], 
              xlab="DC1",ylab="DC2",zlab="DC3",cex.symbols = 3,grid=F, box=FALSE,tick.marks = F,scale.y = T,
              angle = 60,xlim = range(dm_df$DC1),ylim = range(dm_df$DC2),zlim = range(dm_df$DC3),
              color = dm_df$color[dm_df$generation==3],pch = 16)
dev.off()

ggplot(site_climate,aes(x=site,y = bio1,color=bio1))+
  geom_point()+
  scale_color_gradient2(low="#004080", high="#e50000",mid = "white",midpoint = mean(site_climate$bio1))+
  theme(legend.key.size = unit(2, 'cm'))
ggsave(file = paste0("figs/",prefix,"diffusion_map_legend.pdf"),device = "pdf",height = 5,width = 5)




legend_image <- as.raster(matrix(colfunc(20), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 0, 0, 1,1)

### PCA colorred plot

pca_df <- as.data.frame(matrix(nrow=nrow(meta),ncol = 0)) %>%
  mutate(site = meta$site) %>%
  mutate(generation = meta$generation) %>%
  mutate(PC1 = pca$x[,1]) %>%
  mutate(PC2 = pca$x[,2]) %>%
  mutate(PC3 = pca$x[,3]) %>%
  left_join(df_climate_quantile[,c(1,3)],by = "site")

pdf(file = paste0("figs/",prefix,"population_change_pca_generation1.pdf"),useDingbats = F,width = 10,height = 8)
plot(pca_df$PC1[pca_df$generation==1],pca_df$PC2[pca_df$generation==1],xlim=range(pca_df$PC1),ylim=range(pca_df$PC2),col=pca_df$color[pca_df$generation==1],pch=16,cex=3)
dev.off()
pdf(file = paste0("figs/",prefix,"population_change_pca_generation2.pdf"),useDingbats = F,width = 10,height = 8)
plot(pca_df$PC1[pca_df$generation==2],pca_df$PC2[pca_df$generation==2],xlim=range(pca_df$PC1),ylim=range(pca_df$PC2),col=pca_df$color[pca_df$generation==2],pch=16,cex=3)
dev.off()
pdf(file = paste0("figs/",prefix,"population_change_pca_generation3.pdf"),useDingbats = F,width = 10,height = 8)
plot(pca_df$PC1[pca_df$generation==3],pca_df$PC2[pca_df$generation==3],xlim=range(pca_df$PC1),ylim=range(pca_df$PC2),col=pca_df$color[pca_df$generation==3],pch=16,cex=3)
dev.off()

scatterplot3d(pca_df$PC1[pca_df$generation==1],pca_df$PC2[pca_df$generation==1],pca_df$PC3[pca_df$generation==1], 
              xlab="PC1",ylab="PC2",zlab="PC3",cex.symbols = 1.5,grid=TRUE, box=FALSE,
              angle = 60,xlim = range(pca_df$PC1),ylim = range(pca_df$PC2),zlim = range(pca_df$PC3),
              color = pca_df$color[pca_df$generation==1],pch=16)

scatterplot3d(pca_df$PC1[pca_df$generation==2],pca_df$PC2[pca_df$generation==2],pca_df$PC3[pca_df$generation==2], 
              xlab="PC1",ylab="PC2",zlab="PC3",cex.symbols = 1.5,grid=TRUE, box=FALSE,
              angle = 60,xlim = range(pca_df$PC1),ylim = range(pca_df$PC2),zlim = range(pca_df$PC3),
              color = pca_df$color[pca_df$generation==2],pch=16)

scatterplot3d(pca_df$PC1[pca_df$generation==3],pca_df$PC2[pca_df$generation==3],pca_df$PC3[pca_df$generation==3], 
              xlab="PC1",ylab="PC2",zlab="PC3",cex.symbols = 1.5,grid=TRUE, box=FALSE,
              angle = 60,xlim = range(pca_df$PC1),ylim = range(pca_df$PC2),zlim = range(pca_df$PC3),
              color = pca_df$color[pca_df$generation==3],pch = 16)



