rm(list=ls())

library(ggplot2)
library(corrplot)
library(MCMCglmm)
library(r2glmm)
library(genio)
library(MASS)
library(grDevices)
library(dplyr)
library(patchwork)

## local adaptation

## This script fit the stabilizing selection model using all three years data and then focus on generation 1
## and see which bio variables are most explaintory and should be included in the model

prefix <- "LOCAL_ADAPTATION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

## load the era5 stabilizing selection table
mydata <- read.delim("data-intermediate/LOCAL_ADAPTATION_stabilizing_selection_data_2024Jun18_collectionsite_era5_allyears_pt1_pt0.txt",header = T)
mydata$ecotype <- factor(mydata$ecotype,levels = unique(mydata$ecotype))
mydata$site <- factor(mydata$site,levels = unique(mydata$site))
unique_sites <- levels(mydata$site)

## lets first build the model with ALL generations
df <- mydata %>%
  mutate(log_p_ratio = log(pt1/pt0)) %>% 
  filter(!is.infinite(log_p_ratio) & !is.na(log_p_ratio)) %>%
  mutate(bio1_diff_sq = (bio1_ecotype - bio1_site)^2)

plot(df$bio1_ecotype-df$bio1_site,df$pt1-df$pt0)

m1 <- MCMCglmm(log_p_ratio ~ 1 + bio1_diff_sq, random = ~ site + us(1+bio1_diff_sq):ecotype, data = df,
                 family = "gaussian", verbose = T,nitt = 20000,burnin = 3000, pr = TRUE, saveX = TRUE, saveZ = TRUE,thin = 5)
posterior_m1 <- apply(m1$Sol,2,mean)
prd_m1 <- as.numeric(m1$X %*% posterior_m1[1:2] + m1$Z %*% posterior_m1[3:495])
wmax_m1 <- as.numeric(m1$Z[,32:262] %*% posterior_m1[34:264])
vs_m1 <- as.numeric(m1$X[,2] * posterior_m1[2] +  m1$Z[,263:493] %*% posterior_m1[265:495])
wavg_m1 <- as.numeric(m1$Z[,1:31] %*% posterior_m1[3:33])
summary(lm(df$log_p_ratio~wavg_m1))
summary(lm(df$log_p_ratio~vs_m1))
summary(lm(df$log_p_ratio~wmax_m1))
summary(lm(df$log_p_ratio~prd_m1))

## the results are pretty bad indicating selection is not the only / main driver for the change in the population composition


## build the stabilizing selection model for each generation on bio1 only

## generation 1
df1 <- mydata %>%
  filter(generation==1) %>% 
  mutate(log_p_ratio = log(pt1/pt0)) %>% 
  filter(!is.infinite(log_p_ratio) & !is.na(log_p_ratio)) %>%
  mutate(bio1_diff_sq = (bio1_ecotype - bio1_site)^2)

plot(df1$bio1_ecotype-df1$bio1_site,df1$pt1-df1$pt0,ylim = c(-1,1))

m1.1 <- MCMCglmm(log_p_ratio ~ 1 + bio1_diff_sq, random = ~ site + ecotype + us(bio1_diff_sq):ecotype, data = df1,
                 family = "gaussian", verbose = T,nitt = 20000,burnin = 3000, pr = TRUE, saveX = TRUE, saveZ = TRUE,thin = 5)
posterior_m1.1 <- apply(m1.1$Sol,2,mean)
prd_m1.1 <- as.numeric(m1.1$X %*% posterior_m1.1[1:2] + m1.1$Z %*% posterior_m1.1[3:495])
wmax_m1.1 <- as.numeric(m1.1$Z[,32:262] %*% posterior_m1.1[34:264])
vs_m1.1 <- as.numeric(m1.1$X[,2] * posterior_m1.1[2] +  m1.1$Z[,263:493] %*% posterior_m1.1[265:495])
wavg_m1.1 <- as.numeric(m1.1$Z[,1:31] %*% posterior_m1.1[3:33])
summary(lm(df1$log_p_ratio~wavg_m1.1))
summary(lm(df1$log_p_ratio~vs_m1.1))
summary(lm(df1$log_p_ratio~wmax_m1.1))
summary(lm(df1$log_p_ratio~prd_m1.1))



df1 <- df1 %>%
  mutate(prd = prd_m1.1)

R2 <- c()
for(i in unique(df1$ecotype)){
  tmp <- df1[df1$ecotype == i,]
  R2 <- c(R2,summary(lm(tmp$log_p_ratio~tmp$prd))$adj.r.squared)
}
sd(R2)

p1 <- ggplot(df1,aes(x=bio1_ecotype-bio1_site,y=pt1-pt0))+
  geom_point()+
  ylim(c(-1,1))+
  xlim(c(-35,20))+
  ylab("ecotype frequency change \n (generation 1 - generation 0") +
  xlab("climatic different (Annual Mean Temp)")+
  annotate("text",x=-25,y=1,label="stabilizing selection model \n R2: 0.3389") +
  #geom_vline(xintercept = 0,colour = "red")+
  theme_classic()
p1




## generation 2
df2 <- mydata %>%
  filter(generation==2) %>% 
  mutate(log_p_ratio = log(pt1/pt0)) %>% 
  filter(!is.infinite(log_p_ratio) & !is.na(log_p_ratio)) %>%
  mutate(bio1_diff_sq = (bio1_ecotype - bio1_site)^2)

plot(df2$bio1_ecotype-df2$bio1_site,df2$pt1-df2$pt0,ylim = c(-1,1))

m1.2 <- MCMCglmm(log_p_ratio ~ 1 + bio1_diff_sq, random = ~ site + ecotype+ us(bio1_diff_sq):ecotype, data = df2,
                 family = "gaussian", verbose = T,nitt = 20000,burnin = 3000, pr = TRUE, saveX = TRUE, saveZ = TRUE,thin = 5)
posterior_m1.2 <- apply(m1.2$Sol,2,mean)
prd_m1.2 <- as.numeric(m1.2$X %*% posterior_m1.2[1:2] + m1.2$Z %*% posterior_m1.2[3:489])
wmax_m1.2 <- as.numeric(m1.2$Z[,26:256] %*% posterior_m1.2[28:258])
vs_m1.2 <- as.numeric(m1.2$X[,2] * posterior_m1.2[2] +  m1.2$Z[,257:487] %*% posterior_m1.2[259:489])
wavg_m1.2 <- as.numeric(m1.2$Z[,1:25] %*% posterior_m1.2[3:27])
summary(lm(df2$log_p_ratio~wavg_m1.2))
summary(lm(df2$log_p_ratio~vs_m1.2))
summary(lm(df2$log_p_ratio~wmax_m1.2))
summary(lm(df2$log_p_ratio~prd_m1.2))

p2 <- ggplot(df2,aes(x=bio1_ecotype-bio1_site,y=pt1-pt0))+
  geom_point()+
  ylim(c(-1,1))+
  xlim(c(-35,20))+
  ylab("ecotype frequency change \n (generation 2 - generation 1") +
  xlab("climatic different (Annual Mean Temp)")+
  annotate("text",x=-25,y=1,label="stabilizing selection model \n R2: 0.1774") +
  #geom_vline(xintercept = 0,colour = "red")+
  theme_classic()
p2


## generation 3
df3 <- mydata %>%
  filter(generation==3) %>% 
  mutate(log_p_ratio = log(pt1/pt0)) %>% 
  filter(!is.infinite(log_p_ratio) & !is.na(log_p_ratio)) %>%
  mutate(bio1_diff_sq = (bio1_ecotype - bio1_site)^2)

plot(df3$bio1_ecotype-df3$bio1_site,df3$pt1-df3$pt0,ylim = c(-1,1))

m1.3 <- MCMCglmm(log_p_ratio ~ 1 + bio1_diff_sq, random = ~ site + ecotype + us(bio1_diff_sq):ecotype, data = df3,
                 family = "gaussian", verbose = T,nitt = 20000,burnin = 3000, pr = TRUE, saveX = TRUE, saveZ = TRUE,thin = 5)
posterior_m1.3 <- apply(m1.3$Sol,2,mean)
prd_m1.3 <- as.numeric(m1.3$X %*% posterior_m1.3[1:2] + m1.3$Z %*% posterior_m1.3[3:483])
wmax_m1.3 <- as.numeric(m1.3$Z[,20:250] %*% posterior_m1.3[22:252])
vs_m1.3 <- as.numeric(m1.3$X[,2] * posterior_m1.3[2] +  m1.3$Z[,251:481] %*% posterior_m1.3[253:483])
wavg_m1.3 <- as.numeric(m1.3$Z[,1:19] %*% posterior_m1.3[3:21])
summary(lm(df3$log_p_ratio~wavg_m1.3))
summary(lm(df3$log_p_ratio~vs_m1.3))
summary(lm(df3$log_p_ratio~wmax_m1.3))
summary(lm(df3$log_p_ratio~prd_m1.3))

p3 <- ggplot(df3,aes(x=bio1_ecotype-bio1_site,y=pt1-pt0))+
  geom_point()+
  ylim(c(-1,1))+
  xlim(c(-35,20))+
  ylab("ecotype frequency change \n (generation 3 - generation 2") +
  xlab("climatic different (Annual Mean Temp)")+
  annotate("text",x=-25,y=1,label="stabilizing selection model \n R2: 0.1366") +
  #geom_vline(xintercept = 0,colour = "red")+
  theme_classic()
p3

combined_stabilizing_bio1 <- p1 | p2 | p3
combined_stabilizing_bio1
ggsave(plot = combined_stabilizing_bio1,filename =  paste0("figs/",prefix,"stabilizing_selection_per_generation_bio1.pdf"),device = "pdf",width = 20,height = 10,units = "in")



## now run stabilizing selection on the generation1 data using each bio variaible 
## and plot to see which bio varaible explains the highest

df_bio <- mydata %>%
  filter(generation==1) %>% 
  filter(site != 33) %>%
  mutate(log_p_ratio = log(pt1/pt0)) %>% 
  filter(!is.infinite(log_p_ratio) & !is.na(log_p_ratio))

for ( i in 13:19){
  if(i %in% c(4,12,13,14,15,16,17,18,19)){
    df_bio$diff_sq <- (df_bio[,i+6]/10 - df_bio[,i+6+21]/10)^2
    df_bio$diff <- df_bio[,i+6]/10 - df_bio[,i+6+21]/10
  } else{
    df_bio$diff_sq <- (df_bio[,i+6] - df_bio[,i+6+21])^2
    df_bio$diff <- (df_bio[,i+6] - df_bio[,i+6+21])
  }
  
  m <- MCMCglmm(log_p_ratio ~ 1 + diff_sq, random = ~ site + ecotype + us(diff_sq):ecotype, data = df_bio,
                   family = "gaussian", verbose = T,nitt = 20000,burnin = 3000, pr = TRUE, saveX = TRUE, saveZ = TRUE,thin = 5)
  posterior_m <- apply(m$Sol,2,mean)
  prd_m <- as.numeric(m$X %*% posterior_m[1:2] + m$Z %*% posterior_m[3:494])
  R2 <- summary(lm(df_bio$log_p_ratio~prd_m))$adj.r.squared
  
  vs_m <- as.numeric(m$X[,2] * posterior_m[2] +  m$Z[,262:492] %*% posterior_m[264:494])
  R2_vs <- summary(lm(df_bio$log_p_ratio~vs_m))$adj.r.squared

  p <- ggplot(df_bio,aes(x=diff,y=pt1-pt0))+
    geom_point(cex=4)+
    xlim(range(df_bio$diff))+
    ylab("ecotype frequency change") +
    xlab(paste0("climatic difference (BIO",i,")"))+
    #annotate("text",x=10,y=1,label=paste0("R2: ",round(R2,4)))
    ggtitle(paste0("R2_model: ",round(R2,4),"\nR2_Vs:",round(R2_vs,4) ))+
    theme_classic(base_size = 30)
  ggsave(p,filename =  paste0("figs/stabilizing_selection_per_bio/",prefix,"stabilizing_selection_generation_1_bio",i,".pdf"),device = "pdf",width = 7,height = 6,units = "in")
  ggsave(p,filename =  paste0("figs/stabilizing_selection_per_bio/",prefix,"stabilizing_selection_generation_1_bio",i,".png"),device = "png",width = 7,height = 6,dpi = 800)
  
}






