rm(list=ls())

library(ggplot2)
library(boot)
library(dplyr)
library(patchwork)

## local adaptation

## This script will explore any traits are correlated with Vs or Wmax

prefix <- "LOCAL_ADAPTATION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

stabilizing_parameters <- read.delim("data-intermediate/LOCAL_ADAPTATION_ecotype_Wmax_Vs.txt")


founder_climate <- read.delim("data-external/bioclimvars_ecotypes_era5.csv",sep=",") %>%
  filter(ecotype %in% stabilizing_parameters$ecotype) %>%
  arrange(match(ecotype,stabilizing_parameters$ecotype))

cor.test(stabilizing_parameters$Wmax,founder_climate$bio1_new)
df <- tibble(wmax = stabilizing_parameters$Wmax,
             bio1 = founder_climate$bio1_new)
ggplot(df,aes(x=wmax,y=bio1))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_minimal(base_size = 18)

ggsave(file = paste0("figs/",prefix,"Wmax_bio1_relationship"),device = "pdf",width = 5,height = 5)

## load phenotype data
load("data-external/atlas1001_phenotype_matrix_imputed_onlypheno_NEW.rda")

pheno_founder <- pheno %>%
  filter(id %in% stabilizing_parameters$ecotype) %>%
  arrange(match(id,stabilizing_parameters$ecotype))

stabilizing_parameters <- stabilizing_parameters %>%
  filter(ecotype %in% pheno_founder$id) %>%
  arrange(match(ecotype,pheno_founder$id))




table(stabilizing_parameters$ecotype == pheno_founder$id)

cor.test(stabilizing_parameters$Wmax,pheno_founder$FT10,method = "pearson")

correlation <- as.data.frame(matrix(ncol=3,nrow=1883))
colnames(correlation) <- c("phenotype","inv_Vs","Wmax")
correlation$phenotype <- colnames(pheno_founder)[2:1884]
for(i in 1:1883){
  correlation[i,2] <- cor(pheno_founder[,i+1],1/stabilizing_parameters$Vs,)
  correlation[i,3] <- cor(pheno_founder[,i+1],stabilizing_parameters$Wmax)
}

colnames(stabilizing_parameters)[1] <- "id"

combined <- pheno_founder %>%
  left_join(stabilizing_parameters,by = "id")


summary(lm(combined$Wmax~combined$Delta_13C))

p1 <- ggplot(combined,aes(x=Delta_13C,y=Wmax))+
  geom_point(pch=16,cex=3)+
  annotate(geom = "text",label="y = -5.903 - 0.194x \n p-value: 8.837e-07 R2:0.099",x=-34,y=2)+
  geom_smooth(method = "lm")+
  theme_classic()

summary(lm(combined$Wmax~combined$FT10))

p2 <- ggplot(combined,aes(x=FT10,y=Wmax))+
  geom_point(pch=16,cex=3)+
  annotate(geom = "text",label="y = 1.859947 - 0.009934x \n p-value: 2.854e-10 R2: 0.160",x=110,y=2)+
  geom_smooth(method = "lm")+
  theme_classic()

summary(lm(combined$Wmax~combined$X1_FT_Diameter_Field))

p3 <- ggplot(combined,aes(x=X1_FT_Diameter_Field,y=Wmax))+
  geom_point(pch=16,cex=3)+
  annotate(geom = "text",label="y = -0.692462 + 0.055064x \n p-value: 2.997e-11 R2: 0.176",x=25,y=2)+
  geom_smooth(method = "lm")+
  theme_classic()


p4 <- ggplot(combined,aes(x=drought_index,y=1/Vs))+
  geom_point(pch=16,cex=3)+
  annotate(geom = "text",label="y = 5.32239 - 515.55123x \n p-value: 0.0006959 R2: 0.046",x=5e-04,y=7.5)+
  geom_smooth(method = "lm",formula = "y~x")+
  theme_classic()

ggplot(combined,aes(x=FT10,y=1/Vs))+
  geom_point(pch=16,cex=3)+
  geom_smooth(method = "lm",formula = "y~x")+
  theme_classic()

ggplot(combined,aes(x=freezing_survival,y=1/Vs,color=1/Vs))+
  geom_point(pch=16,cex=3)+
  geom_smooth(method = "lm",formula = "y~x",colour = "grey")+
  scale_color_gradientn(colors=c(rev(RColorBrewer::brewer.pal(5,"Reds")),RColorBrewer::brewer.pal(5,"Greens")),
                        values = rescale(c(0,0.001,0.002,0.004,0.006,0.008,0.011,0.015,0.02)))+
  theme_classic(base_size = 18)
ggsave(file = paste0("figs/",prefix,"Vs_DSDS50_relationship.pdf"),device = "pdf",width = 6,height = 5)
cor.test(combined$DSDS50.1,1/combined$Vs)


p <- (p1 | p2) / (p3 | p4)
ggsave(g,file = paste0("figs/",prefix,"Wmax_Vs_relationship"),device = "pdf",width = 5,height = 5)
