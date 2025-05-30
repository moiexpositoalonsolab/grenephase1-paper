rm(list=ls())

library(dplyr)

## local adaptation

## This script generates the vs and wmax phenotypes for 1001g genomes for BSLMM

prefix <- "FUTURE_PREDICTION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

vcf_id <- read.delim("data-external/1001g_grenet_accessions_climate.txt",header=F)$V1
info <- read.delim("data/1001g_regmap_grenet_ecotype_info_corrected_2024May16.csv",sep=",")

## first check the ids are matched
table(info$ecotype_id[info$GrENE.net==T | info$source=="1TG"] %in% vcf_id)
stabilizing_selection <- read.delim("data-intermediate/LOCAL_ADAPTATION_ecotype_Wmax_Vs.txt")

set.seed(100)
validation <- sample(stabilizing_selection$ecotype,30,replace = F)
validation_df <- stabilizing_selection[stabilizing_selection$ecotype %in% validation,]

hist(validation_df$Wmax)
hist(stabilizing_selection$Wmax)
hist(validation_df$Vs)
hist(stabilizing_selection$Vs)

stabilizing_selection_training <- stabilizing_selection[! stabilizing_selection$ecotype %in% validation,]

df <- as.data.frame(matrix("NA",ncol=3,nrow=length(vcf_id)))
colnames(df) <- c("ecotype","Wmax","Vs")
df$ecotype <- vcf_id

for(i in 1:nrow(stabilizing_selection_training)){
  index <- which(df$ecotype == stabilizing_selection_training$ecotype[i])
  df[index,2] <- stabilizing_selection_training$Wmax[i]
  df[index,3] <- stabilizing_selection_training$Vs[i]
}

#write.table(df,"data-intermediate/FUTURE_PREDICTION_1001g_Wmax_Vs.txt",sep="\t",append = F,quote = F,col.names = F,row.names = F)



Vs_prd <- read.delim("data-intermediate/Vs_prd.prdt.txt",header = F)
Wmax_prd <- read.delim("data-intermediate/Wmax_prd.prdt.txt",header = F)

df$Vs_prd <- Vs_prd$V1
df$Wmax_prd <- Wmax_prd$V1

df_prd <- df[df$ecotype %in% validation,]
validation_df <- validation_df %>%
  left_join(.,df_prd[,c(1,4,5)],by="ecotype")



wmax_summary <- summary(lm(Wmax~Wmax_prd,data = validation_df))
vs_summary <- summary(lm(Vs~Vs_prd,data = validation_df))


plot(df_prd$Wmax_prd+wmax_summary$coefficients[1,1],validation_df$Wmax)
abline(a=0,b=1)
text("R-squared 0.4823",x=1.8,y=1)

plot(df_prd$Vs_prd+vs_summary$coefficients[1,1],validation_df$Vs)
abline(a=0,b=1)
text("R-squared 0.5598",x=700,y=400)



df$Vs_prd_adj <- df$Vs_prd*vs_summary$coefficients[2,1] + vs_summary$coefficients[1,1]
df$Wmax_prd_adj <- df$Wmax_prd*wmax_summary$coefficients[2,1] + wmax_summary$coefficients[1,1]


for(i in 1:nrow(df)){
  if(is.na(df$Vs_prd_adj[i])){
    df$Vs_prd_adj[i] <- as.numeric(df$Vs[i])
    df$Wmax_prd_adj[i] <- as.numeric(df$Wmax[i])
  }
}

for(i in 1:nrow(df)){
  if(df$Vs_prd_adj[i] < 0){
    df$Vs_prd_adj[i] <- vs_summary$coefficients[1,1]
  }
}

write.table(df[,c(1,6,7)],"data-intermediate/FUTURE_PREDICTION_1001g_grenet_predicted_Wmax_Vs.txt")




