rm(list=ls())



## This script will evaluate HapFIRE performance at the snp and ecotype level frequency
## reconstruction, and compare it with HAFpipe

library(gdsfmt)
library(SNPRelate)
library(SeqArray)

prefix <- "EXPERIMENTAL_INTRO_"

setwd("~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/data-intermediate/hapfire_performance/")

## 10 x 5 x 2

batch <- paste0("batch",1:10)
num_ecotype <- c(2,5,20,50,150)
depth <- c("1x","10x")

## load the VCF file and convert it to snp matrix
seqVCF2GDS("200_test_chr1.recode.vcf", "200_test_chr1.gds")

genofile <- seqOpen("200_test_chr1.gds")
ecotype_id <- seqGetData(genofile, "sample.id")
snp_matrix <- as.matrix(seqGetData(genofile, "$dosage"))
variant.id <- paste0(1,"_",seqGetData(genofile, "position"))
seqClose(genofile)

dim(snp_matrix)

## HapFIRE

hapFIRE_snp <- as.data.frame(matrix(nrow=100,ncol= 8))
colnames(hapFIRE_snp) <- c("sample","batch","num_ecotype","depth","r2","rmse","mape","source")

hapFIRE_ecotype <- as.data.frame(matrix(nrow=100,ncol= 8))
colnames(hapFIRE_ecotype) <- c("sample","batch","num_ecotype","depth","r2","nrmse","mape","source")

hafpipe_snp <- as.data.frame(matrix(nrow=100,ncol= 8))
colnames(hafpipe_snp) <- c("sample","batch","num_ecotype","depth","r2","rmse","mape","source")


m = 1

for(i in batch){
  for(j in num_ecotype){
    for(k in depth){
      sample_name <- paste(i,j,k,sep = "_")
      ## ground truth
      true_ecotypes <- read.delim(paste0("bam_testing/",i,"/",j,".list"),header = F)
      true_ecotypes <- as.numeric(gsub(".*alignment/(.*).bam","\\1",true_ecotypes$V1))
      snp_matrix_true <- snp_matrix[ecotype_id %in% true_ecotypes,]
      true_snp_freq <- 1-apply(snp_matrix_true,2,function(x) sum(x)/j/2)
      true_ecotypes_freq <- rep(0,200)
      true_ecotypes_freq[match(true_ecotypes,ecotype_id)] <- 1
      true_ecotypes_freq <- true_ecotypes_freq / j

      ## hapFIRE_snp
      hapFIRE_snp[m,1] <- sample_name
      hapFIRE_snp[m,2] <- i
      hapFIRE_snp[m,3] <- j
      hapFIRE_snp[m,4] <- k
      df <- read.delim(paste0("hapFIRE/",sample_name,"_snp_frequency.txt"),header=F)
      df_1 <- df[df$V1==1,]
      hapFIRE_snp[m,5] <- summary(lm(true_snp_freq~df_1$V3))$adj.r.squared
      hapFIRE_snp[m,6] <- sqrt(mean((true_snp_freq-df_1$V3)^2))
      norm_diff <- (true_snp_freq-df_1$V3)/true_snp_freq
      hapFIRE_snp[m,7] <- median(abs(norm_diff[which(! (is.infinite(norm_diff) |is.na(norm_diff)))]))

      # hapFIRE_ecotype
      hapFIRE_ecotype[m,1] <- sample_name
      hapFIRE_ecotype[m,2] <- i
      hapFIRE_ecotype[m,3] <- j
      hapFIRE_ecotype[m,4] <- k
      df_ecotype <- read.delim(paste0("hapFIRE/",sample_name,"_ecotype_frequency.txt"),header=F)
      hapFIRE_ecotype[m,5] <- summary(lm(true_ecotypes_freq~df_ecotype$V2))$adj.r.squared
      hapFIRE_ecotype[m,6] <- sqrt(mean((true_ecotypes_freq-df_ecotype$V2)^2)) / sd(true_ecotypes_freq)
      norm_diff <- (true_ecotypes_freq-df_ecotype$V2)/true_ecotypes_freq
      hapFIRE_ecotype[m,7] <- median(abs(norm_diff[which(! (is.infinite(norm_diff) |is.na(norm_diff)))]))


      ## hafpipe_snp
      hafpipe_snp[m,1] <- sample_name
      hafpipe_snp[m,2] <- i
      hafpipe_snp[m,3] <- j
      hafpipe_snp[m,4] <- k
      df_hafpipe <- read.delim(paste0("haf-pipe/",sample_name,"/",j,"_",k,".bam.1.afSite"),header=T,sep=",")
      hafpipe_snp[m,5] <- summary(lm(true_snp_freq[variant.id %in% paste(1,df_hafpipe$pos,sep = "_")]~df_hafpipe$af))$adj.r.squared
      hafpipe_snp[m,6] <- sqrt(mean((true_snp_freq[variant.id %in% paste(1,df_hafpipe$pos,sep = "_")]-df_hafpipe$af)^2))
      norm_diff <- (true_snp_freq[variant.id %in% paste(1,df_hafpipe$pos,sep = "_")]-df_hafpipe$af)/true_snp_freq[variant.id %in% paste(1,df_hafpipe$pos,sep = "_")]
      hafpipe_snp[m,7] <- median(abs(norm_diff[which(! (is.infinite(norm_diff) |is.na(norm_diff)))]))

      m <- m+1
      print(m)
    }
  }
}

hapFIRE_snp$source <- "hapFIRE"
hafpipe_snp$source <- "HAFpipe"
hapFIRE_ecotype$source <- "hapFIRE"


write.table(hapFIRE_snp,"hapFIRE_snp_summary.txt",append = F,quote = F,sep = "\t",row.names = F)
write.table(hapFIRE_ecotype,"hapFIRE_ecotype_summary.txt",append = F,quote = F,sep = "\t",row.names = F)
write.table(hafpipe_snp,"hafpipe_snp_summary.txt",append = F,quote = F,sep = "\t",row.names = F)

hapFIRE_snp %>%
  group_by(num_ecotype,depth) %>%
  summarise(mean_RMSE = mean(rmse),
            quantile_25 = quantile(rmse,0.25),
            quantile_75 = quantile(rmse,0.75),
            mean_r2 = mean(r2),
            quantile_25_r2 = quantile(r2,0.25),
            quantile_75_r2 = quantile(r2,0.75))

hafpipe_snp %>%
  group_by(num_ecotype,depth) %>%
  summarise(mean_RMSE = mean(rmse),
            quantile_25 = quantile(rmse,0.25),
            quantile_75 = quantile(rmse,0.75),
            mean_r2 = mean(r2),
            quantile_25_r2 = quantile(r2,0.25),
            quantile_75_r2 = quantile(r2,0.75))

hapFIRE_ecotype %>%
  group_by(num_ecotype,depth) %>%
  summarise(mean_RMSE = mean(nrmse),
            quantile_25 = quantile(nrmse,0.25),
            quantile_75 = quantile(nrmse,0.75),
            mean_r2 = mean(r2),
            quantile_25_r2 = quantile(r2,0.25),
            quantile_75_r2 = quantile(r2,0.75))


p1 <- ggplot(hapFIRE_ecotype,aes(x=as.factor(num_ecotype),y = r2,color=depth))+
  geom_boxplot()+
  xlab("number of pooled accessions")+
  ylab("R2")+
  theme_minimal(base_size = 18)
p1
p2 <- ggplot(hapFIRE_ecotype,aes(x=as.factor(num_ecotype),y = nrmse,color=depth))+
  geom_boxplot()+
  xlab("number of pooled accessions")+
  ylab("Normalized Rooted \n Mean Squared Error")+
  theme_minimal(base_size = 18)
p2
ggsave(plot = p1 | p2,filename = paste0("../../figs/",prefix,"HapFIRE_ecotype_frequency_evaluation.pdf"),device = "pdf",width = 10,height = 5,units = "in")


hapFIRE_snp$source <- "hapFIRE"
hafpipe_snp$source <- "HAFpipe"

combined <- rbind(hapFIRE_snp,hafpipe_snp)


p3 <- combined %>%
  filter(depth=="1x") %>%
  ggplot(.,aes(x=as.factor(num_ecotype),y=r2,color=source))+
  geom_boxplot()+
  xlab("number of pooled accessions \n (1x depth)")+
  ylab("R2")+
  theme_minimal(base_size = 18)

p4 <- combined %>%
  filter(depth=="1x") %>%
  ggplot(.,aes(x=as.factor(num_ecotype),y=rmse,color=source))+
  geom_boxplot()+
  xlab("number of pooled accessions \n (1x depth)")+
  ylab("Root Mean Squared Error")+
  theme_minimal(base_size = 18)
p5 <- combined %>%
  filter(depth=="10x") %>%
  ggplot(.,aes(x=as.factor(num_ecotype),y=r2,color=source))+
  geom_boxplot()+
  xlab("number of pooled accessions \n (10x depth)")+
  ylab("R2")+
  theme_minimal(base_size = 18)
p6 <- combined %>%
  filter(depth=="10x") %>%
  ggplot(.,aes(x=as.factor(num_ecotype),y=rmse,color=source))+
  geom_boxplot()+
  xlab("number of pooled accessions \n (10x depth)")+
  ylab("Root Mean Squared Error")+
  theme_minimal(base_size = 18)
p3 | p4
p5 | p6
