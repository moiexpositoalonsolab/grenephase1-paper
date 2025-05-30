rm(list=ls())

library(ggplot2)
library(boot)
library(dplyr)
library(patchwork)
library(zoo)


## local adaptation

## This script will explore the LD block change overtime in GrENE-net

prefix <- "POP_EVOLUTION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

meta <- read.delim("data/merged_sample_table.csv",sep = ",")
terminal_index <- read.delim("data/site_plot_terminal_index.txt")

site_name <- read.delim("data/site_info.txt",check.names = F,encoding = "Latin-1")


## first calcualte the LD blocks from the founder population

founder_0.2 <- read.delim("LD_across_generations/LDblocks_0.2_founder.txt",header=F)
founder_0.2$len <- (founder_0.2$V3-founder_0.2$V2) / 1000
founder_0.2$mid <- (founder_0.2$V3+founder_0.2$V2)/2 / 10^6

founder_0.5 <- read.delim("LD_across_generations/LDblocks_0.5_founder.txt",header=F)
founder_0.5$len <- (founder_0.5$V3-founder_0.5$V2) / 1000
founder_0.5$mid <- (founder_0.5$V3+founder_0.5$V2)/2 / 10^6

founder_0.8 <- read.delim("LD_across_generations/LDblocks_0.8_founder.txt",header=F)
founder_0.8$len <- (founder_0.8$V3-founder_0.8$V2) / 1000
founder_0.8$mid <- (founder_0.8$V3+founder_0.8$V2)/2 / 10^6



## calculate mean LD block length for GrENE-net samples
## all 3 generations
index_gen_3 <-  which(terminal_index$generation1==1 &terminal_index$generation2==2 & terminal_index$generation3==3)


site_plot_index <- c(terminal_index$index1[index_gen_3],terminal_index$index2[index_gen_3],terminal_index$index3[index_gen_3])
all_3_generations_meta <- meta[site_plot_index,]

len0.2 <- c()
len0.5 <- c()
len0.8 <- c()

count0.2 <- c()
count0.5 <- c()
count0.8 <- c()


for (i in 1:nrow(all_3_generations_meta)){
    name <- paste0("LD_across_generations/LDblocks_",0.2,"_",all_3_generations_meta$sample_name[i],".txt")
    tmp <- read.delim(name,header=F)
    len0.2 <- c(len0.2,median(tmp$V3-tmp$V2)/1000)
    count0.2 <- c(count0.2,nrow(tmp))
    
    name <- paste0("LD_across_generations/LDblocks_",0.5,"_",all_3_generations_meta$sample_name[i],".txt")
    tmp <- read.delim(name,header=F)
    len0.5 <- c(len0.5,median(tmp$V3-tmp$V2)/1000)
    count0.5 <- c(count0.8,nrow(tmp))
    
    name <- paste0("LD_across_generations/LDblocks_",0.8,"_",all_3_generations_meta$sample_name[i],".txt")
    tmp <- read.delim(name,header=F)
    len0.8 <- c(len0.8,median(tmp$V3-tmp$V2)/1000)
    count0.8 <- c(count0.8,nrow(tmp))
}

block_summary <- rbind(all_3_generations_meta,all_3_generations_meta,all_3_generations_meta)
block_summary$LD <- c(rep(0.2,426),rep(0.5,426),rep(0.8,426))
block_summary$block_length <- c(len0.2,len0.5,len0.8)
block_summary$count <- c(count0.2,count0.5,count0.8)
block_summary$generation <- as.factor(block_summary$generation)

len_segments <- data.frame(
  xstart = c(1 - 0.3, 2 - 0.3, 3 - 0.3),
  xend   = c(1 + 0.3, 2 + 0.3, 3 + 0.3),
  y      = c(median(founder_0.2$len), median(founder_0.5$len), median(founder_0.8$len))  # the horizontal level you want
)


p1 <- ggplot()+
  geom_boxplot(data=block_summary,aes(x=factor(LD),y=block_length,fill=generation))+
  geom_segment(data=len_segments,aes(x=xstart,xend =xend,y=y,yend=y),col="red",lwd=0.8)+
  xlab("LD (r2)")+
  ylab("Median haplotype block length \n per replicate (kb)")+
  scale_fill_brewer(palette="Dark2")+
  theme_classic(base_size = 18)


p2 <- ggplot()+
  geom_boxplot(data=block_summary[block_summary$LD==0.2,],aes(x=factor(site),y=block_length,fill=generation))+
  geom_hline(yintercept = median(founder_0.2$len),col="red",lwd=0.8)+
  xlab("GrENE-net gardens (with 3 generations of data)")+
  ylab("Median haplotype block length \n per replicate (kb, r2=0.2)")+
  scale_fill_brewer(palette="Dark2")+
  theme_classic(base_size = 18)


p3 <- ggplot()+
  geom_boxplot(data=block_summary[block_summary$LD==0.5,],aes(x=factor(site),y=block_length,fill=generation))+
  geom_hline(yintercept = median(founder_0.5$len),col="red",lwd=0.8)+
  xlab("GrENE-net gardens (with 3 generations of data)")+
  ylab("Median haplotype block length \n per replicate (kb, r2=0.5)")+
  scale_fill_brewer(palette="Dark2")+
  theme_classic(base_size = 18)
p4 <- ggplot()+
  geom_boxplot(data=block_summary[block_summary$LD==0.8,],aes(x=factor(site),y=block_length,fill=generation))+
  geom_hline(yintercept = median(founder_0.8$len),col="red",lwd=0.8)+
  xlab("GrENE-net gardens (with 3 generations of data)")+
  ylab("Median haplotype block length \n per replicate (kb, r2=0.8)")+
  scale_fill_brewer(palette="Dark2")+
  theme_classic(base_size = 18)

p2 / p3 /p4
p1 / p3
ggsave(plot = p1/p3,filename =  paste0("figs/",prefix,"LD_across_generations.pdf"),device = "pdf",width = 10,height = 8)


median(block_summary$block_length[block_summary$LD==0.5])
median(founder_0.5$len)
block_median_1 <- block_summary %>%
  filter(LD==0.5 & generation==1) %>%
  group_by(site) %>%
  summarise(median(block_length))
mean(block_median_1$`median(block_length)`)

block_median_2 <- block_summary %>%
  filter(LD==0.5 & generation==2) %>%
  group_by(site) %>%
  summarise(median(block_length))
mean(block_median_2$`median(block_length)`)
quantile(block_median_2$`median(block_length)`)

block_median_3 <- block_summary %>%
  filter(LD==0.5 & generation==3) %>%
  group_by(site) %>%
  summarise(median(block_length))
mean(block_median_3$`median(block_length)`)
quantile(block_median_3$`median(block_length)`)

table(block_median_1$`median(block_length)` > median(founder_0.5$len))
# block_1 <- read.delim("LD_across_generations/LDblocks_0.2_4_1_3.txt",header=F)
# block_1$len <- (block_1$V3-block_1$V2) / 1000
# block_1$mid <- (block_1$V3+block_1$V2)/2 / 10^6
# 
# block_2 <- read.delim("LD_across_generations/LDblocks_0.2_4_2_3.txt",header=F)
# block_2$len <- (block_2$V3-block_2$V2)/1000
# block_2$mid <- (block_2$V3+block_2$V2)/2 /10^6
# 
# block_3 <- read.delim("LD_across_generations/LDblocks_0.2_4_3_3.txt",header=F)
# block_3$len <- (block_3$V3-block_3$V2) / 1000
# block_3$mid <- (block_3$V3+block_3$V2)/2 / 10^6
# 
# ggplot(block_1,aes(x=mid,y=len))+
#   geom_line(aes(y=rollmean(len, 5, na.pad=T)),cex=1,lty=1,col="red",alpha=0.5)+
#   geom_line(data =block_2, aes(y=rollmean(len, 5, na.pad=T)),cex=1,lty=1,col="blue",alpha=0.5)+
#   geom_line(data =block_3, aes(y=rollmean(len, 5, na.pad=T)),cex=1,lty=1,col="gold",alpha=0.5)+
#   geom_line(data =founder_0.2, aes(y=rollmean(len, 5, na.pad=T)),cex=1,lty=1,col="black",alpha=0.5)+
#   theme_classic()





