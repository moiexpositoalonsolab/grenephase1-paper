rm(list=ls())

## POPULATION EVOLUTION

## This script will use the snp selection LRT1 results to generation a selection atlas for
## the entire grenet experiment. I will use both SNP and haplotypes as units for the atlas

prefix <- "POP_EVOLUTION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)
source("scripts/manhattan_plot_Xing.R")

af <- read.delim("data-intermediate/snp_allele_frequency.txt",header=F)
locations <- read.delim("data-intermediate/chromosome_locations.txt",header=F)
locations <- locations[af$V1>0.05,]
site_list <- read.delim("data-intermediate/list.txt",header=F)



### SNP-based ####

pvalues <- as.data.frame(matrix(nrow=1174429,ncol=nrow(site_list)))

for(i in 1:length(site_list$V1)){
  name <- site_list[i,]
  pvalue <- read.delim(paste("data-intermediate/snp_lrt1/gen1_",name,"_pvalue.txt",sep = ""),header=F)
  pvalue$V1[pvalue$V1==0] <- 10^-20
  pvalues[,i] <- pvalue$V1
}

## altas of selection signatures

count <- apply(pvalues,1,function(x) sum(x< 0.05/1174429/30))

selection_atlas <- as.data.frame(matrix(nrow=1174429,ncol=3))
colnames(selection_atlas) <- c("chr","pos","sig_count")
selection_atlas$chr <- locations[,1]
selection_atlas$pos <- locations[,2]
selection_atlas$sig_count <- count
write.table(selection_atlas,paste0("data-intermediate/",prefix,"selection_atlas_snp_sig_counts.txt"),append = F,quote = F,sep = "\t",row.names = F,col.names = T)

ggplot(selection_atlas,mapping=aes(x=sig_count)) +
  geom_histogram(bins = 100)+
  theme_classic()


### Permutation ### 
count_permutation <- c()
pvalues_permutated <- as.data.frame(matrix(nrow=1174429,ncol=30))

for(i in 1:100){
  pvalues_permutated <- apply(pvalues,2,function(x) sample(x))
  count_permutation <- c(count_permutation,apply(pvalues_permutated,1,function(x) sum(x< 0.05/1174429/30)))
  print(i)
}

hist(count_permutation)
quantile(count_permutation,1-0.000001)
         
sig_snp_locations <- locations[count>5,]


## draw the manhattan plot using the lrt1 parallelism result
source("scripts/manhattan_plot_Xing.R")

data <- as.data.frame(matrix(nrow=nrow(pvalues),ncol=4))
colnames(data) <- c("CHR","BP","SNP","P")
data$CHR <- locations[,1]
data$BP <- locations[,2]
data$SNP <- locations[,3]
data$P <- 10^(-selection_atlas$sig_count)
sig_snps <- data$SNP[count>5]

png(filename =paste0("figs/",prefix,"selection_atlas_snp_sig_counts.png"),width = 8,height = 3,units = "in",res = 1000)
manhattan_new(data,suggestiveline = F,genomewideline = F,highlight = sig_snps,col = c("grey40","grey80"),cex=0.6,ylim=c(0,15),
              xlab = "",xaxt = "n",ylab=name,yaxt="n")
dev.off()



#write.table(sig_snp_locations[,c(1,3,4)],"parallel_shifted_blocks_snps_4.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = F)


## haplotype-based atlas 

## read the haplotype blocks
blocks <- read.delim("data-intermediate/1001g_grenet_genomewide_partition.txt",header = F)
block_pvalue <- as.data.frame(matrix(NA,nrow=nrow(blocks),ncol=nrow(site_list)))

for(i in 1:nrow(blocks)){
  index <- which(locations$V1==blocks$V1[i] & locations$V2 >= blocks$V2[i] & locations$V2 <= blocks$V3[i] )
  if (length(index) > 0){
    block_pvalue[i,] <- apply(pvalues[index,],2,min)
    print(i)
  }
}

sum(is.na(block_pvalue$V2))

haplotype_count <- apply(block_pvalue,1,function(x) sum(x< 0.05/1174429/30))

selection_atlas_haplotype <- as.data.frame(matrix(nrow=nrow(block_pvalue),ncol=4))
colnames(selection_atlas_haplotype) <- c("chr","start","end","sig_count")
selection_atlas_haplotype$chr <- blocks[,1]
selection_atlas_haplotype$start <- blocks[,2]
selection_atlas_haplotype$end <- blocks[,3]
selection_atlas_haplotype$sig_count <- haplotype_count


write.table(selection_atlas_haplotype,paste0("data-intermediate/",prefix,"selection_atlas_haplotype_sig_counts.txt"),append = F,quote = F,sep = "\t",
                                             row.names = F,col.names = T)

### Permutation ### 
haplotype_count_permutation <- c()
block_pvalues_permutated <- as.data.frame(matrix(nrow=nrow(block_pvalue),ncol=nrow(site_list)))
for(i in 1:5000){
  #block_pvalues_permutated <- apply(haplotype_fdr,2,function(x) sample(x))
  #haplotype_count_permutation <- c(haplotype_count_permutation,apply(block_pvalues_permutated,1,function(x) sum(x<haplotype_fdr_thred)))
  block_pvalues_permutated <- apply(block_pvalue,2,function(x) sample(x))
  haplotype_count_permutation <- c(haplotype_count_permutation,apply(block_pvalues_permutated,1,function(x) sum(x<0.05/1174429/30,na.rm = T)))
}

hist(haplotype_count_permutation)
quantile(haplotype_count_permutation,1-0.0000001)
table((haplotype_count >10 & !is.na(haplotype_count)))

haplotype_count_df <- tibble(count = haplotype_count)

ggplot(haplotype_count_df,aes(x=count))+
  geom_histogram(bins = 40,binwidth = 1)+
  geom_vline(xintercept=10,colour = "red")+
  xlab("Number of gardens \n the haplotype block was under selection")+
  theme_classic(base_size = 18)

ggsave(paste0("figs/",prefix,"selection_atlas_haplotype_block_histogram.pdf"),width = 5,height = 4,)

#pdf(paste0("figs/",prefix,"selection_atlas_haplotype_block_histogram.pdf"),width = 5,height = 4,useDingbats = F)
#hist(haplotype_count,xlab="Number of gardens \n the haplotype block was under selection",main=NA)
#dev.off()

significant_blocks <- blocks[haplotype_count >8 & !is.na(haplotype_count),]
write.table(significant_blocks,paste0("data-intermediate/",prefix,"parallel_shifted_haplotype_blocks_8.txt"),
            append = F,quote = F,sep = "\t",row.names = F,col.names = F)
significant_blocks <- blocks[haplotype_count >10 & !is.na(haplotype_count),]
write.table(significant_blocks,paste0("data-intermediate/",prefix,"parallel_shifted_haplotype_blocks_10.txt"),
            append = F,quote = F,sep = "\t",row.names = F,col.names = F)
significant_blocks <- blocks[haplotype_count >11 & !is.na(haplotype_count),]
write.table(significant_blocks,paste0("data-intermediate/",prefix,"parallel_shifted_haplotype_blocks_11.txt"),
            append = F,quote = F,sep = "\t",row.names = F,col.names = F)
significant_blocks <- blocks[haplotype_count >11 & !is.na(haplotype_count),]
write.table(significant_blocks,paste0("data-intermediate/",prefix,"parallel_shifted_haplotype_blocks_11.txt"),
            append = F,quote = F,sep = "\t",row.names = F,col.names = F)
significant_blocks <- blocks[haplotype_count >12 & !is.na(haplotype_count),]
write.table(significant_blocks,paste0("data-intermediate/",prefix,"parallel_shifted_haplotype_blocks_12.txt"),
            append = F,quote = F,sep = "\t",row.names = F,col.names = F)
significant_blocks <- blocks[haplotype_count >13 & !is.na(haplotype_count),]
write.table(significant_blocks,paste0("data-intermediate/",prefix,"parallel_shifted_haplotype_blocks_13.txt"),
            append = F,quote = F,sep = "\t",row.names = F,col.names = F)
significant_blocks <- blocks[haplotype_count >14 & !is.na(haplotype_count),]
write.table(significant_blocks,paste0("data-intermediate/",prefix,"parallel_shifted_haplotype_blocks_14.txt"),
            append = F,quote = F,sep = "\t",row.names = F,col.names = F)
