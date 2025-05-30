rm(list=ls())

library(ggplot2)
library(boot)

## POPULATION EVOLUTION

## This script will analyze the parallelism in GrENE-net using per-generation frequency change
## both ecotype and SNP frequency change

prefix <- "POP_EVOLUTION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)


meta <- read.delim("data/merged_sample_table.csv",sep=",",stringsAsFactors = F)
terminal_index <- read.delim("data/site_plot_terminal_index.txt")

## SNP frequencies
ld_p0 <- read.delim("data/average_seedmix_p0_LDpruned.txt",header=F)
ld_p <- read.delim("data/merged_hapFIRE_allele_frequency_LDpruned.txt",header=T,check.names = F)

## ecotype frequencies
ecotype_p0 <- read.delim("data/founder_ecotype_frequency.txt",header=F)
ecotype_p <- read.delim("data/merged_ecotype_frequency.txt",header=T,check.names = F)


##### calculate the per generation delta frequency change: compare to the previous generation

## generation 1
index_gen_1 <- which(terminal_index$generation1==1)
## generation 2
index_gen_2 <-  which(terminal_index$generation1==1 &terminal_index$generation2==2)
## generation 3
index_gen_3 <-  which(terminal_index$generation1==1 &terminal_index$generation2==2 & terminal_index$generation3==3)


site_plot_index <- c(terminal_index$index1[terminal_index$generation1==1],
                     terminal_index$index2[terminal_index$generation1==1 & terminal_index$generation2==2],
                     terminal_index$index3[terminal_index$generation1==1 & terminal_index$generation2==2 & terminal_index$generation3==3])
experimental_sites <- unique(meta$site[site_plot_index])

experimental_meta <- meta[site_plot_index,]


## define the mean pearson dataframe 
snp_parallelism <- as.data.frame(matrix(nrow=31*3,ncol=5))
colnames(snp_parallelism) <- c("site","mean","lower","upper","generation")
snp_parallelism$site <- c(rep(unique(experimental_meta$site),3))
snp_parallelism$generation <- factor(c(rep(1,31),rep(2,31),rep(3,31)),levels = c(1,2,3))


ecotype_parallelism <- as.data.frame(matrix(nrow=31*3,ncol=5))
colnames(ecotype_parallelism) <- c("site","mean","lower","upper","generation")
ecotype_parallelism$site <- c(rep(unique(experimental_meta$site),3))
ecotype_parallelism$generation <- factor(c(rep(1,31),rep(2,31),rep(3,31)),levels = c(1,2,3))


## define the samplemean function for bootstrap 

samplemean <- function(data,i){return(mean(data[i]))}


## for generation 0 to 1

site_gen1 <- meta$site[terminal_index$index1[index_gen_1]]
unique_site_gen1 <- unique(meta$site[terminal_index$index1[index_gen_1]])
delta_snp_gen1 <- apply(ld_p[,terminal_index$index1[index_gen_1]],2,function(x) x-ld_p0$V3)
delta_ecotype_gen1 <- apply(ecotype_p[,terminal_index$index1[index_gen_1]],2,function(x) x-ecotype_p0$V2)

for(i in 1:length(unique_site_gen1)){
  site_ <- unique_site_gen1[i]
  if (sum(site_gen1 == site_) > 1) {
    #snp correlation
    delta_snp_site_ <- delta_snp_gen1[, which(site_gen1 == site_)]
    snp_mat <- cor(delta_snp_site_,method = "pearson")
    snp_r <- snp_mat[upper.tri(snp_mat)]
    snp_b <- boot(snp_r,samplemean,R=10000)
    snp_ci <- boot.ci(snp_b,type = "basic")
    snp_parallelism[snp_parallelism$site==site_ & snp_parallelism$generation==1,2] <- mean(snp_b$t)
    snp_parallelism[snp_parallelism$site==site_ & snp_parallelism$generation==1,3:4] <- snp_ci$basic[4:5]
    
    #ecotype correlation
    delta_ecotype_site_ <- delta_ecotype_gen1[, which(site_gen1 == site_)]
    ecotype_mat <- cor(delta_ecotype_site_,method = "pearson")
    ecotype_r <- ecotype_mat[upper.tri(ecotype_mat)]
    ecotype_b <- boot(ecotype_r,samplemean,R=10000)
    ecotype_ci <- boot.ci(ecotype_b,type = "basic")
    ecotype_parallelism[ecotype_parallelism$site==site_ & ecotype_parallelism$generation==1,2] <- mean(ecotype_b$t)
    ecotype_parallelism[ecotype_parallelism$site==site_ & ecotype_parallelism$generation==1,3:4] <- ecotype_ci$basic[4:5]
  }
}


## for generation 1 to 2

site_gen2 <- meta$site[terminal_index$index2[index_gen_2]]
unique_site_gen2 <- unique(meta$site[terminal_index$index2[index_gen_2]])
delta_snp_gen2 <- ld_p[,terminal_index$index2[index_gen_2]] - ld_p[,terminal_index$index1[index_gen_2]]
delta_ecotype_gen2 <- ecotype_p[,terminal_index$index2[index_gen_2]] - ecotype_p[,terminal_index$index1[index_gen_2]]

for(i in 1:length(unique_site_gen2)){
  site_ <- unique_site_gen2[i]
  if (sum(site_gen2 == site_) > 1) {
    #snp correlation
    delta_snp_site_ <- delta_snp_gen2[, which(site_gen2 == site_)]
    snp_mat <- cor(delta_snp_site_,method = "pearson")
    snp_r <- snp_mat[upper.tri(snp_mat)]
    snp_b <- boot(snp_r,samplemean,R=10000)
    snp_ci <- boot.ci(snp_b,type = "basic")
    snp_parallelism[snp_parallelism$site==site_ & snp_parallelism$generation==2,2] <- mean(snp_b$t)
    if (!is.null(snp_ci)){
      snp_parallelism[snp_parallelism$site==site_ & snp_parallelism$generation==2,3:4] <- snp_ci$basic[4:5]
    }
    
    #ecotype correlation
    delta_ecotype_site_ <- delta_ecotype_gen2[, which(site_gen2 == site_)]
    ecotype_mat <- cor(delta_ecotype_site_,method = "pearson")
    ecotype_r <- ecotype_mat[upper.tri(ecotype_mat)]
    ecotype_b <- boot(ecotype_r,samplemean,R=10000)
    ecotype_ci <- boot.ci(ecotype_b,type = "basic")
    ecotype_parallelism[ecotype_parallelism$site==site_ & ecotype_parallelism$generation==2,2] <- mean(ecotype_b$t)
    if (!is.null(ecotype_ci)){
      ecotype_parallelism[ecotype_parallelism$site==site_ & ecotype_parallelism$generation==2,3:4] <- ecotype_ci$basic[4:5]
    }
  }
}

## for generation 2 to 3

site_gen3 <- meta$site[terminal_index$index3[index_gen_3]]
unique_site_gen3 <- unique(meta$site[terminal_index$index3[index_gen_3]])
delta_snp_gen3 <- ld_p[,terminal_index$index3[index_gen_3]] - ld_p[,terminal_index$index2[index_gen_3]]
delta_ecotype_gen3 <- ecotype_p[,terminal_index$index3[index_gen_3]] - ecotype_p[,terminal_index$index2[index_gen_3]]

for(i in 1:length(unique_site_gen3)){
  site_ <- unique_site_gen3[i]
  if (sum(site_gen3 == site_) > 1) {
    #snp correlation
    delta_snp_site_ <- delta_snp_gen3[, which(site_gen3 == site_)]
    snp_mat <- cor(delta_snp_site_,method = "pearson")
    snp_r <- snp_mat[upper.tri(snp_mat)]
    snp_b <- boot(snp_r,samplemean,R=10000)
    snp_ci <- boot.ci(snp_b,type = "basic")
    snp_parallelism[snp_parallelism$site==site_ & snp_parallelism$generation==3,2] <- mean(snp_b$t)
    if (!is.null(snp_ci)){
      snp_parallelism[snp_parallelism$site==site_ & snp_parallelism$generation==3,3:4] <- snp_ci$basic[4:5]
    }
    
    #ecotype correlation
    delta_ecotype_site_ <- delta_ecotype_gen3[, which(site_gen3 == site_)]
    ecotype_mat <- cor(delta_ecotype_site_,method = "pearson")
    ecotype_r <- ecotype_mat[upper.tri(ecotype_mat)]
    ecotype_b <- boot(ecotype_r,samplemean,R=10000)
    ecotype_ci <- boot.ci(ecotype_b,type = "basic")
    ecotype_parallelism[ecotype_parallelism$site==site_ & ecotype_parallelism$generation==3,2] <- mean(ecotype_b$t)
    if (!is.null(ecotype_ci)){
      ecotype_parallelism[ecotype_parallelism$site==site_ & ecotype_parallelism$generation==3,3:4] <- ecotype_ci$basic[4:5]
    }
  }
}


combined <- rbind(snp_parallelism,ecotype_parallelism)
combined$source <- c(rep("snp",31*3),rep("ecotype",31*3))


combined_summary <- as.data.frame(matrix(ncol=5,nrow=2*3))
colnames(combined_summary) <-c("generation","mean_r","upper","lower","source")
combined_summary$generation <- rep(c(1,2,3),2)
combined_summary$source <- c(rep("snp",3),rep("ecotype",3))

for(i in 1:3){
  combined_summary[i,2] <- mean(snp_parallelism$mean[snp_parallelism$generation==i],na.rm=T)
  combined_summary[i,3] <- mean(snp_parallelism$mean[snp_parallelism$generation==i],na.rm=T) +  sd(snp_parallelism$mean[snp_parallelism$generation==i],na.rm=T)
  combined_summary[i,4] <- mean(snp_parallelism$mean[snp_parallelism$generation==i],na.rm=T) -  sd(snp_parallelism$mean[snp_parallelism$generation==i],na.rm=T)
  
  combined_summary[i+3,2] <- mean(ecotype_parallelism$mean[ecotype_parallelism$generation==i],na.rm=T)
  combined_summary[i+3,3] <- mean(ecotype_parallelism$mean[ecotype_parallelism$generation==i],na.rm=T) +  sd(ecotype_parallelism$mean[ecotype_parallelism$generation==i],na.rm=T)
  combined_summary[i+3,4] <- mean(ecotype_parallelism$mean[ecotype_parallelism$generation==i],na.rm=T) -  sd(ecotype_parallelism$mean[ecotype_parallelism$generation==i],na.rm=T)
}


ggplot(combined_summary,aes(x=generation,y=mean_r,colour = source)) +
  geom_point(position = position_dodge(0.4),cex=5,pch=16)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1,position=position_dodge(.4))+
  ylim(c(-0.08,0.6))+
  ylab("Parallelism (r)")+
  scale_color_manual(values =c("#004080","#e50000"),labels=c("founder","LD-pruned SNPs"))+
  scale_x_discrete(name ="Generation", limits=c("1","2","3"))+
  theme_classic(base_size = 18)
ggsave(filename =  paste0("figs/",prefix,"per_generation_parallelism_summary.pdf"),device = "pdf",width = 6,height = 4,units = "in")

  
combined_gen1 <- combined %>%
  filter(!is.na(mean)) %>%
  filter(generation==1)
combined_gen1$site
order_site <- combined_gen1$site[order(combined_gen1$mean[combined_gen1$source=="snp"],decreasing = F)]
combined_gen1$site <- factor(combined_gen1$site,levels = order_site)
write.table(combined_gen1,file = "data-intermediate/generation_1_parallelism.txt",append = F,quote = F,sep = "\t",row.names = F)

p <- ggplot(combined_gen1,aes(x=site,y=mean,color=source))+
  coord_flip()+
  geom_point(position = position_dodge(.3),cex=2.5,pch=16)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0,position=position_dodge(.3),linewidth = 0.5)+
  geom_hline(yintercept = 0.294,linetype="dotted")+
  scale_color_manual(values =c("#004080","#e50000"))+
  ylab("Mean parallelism (r)")+
  theme_classic()
p
ggsave(plot = p,filename =  paste0("figs/",prefix,"generation_1_parallelism.pdf"),device = "pdf",width = 4,height = 12,units = "in")
combined_gen1$r <- combined_gen1$mean>0.3
combined_gen1 %>%
  filter(source=="snp") %>%
  ggplot(.,aes(x=site,y=mean,colour = r))+
  geom_point(position = position_dodge(.3),cex=7,pch=16)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0,position=position_dodge(.3),linewidth = 0.5)+
  scale_color_manual(values =c("grey70","#e50000"))+
  ylab("Mean parallelism (r)")+
  theme_classic(base_size = 18)+
  theme(legend.position="none")
## order by Bio 1
site_climate <- read.delim("data-external/bioclimvars_sites_era5_year_2018.csv",sep=",") %>%
  filter(site %in% combined_gen1$site) %>%
  arrange(match(site,combined_gen1$site))

order_site_bio1 <- site_climate$site[order(site_climate$bio1,decreasing = F)]
combined_gen1$site <- factor(combined_gen1$site,levels = order_site_bio1)

p <- ggplot(combined_gen1,aes(x=site,y=mean,color=source))+
  geom_point(position = position_dodge(.3),cex=5,pch=3)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1,position=position_dodge(.3),linewidth = 0.7)+
  geom_hline(yintercept = 0,linetype="dotted")+
  geom_hline(yintercept = mean(combined_gen1$mean[combined_gen1$source=="snp"]),linetype="dotted")+
  scale_color_manual(values =c("#004080","#e50000"))+
  ylab("Mean parallelism (r)")+
  theme_classic()
p
ggsave(plot = p,filename =  paste0("figs/",prefix,"generation_1_parallelism_Bio1_ordered_errorbar.pdf"),device = "pdf",width = 12,height = 4,units = "in")

mean(combined_gen1$mean[combined_gen1$source=="snp"])
mean(combined_gen1$mean[combined_gen1$source=="ecotype"])

combined_snp <- combined %>%
  filter(!is.na(mean)) %>%
  filter(source=="snp")
combined_snp$site <- factor(combined_snp$site,levels = order_site)

p_snp <- ggplot(combined_snp,aes(x=site,y=mean,colour = generation,shape = generation))+
  coord_flip()+
  geom_point(position = position_dodge(.3),cex=3,pch=16)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1,position=position_dodge(.3),linewidth = 0.7)+
  geom_hline(yintercept = 0,linetype="dotted")+
  scale_color_manual(values =c("grey20","grey50","grey80"))+
  ylab("parallelism (r)")+
  ggtitle("snp frequency change")+
  theme_classic()
ggsave(plot = p_snp,filename =  paste0("figs/",prefix,"all_generations_parallelism_snp.pdf"),device = "pdf",width = 8,height = 10,units = "in")


combined_ecotype <- combined %>%
  filter(!is.na(mean)) %>%
  filter(source=="ecotype")
combined_ecotype$site <- factor(combined_ecotype$site,levels = order_site)

p_ecotype <- ggplot(combined_ecotype,aes(x=site,y=mean,colour = generation,shape = generation))+
  coord_flip()+
  geom_point(position = position_dodge(.3),cex=3,pch=16)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1,position=position_dodge(.3),linewidth = 0.7)+
  geom_hline(yintercept = 0,linetype="dotted")+
  scale_color_manual(values =c("grey20","grey50","grey80"))+
  ylab("parallelism (r)")+
  ggtitle("founder frequency change")+
  theme_classic()
ggsave(plot = p_ecotype,filename =  paste0("figs/",prefix,"all_generations_parallelism_ecotype.pdf"),device = "pdf",width = 8,height = 10,units = "in")



p_combined <- p_snp | p_ecotype

ggsave(plot = p_combined,filename =  paste0("figs/",prefix,"all_generations_parallelism_snp_ecotype.pdf"),device = "pdf",width = 8,height = 10,units = "in")




## parallelism for the terminal generation
terminal_generation <- terminal_index$index3
terminal_meta <- meta[terminal_generation,]
terminal_sites <- unique(terminal_meta$site)

delta_snp_terminal <- apply(ld_p[,terminal_generation],2,function(x) x-ld_p0$V3)
delta_ecotype_terminal <- apply(ecotype_p[,terminal_generation],2,function(x) x-ecotype_p0$V2)


terminal_parallelism <- as.data.frame(matrix(nrow=length(terminal_sites)*2,ncol=5))
colnames(terminal_parallelism) <- c("site","mean","upper","lower","source")
terminal_parallelism$site <- rep(terminal_sites,2)
terminal_parallelism$source <- c(rep("snp",31),rep("ecotype",31))

## parallelism for the terminal generation
for(i in 1:length(terminal_sites)){
  site_ <- terminal_sites[i]
  if (sum(terminal_meta$site == site_) > 1) {
    #snp correlation
    delta_snp_site_ <- delta_snp_terminal[, which(terminal_meta$site == site_)]
    snp_mat <- cor(delta_snp_site_,method = "pearson")
    snp_r <- snp_mat[upper.tri(snp_mat)]
    snp_b <- boot(snp_r,samplemean,R=10000)
    snp_ci <- boot.ci(snp_b,type = "basic")
    terminal_parallelism[terminal_parallelism$site==site_ & terminal_parallelism$source=="snp",2] <- mean(snp_b$t)
    if (!is.null(snp_ci)){
      terminal_parallelism[terminal_parallelism$site==site_ & terminal_parallelism$source=="snp",3:4] <- snp_ci$basic[4:5]
    }
    
    #ecotype correlation
    delta_ecotype_site_ <- delta_ecotype_terminal[, which(terminal_meta$site == site_)]
    ecotype_mat <- cor(delta_ecotype_site_,method = "pearson")
    ecotype_r <- ecotype_mat[upper.tri(ecotype_mat)]
    ecotype_b <- boot(ecotype_r,samplemean,R=10000)
    ecotype_ci <- boot.ci(ecotype_b,type = "basic")
    terminal_parallelism[terminal_parallelism$site==site_ & terminal_parallelism$source=="ecotype",2] <- mean(ecotype_b$t)
    if (!is.null(ecotype_ci)){
      terminal_parallelism[terminal_parallelism$site==site_ & terminal_parallelism$source=="ecotype",3:4] <- ecotype_ci$basic[4:5]
    }
  }
}

## get rid of site 33
terminal_parallelism <- terminal_parallelism %>%
  filter(site != 33)

write.table(terminal_parallelism,file=paste0("data-intermediate/",prefix,"terminal_generation_parallelism.txt"),append = F,quote = F,sep = "\t",row.names = F)

