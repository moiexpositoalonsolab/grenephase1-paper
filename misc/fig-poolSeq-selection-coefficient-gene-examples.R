################################################################################
### Goal
### Estimate selection coefficient of top SNPs
################################################################################

library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives/MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses" # GOOGLE DRVIE
setwd(myfolder)

# ################################################################################
# ################################################################################
# We can work on the equations form poolSeq
# ################################################################################
# ################################################################################



# ################################################################################
# # Make a sync file
# 
# # Read the allele dataset of a region in the genomes, subsetted from the 3 million
# alleledata<-read.csv("data-intermediate/merged_hapFIRE_allele_frequency-AT2G27030-CAM5.csv",
#                   header=F)
# # Read the positions
# positiondata<-read.table("data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-AT2G27030-CAM5.csv",
#                     header=F)
# # Get the columns, and parse them
# colunmanmespositions<-c("chr","posID","unknown","pos","ref","alt")
# columnnames<-read.csv("data-intermediate/merged_hapFIRE_allele_frequency-header.csv",
#                       header=F)
#   split_names <- strsplit(as.character(columnnames), split = "_",fixed = T)
#   # Convert the list to a matrix for easier indexing
#   split_matrix <- do.call(rbind, split_names)
#   # Extract the components into separate vectors
#   site <- split_matrix[, 1]
#   year <- split_matrix[, 2]
#   rep  <- split_matrix[, 3]
# # Make the column names with freq
# columnnames<-paste0(columnnames, ".freq")
# 
# # Bind all the dataset
# combineddf<-cbind(positiondata,alleledata)
# colnames(combineddf) <- c(colunmanmespositions,columnnames)
# 
# 
# write.table(file="data-intermediate/merged_hapFIRE_allele_frequency-AT2G27030-CAM5-combined.csv",
#             combineddf,
#             row.names = F, col.names = T, quote = F)
# 
# # alleledatacombined<-read.table("data-intermediate/merged_hapFIRE_allele_frequency-AT2G27030-CAM5.csv",header = T)
# # write.csv(file="data-intermediate/merged_hapFIRE_allele_frequency-AT2G27030-CAM5.csv",
# #           alleledatacombined[,-c(1:6)],
# #             row.names = F, col.names = F, quote = F, sep=',')
# 
# ################################################################################
# library(poolSeq)
# 
# myfile="/Path/to/data.sync"
# myfile="data-intermediate/merged_hapFIRE_allele_frequency-AT2G27030-CAM5.csv"
# 
# mySync <- 
#   read.sync(file="data-intermediate/merged_hapFIRE_allele_frequency-AT2G27030-CAM5.csv", 
#             gen=year,
#             repl=rep)
# 
# 
# load("~/poolSeq-master/data/dmelER.RData")
# totest<-dmelER@alleles %>% data.frame %>%  
#   dplyr::select(-ends_with(".cov")) %>% 
#   dplyr::select(-rising)
# write.table(file="data-intermediate/example-poolSeq.sync",
#             totest,sep = "\t",
#             row.names = F, col.names = T, quote = F)
# 
# mySync <- 
#   read.sync(file="data-intermediate/example-poolSeq.sync", polarization = "minor",
#             gen=c(0, 15, 37, 59, 0, 15, 37, 59,0, 15, 37, 59),
#             repl=c(1, 1, 1,1, 2, 2, 2,2 ,3,3,3,3))
# 
# 
# mySync <- 
#   read.sync(file=myfile, 
#             gen=c(0, 10, 20, 0, 10, 20), repl=c(1, 1, 1, 2, 2, 2), rising = FALSE)
# 


################################################################################

# test Ne estimation
sub<-read_merged_allele_subset(
  subsetfile = "data-intermediate/merged_hapFIRE_allele_frequency-AT2G22540-SVP.csv",
  headfile="data-intermediate/merged_hapFIRE_allele_frequency-header.csv",
  snpfile="data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-AT2G22540-SVP.csv",
  addclimate = T, addflowers = T
)
head(sub)
sub$flower %>% unique
sub$year %>% unique
# View(head(sub))

subsub<-sub %>% 
  dplyr::filter(site==4, year %in% c(1,3)) %>% 
  dplyr::filter(rep==1)

################################################################################
# Get the data
p0=subsub$startfreq
pt=subsub$freq
cov0=rep(1000,nrow(subsub))
covt=rep(1000,nrow(subsub))
Ncensus=5000
poolSize=subsub$flowers %>% unique
t=2
# myne<-
  estimateNe(p0 = p0, 
             pt = pt,
             t=1,
             # method="JR.planI",
             method=c("P.planI", "P.planII",
               "JR.planI", "JR.planII",
               "W.planI", "W.planII",
               "P.alt.1step.planII", "P.alt.2step.planI", "P.alt.2step.planII"),
             cov0 = cov0,
             covt = covt, 
             Ncensus = Ncensus,
             poolSize = poolSize
             ) 

sub$snp %>% unique %>% length

wrapselection<-function(dat=sub,site=4){
  # Underlying function
  wrapselection_<-function(dat=sub, site=4, mysnp="2_9559908"){
    # Step 1: Select only the relevant columns
    sub_selected <- sub %>% 
      dplyr::filter(site==4, snp==mysnp) %>% 
      dplyr::select(rep, year, freq) %>% 
      dplyr::arrange(year)
    # Step 2: Pivot the data to wider format
    sub_pivoted <- sub_selected %>%
      pivot_wider(
        names_from = year,    # Columns will be timepoints (years)
        values_from = freq    # Values will be the 'freq' column
      ) %>% 
      dplyr::select(-rep)
    # View the pivoted dataset
    # print(sub_pivoted)
    estimateSH(traj = sub_pivoted,t = c(1:4), Ne=1,haploid = TRUE,h = 0.5)$s
  }
  # Go through the dataset
  res<-lapply(unique(sub$snp), function(i) wrapselection_(sub,site,i) ) %>% unlist
  result<-data.frame(snp=unique(sub$snp), site=site, s=res)
}
# Calculate seleciton
selectionres<-wrapselection(sub,4)
  # selectionresbackup<-selectionres
# Explore
plot(selectionres$s)
hist(selectionres$s)
# Change column
selectionres<-
selectionres %>% 
  data.frame %>% 
  dplyr::mutate(snps=snp) %>% 
  separate(snps, into = c("chr", "pos"), sep='_') %>% 
  dplyr::mutate(chr=as.numeric(chr),pos=as.numeric(pos)) %>% 
  merge(dplyr::select(p0,startfreq,snp), by='snp') %>% 
  dplyr::mutate(p0fold= ifelse(startfreq>=0.5,1-startfreq, startfreq)) %>% 
  dplyr::rename(p0=startfreq)

ggplot(selectionres)+
  geom_point(aes(y=s, x=pos))+
  theme_minimal()

ggplot(selectionres)+
  geom_point(aes(y=(s), x=p0))+
  theme_minimal()

ggplot(selectionres)+
  geom_point(aes(y=abs(s), x=p0))+
  theme_minimal()
ggplot(selectionres)+
  geom_point(aes(y=abs(s), x=p0fold))+
  scale_x_log10()+
  theme_minimal()
