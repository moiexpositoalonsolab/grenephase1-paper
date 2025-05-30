### Goal: to explore frequency trajectories in genes of possible importance
# The idea is to unpack the LFMM or others into something that is
# visually clear

################################################################################
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
library(topr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

#####*********************************************************************######
### FUNCTIONS ####


read_merged_allele_subset<-function(headfile="data-intermediate/merged_hapFIRE_allele_frequency-header.csv",
                                    subsetfile="data-intermediate/merged_hapFIRE_allele_frequency-AT5G10140-FLC.csv",
                                    snpfile="data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-AT5G10140-FLC.csv",
                                    frequencyp0file="data/average_seedmix_p0.txt",
                                    worldclimsitesdatafile="grene/data/worldclim_sitesdata.rda",
                                    sampledatafile="grene/data/samples_data.rda",
                                    addclimate=TRUE,
                                    addflowers=FALSE
){


  # Read header
  ahead<-read.csv(headfile,header = F)
  # Read SNP locations
  alocations<-read.table(snpfile,header = F)
  # Read the allele frequencies
  a<-read.csv(subsetfile,header = F)
  # Add columns needed
  colnames(a)<-ahead
  a$snp=alocations$V2
  a$chr=alocations$V1
  a$pos=alocations$V4

  # Pivot to add other things
  a<-
    a %>%
    # head %>%
    pivot_longer(
      cols = -c(snp, chr, pos),
      names_to = c("site", "year", "rep"),
      # names_pattern = "X(\\d+)_(\\d+)_(\\d+)",
      names_pattern = "(\\d+)_(\\d+)_(\\d+)",
      values_to = "freq"
    )

  # Merge with starting frequencies
  p0<-read.table(frequencyp0file, header=F)
  p0$snp<-paste0(p0$V1,"_",p0$V2)
  p0<-p0 %>% dplyr::rename(chr=V1,pos=V2,startfreq=V3)
  a<-
    p0 %>%
    dplyr::select(startfreq, snp) %>%
    merge(.,a,by='snp')
  # add the year zero
  azero<-dplyr::filter(a, year==1) %>%
    mutate(year=0, freq=startfreq)
  a<-rbind(a,azero)
  # Merge with climate
  if(addclimate){
    # Read climate
    load(worldclimsitesdatafile)
    a<-
      worldclim_sitesdata %>%
      dplyr::select(site,starts_with("bio")) %>%
      merge(.,a,by='site')
  }
  # merge with flower numbers
  if(addflowers){
    # Read sample size
    load("grene/data/samples_data.rda")
    samples_tmp<-samples_data %>%
      dplyr::filter(year<2021) %>%
      dplyr::mutate(year= year-2017) %>%
      dplyr::rename(flowers=flowerscollected) %>%
      dplyr::rename(rep=plot) %>%
      dplyr::select(site, year, rep, flowers) %>%
      group_by(site,year,rep) %>%
      summarise(flowers=sum(flowers))
    # dplyr::select(site,plot,year, date,month,day,flowers) %>%
    # dplyr::mutate(julian_day = yday(ymd(date)))
    a<-
      a %>%
      merge(.,
            samples_tmp, by.x=c("site","year","rep"), by.y=c("site","year","rep"), all.x=T
      )
  }
  # END
  return(a)
}



#####*********************************************************************######

#### Get intermedaite datasets ####

# Intermediate datasets
p0<-read.table("data/average_seedmix_p0_LDpruned.txt", header=F)
p0<-read.table("data/average_seedmix_p0.txt", header=F)
p0$snp<-paste0(p0$V1,"_",p0$V2)
p0<-p0 %>% dplyr::rename(chr=V1,pos=V2,startfreq=V3)


# Site cliamte
load("grene/data/worldclim_sitesdata.rda")

# TAIR
tair<-read.csv("data-external/TAIR10_parsedgenes.csv",header=T)
head(tair)

# Annotated SNPs
load("data-intermediate/alleles-TAIR10_parsedgenes.rda")
annotated_snps %>% dim
annotated_snps %>% head

#####*********************************************************************######
##### GET GENE SUBSETS #####
###################################
# Extract FLC
# Define your start and end rows as variables
atname="AT1G45249"
commonname="ABF2"
# Define start and end
buffer= 100
myalleles <- annotated_snps %>% dplyr::filter(closest_gene==atname)
start_row <- which(annotated_snps$snp==head(myalleles$snp,1)) - buffer
end_row <- which(annotated_snps$snp==tail(myalleles$snp,1)) + buffer

# start_row
# end_row
# start_row <- 100
# end_row <- 500

# Construct the sed command dynamically

subsetfile<-paste0("data-intermediate/merged_hapFIRE_allele_frequency-",atname,"-",commonname,".csv")
snpfile<-paste0("data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-",atname,"-",commonname,".csv")

sed_command <- paste("sed -n '", start_row, ",", end_row, "p' data/merged_hapFIRE_allele_frequency.csv > ", subsetfile,sep="")
sed_command_locations <- paste("sed -n '", start_row, ",", end_row, "p' data/greneNet_final_v1.1.bim > ",snpfile,atname,"-",commonname,".csv", sep="")
cat(sed_command)
cat(sed_command_locations)

run=FALSE
if(run) system(paste(sed_command_locations," &"), intern=TRUE)
if(run) system(paste(sed_command_locations," &"), intern=TRUE)


#####*********************************************************************######
##### Read ABF2 #####

tairraw<-read.table("data-external/TAIR10_GFF3_genes_transposons.gff")

orientation="reverse"
gene_start=17165125
gene_end=17167741
buffer=1300
chr="Chr1"

tairsub<-
  tairraw %>%
  dplyr::filter(grepl(atname,V9))

tairsub<- # just the region
  tairraw %>%
  dplyr::filter(V1==chr) %>%
  dplyr::filter(V4 > gene_start-1000) %>%
  dplyr::filter(V5 < gene_end+1300)

exon_data <-
  tairsub %>%
  filter(V3 == "exon") %>%
  mutate(start = as.numeric(V4),
         end = as.numeric(V5),
         isoform = V9)  # Extract exon information
#####*********************************************************************######
##### Read ABF2 #####
##### Read gene #####
x
s
ub<-read_merged_allele_subset(
  subsetfile = subsetfile,
  headfile="data-intermediate/merged_hapFIRE_allele_frequency-header.csv",
  snpfile=snpfile
)

head(sub)
dim(sub)

# get LMs just for corroboration
subres<-
  sub %>%
  dplyr::group_by(snp) %>%
  dplyr::summarize(r=cor.test(freq-startfreq,bio1, method='k')$estimate,
                   p=cor.test(freq-startfreq,bio1, method='k')$p.value) %>%
  merge(.,sub,by="snp")
head(subres)
dim(subres)


xx# minimanhattan
buffer=1300
minimanhattan<-
subres %>%
  ggplot(.)+
  geom_point(aes(y=-log10(p), x=pos))+
  geom_vline( xintercept = gene_start)+
  geom_vline( xintercept = gene_end)+
  xlim(c(gene_start-buffer, gene_end+buffer))+
  theme_minimal()

minimanhattan+
geom_rect(data=exon_data,  # adding the genes
          aes(xmin = start, xmax = end,
              ymin = -as.numeric(as.factor(isoform)) - 0.4, ymax = -as.numeric(as.factor(isoform)) -0.1
          ), fill = "grey")


top<-
  subres %>%
  dplyr::filter(pos >gene_start, pos<gene_end) %>%
  dplyr::filter(p==min(p)) %>%
  arrange(p) %>%
  head
top

# top snp
#1_17166232

sub %>%
  dplyr::filter(snp=="1_17166232") %>%
  ggplot(.)+
  geom_point(aes(y= freq*100, x=bio1, color=bio1))+
  stat_smooth(aes(y= freq*100, x=bio1, color=bio1),method='glm',color='grey')+
  scale_color_gradientn(name = "",colours = rev(redblue))+
  scale_y_log10()+
  ylab("Allele frequency (%)")+
  xlab("Annual Temperature (C)")+
  theme_minimal()
