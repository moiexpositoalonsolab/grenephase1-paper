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
atname="AT5G10140"
commonname="FLC"
# Define start and end
myalleles <- annotated_snps %>% dplyr::filter(closest_gene==atname)
  start_row <- which(annotated_snps$snp==head(myalleles$snp,1)) - 500
  end_row <- which(annotated_snps$snp==tail(myalleles$snp,1)) + 500

# start_row
# end_row
# start_row <- 100
# end_row <- 500

# Construct the sed command dynamically
sed_command <- paste("sed -n '", start_row, ",", end_row, "p' data/merged_hapFIRE_allele_frequency.csv > data-intermediate/merged_hapFIRE_allele_frequency-",atname,"-",commonname,".csv", sep="")
sed_command_locations <- paste("sed -n '", start_row, ",", end_row, "p' data/greneNet_final_v1.1.bim > data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-",atname,"-",commonname,".csv", sep="")
cat(sed_command)

# Use system() to execute the command
# system(sed_command, intern=TRUE)

###################################
# Extract DOG1
# Define your start and end rows as variables
atname="AT5G45830"
commonname="DOG1"
# Define start and end
myalleles <- annotated_snps %>% dplyr::filter(closest_gene==atname)
start_row <- which(annotated_snps$snp==head(myalleles$snp,1)) - 500
end_row <- which(annotated_snps$snp==tail(myalleles$snp,1)) + 500
# start_row <- 100
# end_row <- 500


# Construct the sed command dynamically
sed_command <- paste("sed -n '", start_row, ",", end_row, "p' data/merged_hapFIRE_allele_frequency.csv > data-intermediate/merged_hapFIRE_allele_frequency-",atname,"-",commonname,".csv", sep="")
sed_command_locations <- paste("sed -n '", start_row, ",", end_row, "p' data/greneNet_final_v1.1.bim > data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-",atname,"-",commonname,".csv", sep="")
cat(sed_command)
cat(sed_command_locations)

# Use system() to execute the command
# system(sed_command, intern=TRUE)


###################################
# Extract CAM5
# Define your start and end rows as variables
atname="AT2G27030"
commonname="CAM5"
# Define start and end
myalleles <- annotated_snps %>% dplyr::filter(closest_gene==atname)
start_row <- which(annotated_snps$snp==head(myalleles$snp,1)) - 500
end_row <- which(annotated_snps$snp==tail(myalleles$snp,1)) + 500
# start_row <- 100
# end_row <- 500


# Construct the sed command dynamically
sed_command <- paste("sed -n '", start_row, ",", end_row, "p' data/merged_hapFIRE_allele_frequency.csv > data-intermediate/merged_hapFIRE_allele_frequency-",atname,"-",commonname,".csv", sep="")
sed_command_locations <- paste("sed -n '", start_row, ",", end_row, "p' data/greneNet_final_v1.1.bim > data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-",atname,"-",commonname,".csv", sep="")
cat(sed_command)
cat(sed_command_locations)

# Use system() to execute the command
# system(sed_command, intern=TRUE)

###################################
# Extract SVP
atname="AT2G22540"
commonname="SVP"
	
# Define start and end
myalleles <- annotated_snps %>% dplyr::filter(closest_gene==atname)
start_row <- which(annotated_snps$snp==head(myalleles$snp,1)) - 500
end_row <- which(annotated_snps$snp==tail(myalleles$snp,1)) + 500

# Construct the sed command dynamically
sed_command <- paste("sed -n '", start_row, ",", end_row, "p' data/merged_hapFIRE_allele_frequency.csv > data-intermediate/merged_hapFIRE_allele_frequency-",atname,"-",commonname,".csv", sep="")
sed_command_locations <- paste("sed -n '", start_row, ",", end_row, "p' data/greneNet_final_v1.1.bim > data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-",atname,"-",commonname,".csv", sep="")
cat(sed_command)
cat(sed_command_locations)

###################################
###### Extract SVP #####
atname="AT4G20370"
commonname="TSF"

# Define start and end
myalleles <- annotated_snps %>% dplyr::filter(closest_gene==atname)
start_row <- which(annotated_snps$snp==head(myalleles$snp,1)) - 500
end_row <- which(annotated_snps$snp==tail(myalleles$snp,1)) + 500

# Construct the sed command dynamically
sed_command <- paste("sed -n '", start_row, ",", end_row, "p' data/merged_hapFIRE_allele_frequency.csv > data-intermediate/merged_hapFIRE_allele_frequency-",atname,"-",commonname,".csv", sep="")
sed_command_locations <- paste("sed -n '", start_row, ",", end_row, "p' data/greneNet_final_v1.1.bim > data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-",atname,"-",commonname,".csv", sep="")
cat(sed_command)
cat(sed_command_locations)

#### Extract CAM5 #####

# Extract CAM5
# Define your start and end rows as variables
atname="AT1G52180"
commonname="CAM5"
# Define start and end
myalleles <- annotated_snps %>% dplyr::filter(closest_gene==atname)
start_row <- which(annotated_snps$snp==head(myalleles$snp,1)) - 500
end_row <- which(annotated_snps$snp==tail(myalleles$snp,1)) + 500

sed_command <- paste("sed -n '", start_row, ",", end_row, "p' data/merged_hapFIRE_allele_frequency.csv > data-intermediate/merged_hapFIRE_allele_frequency-",atname,"-",commonname,".csv", sep="")
sed_command_locations <- paste("sed -n '", start_row, ",", end_row, "p' data/greneNet_final_v1.1.bim > data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-",atname,"-",commonname,".csv", sep="")
cat(sed_command)
cat(sed_command_locations)


#### Extract Aqua #####
# Define your start and end rows as variables
atname="AT1G52180"
commonname="Aquaporinlike"
# Define start and end
myalleles <- annotated_snps %>% dplyr::filter(closest_gene==atname)
start_row <- which(annotated_snps$snp==head(myalleles$snp,1)) - 500
end_row <- which(annotated_snps$snp==tail(myalleles$snp,1)) + 500

sed_command <- paste("sed -n '", start_row, ",", end_row, "p' data/merged_hapFIRE_allele_frequency.csv > data-intermediate/merged_hapFIRE_allele_frequency-",atname,"-",commonname,".csv", sep="")
sed_command_locations <- paste("sed -n '", start_row, ",", end_row, "p' data/greneNet_final_v1.1.bim > data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-",atname,"-",commonname,".csv", sep="")
cat(sed_command)
cat(sed_command_locations)


#### Extract Dormancy guy #####
# Define your start and end rows as variables
atname="AT4G19230"
commonname="CYP707A1"
# Define start and end
myalleles <- annotated_snps %>% dplyr::filter(closest_gene==atname)
start_row <- which(annotated_snps$snp==head(myalleles$snp,1)) - 500
end_row <- which(annotated_snps$snp==tail(myalleles$snp,1)) + 500

sed_command <- paste("sed -n '", start_row, ",", end_row, "p' data/merged_hapFIRE_allele_frequency.csv > data-intermediate/merged_hapFIRE_allele_frequency-",atname,"-",commonname,".csv", sep="")
sed_command_locations <- paste("sed -n '", start_row, ",", end_row, "p' data/greneNet_final_v1.1.bim > data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-",atname,"-",commonname,".csv", sep="")
cat(sed_command)
cat(sed_command_locations)

#####*********************************************************************######
#####*********************************************************************######
##### Read flc #####
# flc<-read.table

sub<-read_merged_allele_subset(
  subsetfile = "data-intermediate/merged_hapFIRE_allele_frequency-AT5G10140-FLC.csv",
  headfile="data-intermediate/merged_hapFIRE_allele_frequency-header.csv",
  snpfile="data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-AT5G10140-FLC.csv"
)

head(sub)
dim(sub)

# get LMs just for corroboration
subres<-
  sub %>% 
  dplyr::group_by(snp) %>% 
  dplyr::summarize(r=cor.test(freq-startfreq,bio1)$estimate,
                   p=cor.test(freq-startfreq,bio1)$p.value) %>% 
  merge(.,sub,by="snp")
subres %>% head


# minimanhattan  
flc_start=3172562
flc_end=3179429
subres %>% 
  ggplot(.)+
  geom_point(aes(y=-log10(p), x=pos))+
  geom_vline( xintercept = 3172562)+
  geom_vline( xintercept = 3179429)+
  # xlim(c(3172562-500,3179429+500))+
  theme_minimal()


#####*********************************************************************######
##### Read dog #####

sub<-read_merged_allele_subset(
  subsetfile = "data-intermediate/merged_hapFIRE_allele_frequency-AT5G45830-DOG1.csv",
  headfile="data-intermediate/merged_hapFIRE_allele_frequency-header.csv",
  snpfile="data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-AT5G45830-DOG1.csv"
)
head(sub)
dim(sub)

# get LMs just for corroboration
subres<-
  sub %>% 
  dplyr::group_by(snp) %>% 
  dplyr::summarize(r=cor.test(freq-startfreq,bio1)$estimate,
                   p=cor.test(freq-startfreq,bio1)$p.value) %>% 
  merge(.,sub,by="snp")
subres %>% head
# minimanhattan  
genstart=18588545
genend=18591687
subres %>% 
  ggplot(.)+
  geom_point(aes(y=-log10(p), x=pos))+
  geom_vline( xintercept = genstart)+
  geom_vline( xintercept = genend)+
  # xlim(c(18588545-500,genend+500))+
  theme_minimal()


#####*********************************************************************######
##### Read cam5 #####
cam<-read_merged_allele_subset(
          subsetfile = "data-intermediate/merged_hapFIRE_allele_frequency-AT2G27030-CAM5.csv",
          headfile="data-intermediate/merged_hapFIRE_allele_frequency-header.csv",
          snpfile="data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-AT2G27030-CAM5.csv"
          )

head(cam)
dim(cam)

# get LMs just for corroboration
camres<-
cam %>% 
  dplyr::group_by(snp) %>% 
  dplyr::summarize(r=cor.test(freq-startfreq,bio1)$estimate,
                   p=cor.test(freq-startfreq,bio1)$p.value) %>% 
  merge(.,cam,by="snp")
# minimanhattan  
camres %>% 
  ggplot(.)+
  geom_point(aes(y=-log10(p), x=pos))+
  geom_vline( xintercept = 11531967)+
  geom_vline( xintercept = 11534358)+
  xlim(c(11531967-500,11534358+500))+
  theme_minimal()
# create top snps and neutral
camrestop<-
camres %>% 
  dplyr::filter(-log10(p)>14) %>% 
  dplyr::filter(startfreq<0.5) 
camresmax<-
  camres %>% 
  dplyr::filter(-log10(p)  == max(-log10(p)))
camresneutral<-
  camres %>% 
  dplyr::mutate(roundfreq=round(freq*1000)) %>% 
  dplyr::filter(roundfreq %in% round(camrestop$freq*1000)) %>% 
  dplyr::filter(!(snp%in%camrestop$snp), -log10(p)<5) %>% 
  dplyr::filter(snp %in% sample(unique(snp), 10))

  
#plots
topsnpplot_cam5<-
camresmax %>% 
  dplyr::filter(year==1) %>% 
ggplot(.)+
  geom_point(aes((freq-startfreq),x=bio1),alpha=0.5,color="#75C376")+
  stat_summary(aes((freq-startfreq),x=bio1), size=1, alpha=1,color="#75C376")+
  stat_smooth(aes((freq-startfreq),x=bio1), method='glm', color="#75C376")+
  geom_hline(yintercept = 0,lty="dotted")+
  # scale_y_log10()+
  ylab("p1-p0")+
  xlab("Annual temperature (C)")+
  ylim(c(-0.06,+0.14))+
  theme_minimal()
topsnpplot_cam5
save_plot(plot=topsnpplot_cam5, file="figs/fig-allele_climate_lfmm_adaptive-AT2G27030-CAM5.pdf",base_height = 3.5,base_width = 3.5)
save_plot(plot=topsnpplot_cam5, file="figs/fig-allele_climate_lfmm_adaptive-AT2G27030-CAM5.png",base_height = 3.5,base_width = 3.5)
# #plots ofneutral SNps
neutralsnpplot_cam5<-
  camresneutral %>% 
  dplyr::filter(year==1) %>% 
ggplot(.)+
  geom_point(aes((freq-startfreq),x=bio1),alpha=0.5,color="grey")+
  stat_summary(aes((freq-startfreq),x=bio1), size=1, alpha=1,color="grey")+
  stat_smooth(aes((freq-startfreq),x=bio1), method='glm', color="grey")+
  geom_hline(yintercept = 0,lty="dotted")+
  ylab("p1-p0")+
  xlab("Annual temperature (C)")+
  ylim(c(-0.06,+0.14))+
  theme_minimal()
neutralsnpplot_cam5
save_plot(plot=neutralsnpplot_cam5, file="figs/fig-allele_climate_lfmm_adaptive-AT2G27030-CAM5-neutral.pdf",base_height = 3.5,base_width = 3.5)
save_plot(plot=neutralsnpplot_cam5, file="figs/fig-allele_climate_lfmm_adaptive-AT2G27030-CAM5-neutral.png",base_height = 3.5,base_width = 3.5)

# Do the temporal analysis

cam5snpstemporal<-
camresmax %>% 
  dplyr::filter(bio1>15 | bio1<8 ) %>% # only look at extremes
  dplyr::mutate(warmvscold= bio1>12) %>% 
  # dplyr::filter(year>0) %>% 
  ggplot(.)+
  geom_jitter(aes(y=freq-startfreq, x=year), alpha=0.2,width = 0.025)+
  stat_summary(aes(y=freq-startfreq, x=year, group=site))+
  stat_smooth(aes(y=freq-startfreq, x=year, group=warmvscold, color=warmvscold), method='glm')+
  scale_color_manual(values = c( "#4E99C5", "#C4403C"))+
  geom_hline(lty='dotted',yintercept = 0)+
  # geom_line(aes(y=freq, x=year, group=site))+
  facet_wrap(~warmvscold, ncol=1)+  
  theme_minimal()+
  ylab("p1-p0")+
  xlab("Year")+
  coord_cartesian(y=c(-0.15,+0.25))+
  theme(legend.position = "none")
cam5snpstemporal

save_plot(plot=cam5snpstemporal, file="figs/fig-allele_climate_lfmm_adaptive-AT2G27030-CAM5-warm-vs-cold-trajectories.pdf",base_height = 7,base_width = 3.5)
save_plot(plot=cam5snpstemporal, file="figs/fig-allele_climate_lfmm_adaptive-AT2G27030-CAM5-warm-vs-cold-trajectories.png",base_height = 7,base_width = 3.5)


### Add the gene models below
# Cam5 gene models
tairraw<-read.table("data-external/TAIR10_GFF3_genes_transposons.gff")
tairsub<-
  tairraw %>% 
  dplyr::filter(grepl("AT2G27030",V9))
tairsub
exon_data <- tairsub %>%
  filter(V3 == "exon") %>%
  mutate(start = as.numeric(V4),
         end = as.numeric(V5),
         isoform = V9)  # Extract exon information
camres %>% 
  ggplot(.)+
  geom_point(aes(y=-log10(p), x=pos))+
  geom_vline( xintercept = 11531967)+
  geom_vline( xintercept = 11534358)+
  xlim(c(11531967-500,11534358+500))+
  geom_rect(data=exon_data,  # adding the genes
            aes(xmin = start, xmax = end,
                ymin = -as.numeric(as.factor(isoform)) - 0.4, ymax = -as.numeric(as.factor(isoform)) -0.1
              ), fill = "grey") +
  theme_minimal()


#####*********************************************************************######
##### Read SVP #####
myatgene="AT2G22540"
mygene="SVP"

sub<-read_merged_allele_subset(
  subsetfile = "data-intermediate/merged_hapFIRE_allele_frequency-AT2G22540-SVP.csv",
  headfile="data-intermediate/merged_hapFIRE_allele_frequency-header.csv",
  snpfile="data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-AT2G22540-SVP.csv"
)
head(sub)
dim(sub)

# get LMs just for corroboration
subres<-
  sub %>% 
  dplyr::group_by(snp) %>% 
  dplyr::filter(year==1) %>% 
  dplyr::summarize(r=cor.test(freq-startfreq,bio1,method = "s")$estimate,
                   p=cor.test(freq-startfreq,bio1,method = "s")$p.value) %>% 
  merge(.,sub,by="snp")
subres %>% head
### get the SVP gene
  tairsub<-
    tairraw %>% 
    dplyr::filter(grepl(myatgene,V9))
  tairsub
  exon_data <- tairsub %>%
    filter(V3 == "exon") %>%
    mutate(start = as.numeric(V4),
           end = as.numeric(V5),
           isoform = V9)  # Extract exon information
# minimanhattan  
# genstart=18588545
genstart=tair[,"Position_start"][tair[,"Gene"]==myatgene]
genend=tair[,"Position_end"][tair[,"Gene"]==myatgene]
minimanhattan<-
subres %>% 
  ggplot(.)+
  geom_point(aes(y=-log10(p), x=pos))+
  geom_vline( xintercept = genstart)+
  geom_vline( xintercept = genend)+
  xlim(c(genstart-5000,genend+5000))+
  theme_minimal()
minimanhattan
#minimanhatatn and gene model
minimanhattan+
  geom_rect(data=exon_data,  # adding the genes
            aes(xmin = start, xmax = end,
                ymin = -as.numeric(as.factor(isoform)) - 0.4, ymax = -as.numeric(as.factor(isoform)) -0.1
            ), fill = "grey") 
    


# select to the maximum
submax<-
  subres %>% 
  dplyr::filter(snp=="2_9583019")
  # dplyr::filter(pos>genstart-500 ,pos<genend+500) %>% 
  # dplyr::filter(-log10(p)>11)
  # dplyr::filter(-log10(p) == max(-log10(p)) )
# top Snp from LFMM 2_9583019

# frequency trajectory
submax %>% 
  dplyr::filter(year==1) %>% 
  ggplot(.)+
  geom_point(aes((freq-startfreq),x=bio1),alpha=0.5,color="#75C376")+
  stat_summary(aes((freq-startfreq),x=bio1), size=1, alpha=1,color="#75C376")+
  stat_smooth(aes((freq-startfreq),x=bio1), method='glm', color="#75C376")+
  geom_hline(yintercept = 0,lty="dotted")+
  # scale_y_log10()+
  ylab("p1-p0")+
  xlab("Annual temperature (C)")+
  ylim(c(-0.06,+0.14))+
  theme_minimal()

submax %>% 
  dplyr::filter(bio1>17 | bio1<12 ) %>% # only look at extremes
  dplyr::mutate(warmvscold= bio1>12) %>% 
  # dplyr::filter(year>0) %>% 
  ggplot(.)+
  geom_jitter(aes(y=freq-startfreq, x=year), alpha=0.2,width = 0.025)+
  stat_summary(aes(y=freq-startfreq, x=year, group=site))+
  stat_smooth(aes(y=freq-startfreq, x=year, group=warmvscold, color=warmvscold), method='glm')+
  # scale_color_manual(values = c( "#4E99C5", "#C4403C"))+
  # geom_hline(lty='dotted',yintercept = 0)+
  # # geom_line(aes(y=freq, x=year, group=site))+
  facet_wrap(~warmvscold, ncol=1)+
  # theme_minimal()+
  ylab("p1-p0")+
  xlab("Year")+
  coord_cartesian(y=c(-0.15,+0.25))+
  # theme(legend.position = "none")
  theme_minimal()


subres %>% 
  dplyr::filter(pos>genstart-5000,pos<genend+5000) %>% 
  dplyr::filter(-log10(p) == max(-log10(p))) %>% 
  merge(.,annotated_snps,by="snp")



#####*********************************************************************######
##### Read TSF #####
myatgene="AT4G20370"
mygene="TSF"

sub<-read_merged_allele_subset(
  subsetfile = "data-intermediate/merged_hapFIRE_allele_frequency-AT4G20370-TSF.csv",
  headfile="data-intermediate/merged_hapFIRE_allele_frequency-header.csv",
  snpfile="data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-AT4G20370-TSF.csv"
)

head(sub)
dim(sub)

# get LMs just for corroboration
subres<-
  sub %>% 
  dplyr::group_by(snp) %>% 
  dplyr::filter(year==3) %>%  # TSF was year 3
  dplyr::summarize(r=cor.test(freq-startfreq,bio1,method = "s")$estimate,
                   p=cor.test(freq-startfreq,bio1,method = "s")$p.value) %>% 
  merge(.,sub,by="snp")
subres %>% head

### get the TSF gene
tairsub<- # just the region
  tairraw %>% 
  dplyr::filter(V1=="4") %>% 
  dplyr::filter(V4 <genstart-1000) %>% 
  dplyr::filter(V5 <genend+5000)

tairsub<-
  tairraw %>% 
  dplyr::filter(grepl(myatgene,V9))

exon_data <- tairsub %>%
  filter(V3 == "exon") %>%
  mutate(start = as.numeric(V4),
         end = as.numeric(V5),
         isoform = V9)  # Extract exon information
# minimanhattan  
# genstart=18588545
genstart=tair[,"Position_start"][tair[,"Gene"]==myatgene]
genend=tair[,"Position_end"][tair[,"Gene"]==myatgene]

# Mini manhattan

minimanhattan<-
  subres %>% 
  ggplot(.)+
  geom_point(aes(y=-log10(p), x=pos))+
  geom_vline( xintercept = genstart)+
  geom_vline( xintercept = genend)+
  xlim(c(genstart-1000,genend+5000))+
  theme_minimal()
minimanhattan

#minimanhatatn and gene model
minimanhattan+
  geom_rect(data=exon_data,  # adding the genes
            aes(xmin = start, xmax = end,
                ymin = -as.numeric(as.factor(isoform)) - 0.4, ymax = -as.numeric(as.factor(isoform)) -0.1
            ), fill = "grey") 


# select to the maximum
submax<-
  subres %>% 
  dplyr::filter(-log10(p) == max(-log10(p)) )
submax %>% head  

# select the maximum within the gene
submax<-
  subres %>% 
  dplyr::filter(pos>genstart ,pos<genend) %>%
  # dplyr::filter(-log10(p)>8)
  dplyr::filter(-log10(p) == max(-log10(p)) )
submax %>% head  

# frequency trajectory
submax %>% 
  dplyr::filter(year==1) %>% 
  ggplot(.)+
  geom_point(aes((freq-startfreq),x=bio1),alpha=0.5,color="#75C376")+
  stat_summary(aes((freq-startfreq),x=bio1), size=1, alpha=1,color="#75C376")+
  stat_smooth(aes((freq-startfreq),x=bio1), method='glm', color="#75C376")+
  geom_hline(yintercept = 0,lty="dotted")+
  # scale_y_log10()+
  ylab("p1-p0")+
  xlab("Annual temperature (C)")+
  ylim(c(-0.06,+0.14))+
  theme_minimal()

# submax %>% 
#   dplyr::filter(bio1>17 | bio1<12 ) %>% # only look at extremes
#   dplyr::mutate(warmvscold= bio1>12) %>% 
#   # dplyr::filter(year>0) %>% 
#   ggplot(.)+
#   geom_jitter(aes(y=freq-startfreq, x=year), alpha=0.2,width = 0.025)+
#   stat_summary(aes(y=freq-startfreq, x=year, group=site))+
#   stat_smooth(aes(y=freq-startfreq, x=year, group=warmvscold, color=warmvscold), method='glm')+
#   # scale_color_manual(values = c( "#4E99C5", "#C4403C"))+
#   # geom_hline(lty='dotted',yintercept = 0)+
#   # # geom_line(aes(y=freq, x=year, group=site))+
#   facet_wrap(~warmvscold, ncol=1)+
#   # theme_minimal()+
#   ylab("p1-p0")+
#   xlab("Year")+
#   coord_cartesian(y=c(-0.15,+0.25))+
#   # theme(legend.position = "none")
#   theme_minimal()


subres %>% 
  dplyr::filter(pos>genstart-5000,pos<genend+5000) %>% 
  dplyr::filter(-log10(p) == max(-log10(p))) %>% 
  merge(.,annotated_snps,by="snp")



save_plot(plot=topsnpplot_cam5, 
          file="figs/fig-allele_climate_lfmm_adaptive-AT2G27030-CAM5.pdf",
          base_height = 3.5,base_width = 3.5)
save_plot(plot=topsnpplot_cam5, 
          file="figs/fig-allele_climate_lfmm_adaptive-AT2G27030-CAM5.png",
          base_height = 3.5,base_width = 3.5)




#####*********************************************************************######
#### Aqua #####
myatgene="AT1G52180"
mygene="Aquaporinlike"

sub<-read_merged_allele_subset(
  subsetfile = "data-intermediate/merged_hapFIRE_allele_frequency-AT1G52180-Aquaporinlike.csv",
  headfile="data-intermediate/merged_hapFIRE_allele_frequency-header.csv",
  snpfile="data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-AT1G52180-Aquaporinlike.csv"
)
head(sub)
dim(sub)

# get LMs just for corroboration
subres<-
  sub %>% 
  dplyr::group_by(snp) %>% 
  dplyr::summarize(r=cor.test(freq-startfreq,bio1)$estimate,
                   p=cor.test(freq-startfreq,bio1)$p.value) %>% 
  merge(.,sub,by="snp")
subres %>% head
# positions of gene  
genstart=tair[,"Position_start"][tair[,"Gene"]==myatgene]
genend=tair[,"Position_end"][tair[,"Gene"]==myatgene]

# minimanhattan
subres %>% 
  ggplot(.)+
  geom_point(aes(y=-log10(p), x=pos))+
  geom_vline( xintercept = genstart)+
  geom_vline( xintercept = genend)+
  xlim(c(genstart-500,genend+500))+
  theme_minimal()

# create top snps and neutral
subrestop<-
  subres %>% 
  dplyr::filter(-log10(p) == max(-log10(p))) #%>% 
  # dplyr::filter(-log10(p)>18) #%>% 
  # dplyr::filter(startfreq<0.5) 

topsnpplot<-
subrestop %>% 
  dplyr::filter(year==1) %>% 
  ggplot(.)+
  geom_point(aes((freq-startfreq),x=bio18, color=bio18),alpha=0.5)+
  stat_summary(aes((freq-startfreq),x=bio18,color=bio18), size=1, alpha=1)+
  stat_smooth(aes((freq-startfreq),x=bio18), method='glm', color="grey")+
  scale_color_gradientn(colors = brewer.pal(9,"BrBG"))+
  geom_hline(yintercept = 0,lty="dotted")+
  # scale_y_log10()+
  ylab("Frequency change (pt+1-pt)")+
  xlab("Summer precipitation (mm)")+
  ylim(c(-0.06,+0.14))+
  theme_minimal()
topsnpplot
save_plot(plot=topsnpplot, 
          file="figs/fig-allele_climate_lfmm_adaptive-AT1G52180-Aquaporinlike.png",
          base_height = 3.5,base_width = 4)
save_plot(plot=topsnpplot, 
          file="figs/fig-allele_climate_lfmm_adaptive-AT1G52180-Aquaporinlike.pdf",
          base_height = 3.5,base_width = 4)

temporalpolot<-
subrestop %>% 
  dplyr::mutate(drywet=bio18>120) %>% 
  # dplyr::filter(year==1) %>%
  ggplot(.)+
  geom_point(aes(y=freq-startfreq,x=year,color=bio18),
             alpha=0.5)+
  stat_summary(aes(y=freq-startfreq,x=year,color=bio18),
               size=1, alpha=1)+
  stat_smooth(aes(y=freq-startfreq,x=year,group=site,color=bio18),
              formula=y~x,se=F,
               size=1, alpha=1, method="glm")+
  scale_color_gradientn(colors = brewer.pal(9,"BrBG"))+
  geom_hline(yintercept = 0,lty="dotted")+
  facet_wrap(~drywet)+
  ylab("Frequency change (pt+1-pt)")+
  xlab("Year")+
  ylim(c(-0.06,+0.14))+
  theme_minimal()
temporalpolot
save_plot(plot=temporalpolot, 
          file="figs/fig-allele_climate_lfmm_adaptive-AT1G52180-Aquaporinlike-temporal.png",
          base_height = 3.5,base_width = 2*3.5)
save_plot(plot=temporalpolot, 
          file="figs/fig-allele_climate_lfmm_adaptive-AT1G52180-Aquaporinlike-temporal.pdf",
          base_height = 3.5,base_width = 2*3.5)


#####*********************************************************************######
#### Dormancy guy #####
myatgene="AT4G19230"
mygene="CYP707A1"

sub<-read_merged_allele_subset(
  subsetfile = paste("data-intermediate/merged_hapFIRE_allele_frequency-",atname,"-",commonname,".csv",sep=""),
  headfile="data-intermediate/merged_hapFIRE_allele_frequency-header.csv",
  snpfile=paste("data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-",atname,"-",commonname,".csv", sep="")
)
head(sub)
dim(sub)

# get LMs just for corroboration
subres<-
  sub %>% 
  dplyr::group_by(snp) %>% 
  dplyr::summarize(r=cor.test(freq-startfreq,bio1)$estimate,
                   p=cor.test(freq-startfreq,bio1)$p.value) %>% 
  merge(.,sub,by="snp")
subres %>% head
# positions of gene  
genstart=tair[,"Position_start"][tair[,"Gene"]==myatgene]
genend=tair[,"Position_end"][tair[,"Gene"]==myatgene]

# minimanhattan
subres %>% 
  ggplot(.)+
  geom_point(aes(y=-log10(p), x=pos))+
  geom_vline( xintercept = genstart)+
  geom_vline( xintercept = genend)+
  xlim(c(genstart-500,genend+500))+
  theme_minimal()

# create top snps and neutral
subrestop<-
  subres %>% 
  dplyr::filter(-log10(p) == max(-log10(p))) #%>% 
# dplyr::filter(-log10(p)>18) #%>% 
# dplyr::filter(startfreq<0.5) 

topsnpplot<-
  subrestop %>% 
  dplyr::filter(year==1) %>% 
  ggplot(.)+
  geom_point(aes((freq-startfreq),x=bio18, color=bio18),alpha=0.5)+
  stat_summary(aes((freq-startfreq),x=bio18,color=bio18), size=1, alpha=1)+
  stat_smooth(aes((freq-startfreq),x=bio18), method='glm', color="grey")+
  scale_color_gradientn(colors = brewer.pal(9,"BrBG"))+
  geom_hline(yintercept = 0,lty="dotted")+
  # scale_y_log10()+
  ylab("Frequency change (pt+1-pt)")+
  xlab("Summer precipitation (mm)")+
  ylim(c(-0.06,+0.14))+
  theme_minimal()
topsnpplot

save_plot(plot=topsnpplot,
          file=paste0("figs/fig-allele_climate_lfmm_adaptive-",myatname,"-",mygene,".png"),
          base_height = 3.5,base_width = 4)
save_plot(plot=topsnpplot,
          file=paste0("figs/fig-allele_climate_lfmm_adaptive-",myatname,"-",mygene,".pdf"),
          base_height = 3.5,base_width = 4)

temporalpolot<-
  subrestop %>% 
  dplyr::mutate(drywet=bio18>80) %>% 
  # dplyr::filter(year==1) %>%
  ggplot(.)+
  geom_point(aes(y=freq-startfreq,x=year,color=bio18),
             alpha=0.5)+
  stat_summary(aes(y=freq-startfreq,x=year,color=bio18),
               size=1, alpha=1)+
  stat_smooth(aes(y=freq-startfreq,x=year,group=site,color=bio18),
              formula=y~x,se=F,
              size=1, alpha=1, method="glm")+
  scale_color_gradientn(colors = brewer.pal(9,"BrBG"))+
  geom_hline(yintercept = 0,lty="dotted")+
  facet_wrap(~drywet)+
  ylab("Frequency change (pt+1-pt)")+
  xlab("Year")+
  ylim(c(-0.06,+0.14))+
  theme_minimal()
temporalpolot
save_plot(plot=temporalpolot, 
          file=paste0("figs/fig-allele_climate_lfmm_adaptive-",myatname,"-",mygene,"-temporal.png"),
          base_height = 3.5,base_width = 2*3.5)
save_plot(plot=temporalpolot, 
          file=paste0("figs/fig-allele_climate_lfmm_adaptive-",myatname,"-",mygene,"-temporal.pdf"),
          base_height = 3.5,base_width = 2*3.5)
