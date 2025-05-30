# Goal, visualize TSF
# Dormancy guy
# AT4G20370

# myatgene="AT4G20370"
# mygene="TSF"
#### Dormancy guy #####
myatgene<-myatname<-"AT4G19230"
mygene<-commonname<-"CYP707A1"

################################################################################
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)
################################################################################

# library(grene)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
theme_set(theme_minimal())
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
            samples_tmp, by.x=c("site","year","rep"), 
            by.y=c("site","year","rep"), all.x=T
      )
  }
  # END
  return(a)  
}
#####*********************************************************************######
##### EXTRACT ##$#####
# Annotated SNPs
load("data-intermediate/alleles-TAIR10_parsedgenes.rda")
annotated_snps %>% dim
nrow(annotated_snps)

# Define start and end
myalleles <- annotated_snps %>% dplyr::filter(closest_gene==myatname)
start_row <- which(annotated_snps$snp==head(myalleles$snp,1)) - 1000
end_row <- which(annotated_snps$snp==tail(myalleles$snp,1)) + 1000

sed_command <- paste("sed -n '", start_row, ",", end_row, "p' data/merged_hapFIRE_allele_frequency.csv > data-intermediate/merged_hapFIRE_allele_frequency-",atname,"-",commonname,".csv", sep="")
sed_command_locations <- paste("sed -n '", start_row, ",", end_row, "p' data/greneNet_final_v1.1.bim > data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-",atname,"-",commonname,".csv", sep="")
cat(sed_command)
cat(sed_command_locations)


"grep data/merged_hapFIRE_allele_frequency.csv > data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-AT4G19230-CYP707A1-X4_10523654.csv"
#####*********************************************************************######
##### Frequencies #####

load("data-intermediate/alleles-LDpruned-TAIR10_parsedgenes.rda")

tair<-read.csv("data-external/TAIR10_parsedgenes.csv",header=T)

sub<-read_merged_allele_subset(
  subsetfile = paste("data-intermediate/merged_hapFIRE_allele_frequency-",myatname,"-",commonname,".csv",sep=""),
  headfile="data-intermediate/merged_hapFIRE_allele_frequency-header.csv",
  snpfile=paste("data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-",myatname,"-",commonname,".csv", sep="")
)

head(sub)
dim(sub)

table(sub$snp=="X4_10523654")# check whether top snp is here

onlytop<-read.csv("data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-AT4G19230-CYP707A1-4_10523654.csv")

# get LMs just for corroboration
subres<-
  sub %>% 
  dplyr::group_by(snp) %>% 
  dplyr::filter(year==3) %>%  # TSF was year 3
  dplyr::summarize(r=cor.test(freq-startfreq,bio18,method = "k")$estimate,
                   p=cor.test(freq-startfreq,bio18,method = "k")$p.value) %>% 
  merge(.,sub,by="snp")
subres %>% head
subres %>% 
  arrange(p) %>% 
  head
subres %>% 
  dplyr::filter(p==min(p)) %>% 
  head
genstart=tair[,"Position_start"][tair[,"Gene"]==myatgene]
genend=tair[,"Position_end"][tair[,"Gene"]==myatgene]

subres$snp %>% unique

# get the top 10 snps within the gene
subsummary<-
  subres %>% 
  dplyr::group_by(snp) %>% 
  dplyr::summarise(p=min(p)) %>% 
  arrange(p) %>% 
  head()

# minimanhattan
minimanhattan<-
subres %>% 
  ggplot(.)+
  geom_point(aes(y=-log10(p), x=pos))+
  geom_vline( xintercept = genstart)+
  geom_vline( xintercept = genend)+
  xlim(c(genstart-5000,genend+5000))+
  theme_minimal()
minimanhattan

# create top snps and neutral
subrestop<-
  subres %>% 
  dplyr::mutate(SNP=snp) %>% 
  tidyr::separate(SNP,c("Chr","Pos")) %>% 
  dplyr::filter(Pos>genstart-2800 ,Pos<genend+5000) %>%  # also around
  dplyr::filter(Pos>genstart ,Pos<genend) %>%            # only within
  dplyr::filter(-log10(p) == max(-log10(p))) #%>% 
# dplyr::filter(-log10(p)>18) #%>% 
# dplyr::filter(startfreq<0.5) 
# 4_10521786

# ANNUALPREC
topsnpplot<-
  subrestop %>% 
  group_by(site) %>% 
  dplyr::filter(year==max(year)) %>% 
  # dplyr::filter(year==3) %>% 
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

# ANNUAL TEMP
subrestop %>% 
  group_by(site) %>% 
  dplyr::filter(year==max(year)) %>% 
  # dplyr::filter(year==3) %>% 
  ggplot(.)+
  geom_point(aes((freq-startfreq),x=bio1, color=bio1),alpha=0.5)+
  stat_summary(aes((freq-startfreq),x=bio1,color=bio1), size=1, alpha=1)+
  stat_smooth(aes((freq-startfreq),x=bio1), method='glm', color="grey")+
  scale_color_gradientn(colors = brewer.pal(9,"Reds"))+
  geom_hline(yintercept = 0,lty="dotted")+
  # scale_y_log10()+
  ylab("Frequency change (pt+1-pt)")+
  xlab("Annual temp (mm)")+
  ylim(c(-0.06,+0.14))+
  theme_minimal()

subresmax<-
  subres %>% 
  dplyr::filter(p==min(p))
subresmax$snp %>% unique
# 4:10513686
# 4_10514658


subres %>% 
  dplyr::filter(snp=="4_10514658") %>% 
  # group_by(site) %>% 
  # dplyr::filter(year==max(year)) %>% 
  # dplyr::filter(year==3) %>% 
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



temporalpolot<-
  subrestop %>% 
  dplyr::mutate(drywet=bio18>80) %>% 
  dplyr::filter(snp=="4_10523654") %>% 
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



#####*********************************************************************######
##### VCF ######

# Accessions info
load("grene/data/worldclim_ecotypesdata.rda")
g231<-worldclim_ecotypesdata$ecotypeid
worldclim_ecotypesdata<-read.csv("data-external/bioclimvars_ecotypes_era5.csv") %>% 
  rename(ecotypeid=ecotype) %>% 
  dplyr::filter(ecotypeid%in%g231)
worldclim_ecotypesdata %>% head
# Define individuals to keep (use actual sample names from your VCF file)
individuals_to_keep <- paste(worldclim_ecotypesdata$ecotypeid,worldclim_ecotypesdata$ecotypeid,sep = "_")
length(individuals_to_keep)

# Read VCF
library("vcfR")
vcf_file<-"data-external/arabidopsistranscriptomes/subset_CYP707A1.recode.vcf"
  # old one from 1001G
vcf_file<-"data-intermediate//subset_CYP707A1.recode.vcf"
  # the assembled VCF for grene-net
vcf<-my_vcf <- read.vcfR(vcf_file)
# vcf_subset <- vcf 
# vcf_subset@gt <- vcf@gt[, c(1,which(colnames(vcf@gt) %in% individuals_to_keep))  ] # the 1 is to keep the first column format

as.numeric(vcf@fix[, "POS"]) %>% unique
as.numeric(vcf@fix[, "CHROM"]) %>% unique
range(positions, na.rm = TRUE)

# gene model
source("functions-gene-model-plot.R")
genmod<-gene_model_plot(myatname)

# # Sorted samples based on bio1
# our_popmap<-data.frame( ind=colnames(vcf_subset@gt)[-1],
#                         pop= NA) %>% 
#   # mutate(ecotypes = ind)%>%
#   separate(col = ind, into = c("id2","ecotypeid"),sep = "_") %>% 
#   merge(.,worldclim_ecotypesdata, by.x="ecotypeid",by.y="ecotypeid", all.x=T) %>% 
#   # mutate(pop=rank(bio1)) %>%
#   mutate(pop=round(bio1_new)) %>%
#   mutate(ind = paste0(id2,"_",id2)) %>% 
#   select(ind, pop)
# our_popmap <- our_popmap[order(our_popmap$pop, our_popmap$ind), ]
# head(our_popmap)

# Sorted samples based on bio1
our_popmap<-data.frame( ind=colnames(vcf@gt)[-1],
                        pop= NA) %>% 
  # mutate(ecotypes = ind)%>%
  # separate(col = ind, into = c("id2","ecotypeid"),sep = "_") %>% 
  merge(.,worldclim_ecotypesdata, by.x="ind",by.y="ecotypeid", all.x=T) %>% 
  # mutate(pop=rank(bio1)) %>%
  mutate(pop=round(bio1_new)) %>%
  # mutate(ind = paste0(id2,"_",id2)) %>% 
  select(ind, pop)
our_popmap <- our_popmap[order(our_popmap$pop, our_popmap$ind), ]
head(our_popmap)

# Plot VCF
library(GenotypePlot)
new_plot2 <- genotype_plot(vcf_object  =  vcf_subset,
                           popmap = our_popmap,
                           snp_label_size = 10000,
                           plot_allele_frequency=F,
                           cluster        = FALSE,
                           colour_scheme=c("#d4b9da","#e7298a","#980043")
)


# Make a grid
vcfplot<-
  plot_grid(
    genmod,
    new_plot2$positions,
    new_plot2$genotypes,
    ncol=1,align = "v",
    rel_heights = c(1,1,5)
  )
vcfplot

save_plot(
  paste0("figs/fig-CYP-vcf-genemodel.pdf"),
  vcfplot,
  base_width = 7,base_height = 9 )

save_plot(
  paste0("figs/fig-CYP-vcf-genemodel.png"),
  vcfplot,
  base_width = 7,base_height = 9 )


#####*********************************************************************######
# Extract genotype matrix
source("functions-vcf-genotypes-conveniences.R")
G<-extract_genotypes_vcf(vcf_subset)
G[1:5,1:5]
G$id<- rownames(G)
# G<-G %>% 
#   mutate(idtmp=id) %>% 
#   separate(col = "idtmp",into = c("ecotypeid","fam"))
# G$ecotypeid

G$X4_10514658
G$X4_10514658 %>% table
G$id  
G$ecotypeid

#####*********************************************************************######
#### Phenotypes #######

pheno<-read.csv("data-external/atlas1001_phenotype_matrix_imputed_onlypheno_NEW.csv")
pheno<-read.csv("~/natvar/atlas_phenotype_matrix_withid.csv")
interstcolumns<-c("id","FT16","FT10")
pheno$DSDS90
interstcolumns<-c("id",
                  # "ABA",
                  "ABA_96h_low_water_potential",
                  "DSDS50","FT10",
                  "DSDS10","DSDS90","base_perc",
                  "d4_4C_perc","d4_10C_perc",
                  "d32_10C_perc","d32_4C_perc"
                  ) 
phenosub<-pheno[,interstcolumns] %>% 
  rename(ABA=ABA_96h_low_water_potential)

# Merge VCF with phenotypes  
Gm<-merge(G, phenosub,by.x="id",by.y='id')
G$X4_10514658 %>% table

G$X4_10523654 %>% table
# this is from tati
ggplot(Gm)+
  geom_boxplot(aes(y=base_perc, x=factor(X4_10523654)))+
  geom_point(aes(y=base_perc, x=factor(X4_10523654)))
ggplot(Gm)+
  geom_boxplot(aes(y=d32_4C_perc, x=factor(X4_10523654)))+
  geom_point(aes(y=d32_4C_perc, x=factor(X4_10523654)))


 "#8C510A" "#BF812D" "#DFC27D" "#F6E8C3" "#F5F5F5" "#C7EAE5" "#80CDC1" "#35978F" "#01665E"

 fig_cypdormancy<-
  ggplot(Gm)+
    # geom_violin(aes(y=FT10,x=factor(X4_10999638)))+
    geom_boxplot(aes(y=base_perc,x=factor(X4_10523654)))+
    geom_jitter(aes(y=base_perc,x=factor(X4_10523654), color=factor(X4_10523654)), alpha=0.5, shape=16)+
    # colour_scheme=c("#d4b9da","#e7298a","#980043")
    scale_color_manual("",values = c(
      "#35978F",
      "#8C510A"
      # "#d4b9da",
      # "#e7298a",
      # "#980043"
      ), guide="none")+
    ylab("Germination proportion (%)")+
    xlab("SNP 4_10523654")+
    theme_minimal()+
    theme( guide = "none")
    
  fig_cypdormancy
save_plot("figs/fig-CYP-snps-dormancy.pdf",fig_cypdormancy, 
          base_height = 5,base_width = 6)
save_plot("figs/fig-CYP-snps-dormancy.png",fig_cypdormancy, 
          base_height = 5,base_width = 6)



#### Climate #######
load('grene/data/worldclim_ecotypesdata.rda')
load('grene/data/ecotypes_data.rda')

Gc<-merge(G, worldclim_ecotypesdata,by.x="id",by.y='ecotypeid') %>% 
  merge(ecotypes_data,by.x="id",by.y='ecotypeid')
  
fig_cypdormancylocation<-
ggplot(Gc)+
  geom_boxplot(aes(y=latitude,x=factor(X4_10523654)))+
  geom_jitter(aes(y=latitude,x=factor(X4_10523654), color=factor(X4_10523654)), alpha=0.5, shape=16)+
  scale_color_manual("",values = c(
    "#35978F",
    "#8C510A"
    # "#d4b9da",
    # "#e7298a",
    # "#980043"
  ), guide="none")+
  ylab("Latitude (N)")+
  xlab("SNP 4_10523654")+
  theme_minimal()+
  theme( guide = "none")
  # colour_scheme=c("#d4b9da","#e7298a","#980043")

fig_cypdormancyclimate<-
ggplot(Gc)+
  geom_boxplot(aes(y=bio18,x=factor(X4_10523654)))+
  geom_jitter(aes(y=bio18,x=factor(X4_10523654), color=factor(X4_10523654)), alpha=0.5, shape=16)+
  scale_color_manual("",values = c(
    "#35978F",
    "#8C510A"
    # "#d4b9da",
    # "#e7298a",
    # "#980043"
  ), guide="none")+
  ylab("Summer precipitation (mm)")+
  xlab("SNP 4_10523654")+
  theme_minimal()+
  theme( guide = "none")

sink("tables/test-correlation-CYP-SNP-with-dormancy.txt")

t.test(Gm$base_perc[Gm$X4_10523654==1], Gm$base_perc[Gm$X4_10523654==0])
wilcox.test(Gm$base_perc[Gm$X4_10523654==1], Gm$base_perc[Gm$X4_10523654==0])


t.test(Gc$latitude[Gc$X4_10523654==1], Gc$latitude[Gc$X4_10523654==0])
wilcox.test(Gc$latitude[Gc$X4_10523654==1], Gc$latitude[Gc$X4_10523654==0])

t.test(Gc$bio18[Gc$X4_10523654==1], Gc$bio18[Gc$X4_10523654==0])
wilcox.test(Gc$bio18[Gc$X4_10523654==1], Gc$bio18[Gc$X4_10523654==0])

sink()


figdorm<-plot_grid(fig_cypdormancyclimate,fig_cypdormancylocation,fig_cypdormancy, ncol=3)
figdorm
save_plot("figs/fig-CYP-snps-dormancy-dormancy-climate-location-evidence.pdf",figdorm, 
          base_height = 5,base_width = 7)
save_plot("figs/fig-CYP-snps-dormancy-dormancy-climate-location-evidence.png",figdorm, 
          base_height = 5,base_width = 7)


# # TOP SNP visualiy
# ggplot(Gm)+
#   geom_boxplot(aes(y=ABA, x=factor(X4_10514658)))+
#   geom_point(aes(y=ABA, x=factor(X4_10514658)))
# ggplot(Gm)+
#   geom_boxplot(aes(y=d32_4C_perc, x=factor(X4_10514658)))+
#   geom_point(aes(y=d32_4C_perc, x=factor(X4_10514658)))
# ggplot(Gm)+
#   geom_boxplot(aes(y=d4_4C_perc, x=factor(X4_10514658)))+
#   geom_point(aes(y=d4_4C_perc, x=factor(X4_10514658)))
# 
# ggplot(Gm)+
#   geom_point(aes(y=d4_4C_perc, x=d32_4C_perc, color=factor(X4_10514658)))
# 
#   
# ggplot(Gm)+
#   geom_boxplot(aes(y=DSDS50, x=factor(X4_10514658)))+
#   geom_point(aes(y=DSDS50, x=factor(X4_10514658)))
# ggplot(Gm)+
#   geom_boxplot(aes(y=DSDS10, x=factor(X4_10514658)))+
#   geom_point(aes(y=DSDS10, x=factor(X4_10514658)))
# ggplot(Gm)+
#   geom_boxplot(aes(y=DSDS90, x=factor(X4_10514658)))+
#   geom_point(aes(y=DSDS90, x=factor(X4_10514658)))
# ggplot(Gm)+
#   geom_boxplot(aes(y=base_perc, x=factor(X4_10514658)))+
#   geom_point(aes(y=base_perc, x=factor(X4_10514658)))
# ggplot(Gm)+
#   geom_boxplot(aes(y=FT10, x=factor(X4_10514658)))+
#   geom_point(aes(y=FT10, x=factor(X4_10514658)))


t.test(Gm$ABA[Gm$X4_10514658==1], Gm$ABA[Gm$X4_10514658==0])
wilcox.test(Gm$ABA[Gm$X4_10514658==1], Gm$ABA[Gm$X4_10514658==0])


t.test(Gm$base_perc[Gm$X4_10514658==1], Gm$base_perc[Gm$X4_10514658==0])
wilcox.test(Gm$base_perc[Gm$X4_10514658==1], Gm$base_perc[Gm$X4_10514658==0])

t.test(Gm$DSDS50[Gm$X4_10514658==1], Gm$DSDS50[Gm$X4_10514658==0])
wilcox.test(Gm$DSDS50[Gm$X4_10514658==1], Gm$DSDS50[Gm$X4_10514658==0])


# [1] "4_10514108" ,"4_10514220", "4_10514456" ,"4_10515806" ,"4_10517772", "4_10518057", "4_10519863" ,"4_10520431" ,"4_10524150"
c("4_10514108" ,"4_10514220", "4_10514456" ,"4_10515806" ,"4_10517772", "4_10518057", "4_10519863" ,"4_10520431" ,"4_10524150")

Gm$poly<-
  apply(Gm[,c("X4_10514108" ,"X4_10514220", "X4_10514456" ,"X4_10515806" ,"X4_10517772", "X4_10518057", "X4_10519863" ,"X4_10520431" ,"X4_10524150")], 1, sum)
          

ggplot(Gm)+
  geom_boxplot(aes(y=DSDS50, x=factor(X4_10514658)))+
  geom_point(aes(y=DSDS50, x=factor(X4_10514658)))

ggplot(Gm)+
  geom_boxplot(aes(y=DSDS50, x=factor(poly)))+
  geom_point(aes(y=DSDS50, x=(poly)))




