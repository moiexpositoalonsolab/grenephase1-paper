# Goal, visualize TSF
# TWIN SISTER OF FT (TSF)
# AT4G20370

myatgene="AT4G20370"
mygene="TSF"

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
##### Frequencies #####

load("data-intermediate/alleles-LDpruned-TAIR10_parsedgenes.rda")

tair<-read.csv("data-external/TAIR10_parsedgenes.csv",header=T)


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
  # dplyr::filter(year==3) %>%  # TSF was year 3
  dplyr::filter(year==1) %>%  # TSF was year 3
  dplyr::summarize(r=cor.test(freq-startfreq,bio1,method = "k")$estimate,
                   p=cor.test(freq-startfreq,bio1,method = "k")$p.value) %>%
  merge(.,sub,by="snp") %>%
  mutate(snpchange=snp) %>%
  separate(snpchange,c("chr","pos"))

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

subrestop<-
  subres %>%
  mutate(snpchange=snp) %>%
  separate(snpchange,c("chr","pos"))


  dplyr::filter(pos<11000000) %>%
  dplyr::filter(pos>10997000) %>%
  dplyr::filter(p==min(p))

subrestop %>% dim
table(subrestop$pos>10997000)
table(subrestop$pos>10997000 & subrestop$pos<11000000)

# Minimanhattan
ggplot(subrestop)+
  geom_point(aes(y=-log10(p), x=pos))

# get the top 10 snps within the region, but downstream of the gene
subsummary<-
  subres %>%
  dplyr::filter(pos<11002500) %>%
  # dplyr::filter(pos<11000000) %>%
  dplyr::group_by(snp) %>%
  dplyr::dplyr(p==min(p,na.rm=T)) %>%
  arrange(p) %>%
  head()


# Minihamnhattan
subres %>%
  ggplot(.)+
  geom_point(aes(y=-log10(p), x=as.numeric(pos)))+
  # geom_vline( xintercept = genstart)+
  # geom_vline( xintercept = genend)+
  xlim(c(genstart-5000,genend+5000))+
  theme_minimal()

subres %>%
  dplyr::filter(snp=="4_10999396") %>%  # identified by tati
  ggplot()+
  geom_point(aes(y=freq,x=year,color=bio1))+
  stat_smooth(aes(y=freq,x=year,color=bio1, group=snp),method="glm")+
  scale_color_gradientn(colours = rev(brewer.pal(9,"RdBu")))+
  facet_grid(~bio1>10)

subres %>%
  dplyr::filter(snp=="4_10999396") %>% # identified by tati
  ggplot()+
  geom_point(aes(y=freq-startfreq,x=bio1,color=bio1))

# The one from Leet et al
subres %>%
  dplyr::filter(snp=="4_10999188") %>%
  ggplot()+
  geom_point(aes(y=freq,x=year,color=bio1))+
  stat_smooth(aes(y=freq,x=year,color=bio1, group=snp),method="glm")+
  scale_color_gradientn(colours = rev(brewer.pal(9,"RdBu")))+
  facet_grid(~bio1>10)



# The one in aragwas
# https://aragwas.1001genomes.org/#/gene/AT4G20370

subres %>%
  # dplyr::filter(snp=="4_10999096") %>%
  dplyr::filter(snp=="4_10999188") %>%
  ggplot()+
  geom_point(aes(y=freq,x=year,color=bio1))+
  stat_smooth(aes(y=freq,x=year,color=bio1, group=snp),method="glm")+
  scale_color_gradientn(colours = rev(brewer.pal(9,"RdBu")))+
  facet_grid(~bio1>12)

subres %>%
  # dplyr::filter(snp=="4_10999096") %>%
  # dplyr::filter(p==min(p)) %>%
  # dplyr::filter(snp=="4_10999188") %>%
  dplyr::filter(snp=="4_11004863") %>%
  ggplot()+
  geom_point(aes(y=1-(freq-startfreq),x=bio1,color=bio1))+
  stat_smooth(aes(y=1-(freq-startfreq),x=bio1,color=bio1), method="glm")+
  scale_color_gradientn(colours = rev(brewer.pal(9,"RdBu")))#+


# subres %>%  # For visualization the gradient
#   dplyr::filter(pos<11000000) %>%
#   ggplot()+
#   geom_point(aes(y=freq-startfreq,x=bio1,color=bio1))+
#   # stat_smooth(aes(y=freq-startfreq,x=bio1),method="glm")+
#   scale_color_gradientn(colours = rev(brewer.pal(9,"RdBu")))#+
#   # facet_grid(~bio1>10)

subrestop %>%  # For visualization the gradient
  ggplot()+
  geom_point(aes(y=freq-startfreq,x=bio1,color=bio1))+
  stat_smooth(aes(y=freq-startfreq,x=bio1),method="glm")+
  scale_color_gradientn(colours = rev(brewer.pal(9,"RdBu")))#+
# facet_grid(~bio1>10)


# This one had pretty high p value 4_10992697

#####*********************************************************************######
##### TRANSCRIPTOMES #####
################################################################################
# Read backbone fam
# fam<-read.table("/2029g/2029g.fam",header=F)
# Load normalized expression data
load("data-external/arabidopsistranscriptomes/1001t/TG_data_20180606.Rdata")
load("data-external/arabidopsistranscriptomes/1001t/gene_infoV2.Rdata")

# remove batch MU, which was the original low quality
ind2use<- which(TG.meta$batch_comb != "MU")
# subset log2 expression for all good individuals
d<-TG.genes$d_log2_batch[,ind2use]
# Generate an expression matrix (individuals=rows, genes=cols)
d<-data.frame(TG.meta$index[ind2use], t(d))
colnames(d)<-c("id",gene_infoV2$Name)

##### Ectract
d2<-TG.trans$d_log2[,ind2use]
# d2<-TG.trans$raw[,ind2use]
# d<-TG.genes$d_log2_batch[,ind2use]
# Generate an expression matrix (individuals=rows, genes=cols)
d2<-data.frame(TG.meta$index[ind2use], t(d2))
colnames(d2)<-c("id",TG.trans$names)

#####*********************************************************************######
##### Try to find correlation ########
# Try to find CAM5 in the expression data
# genes
genesubset<-d %>% dplyr::select(id,AT4G20370,AT1G65480) # TSF and FT
genesubset
ggplot(genesubset)+
  geom_point(aes(y=AT4G20370, x=AT1G65480))+
  xlab("FT")+ylab("TSF")
genesubset[,2] %>% is.na %>% table

# transcripts
transcriptsubset<-d2 %>% data.frame() %>%  dplyr::select(id,AT4G20370.1, AT1G65480.1,AT1G65480.2 )
  head(transcriptsubset)
ggplot(transcriptsubset)+
  geom_point(aes(y=AT4G20370.1, x=AT1G65480.1))+
  xlab("FT")+ylab("TSF")

#####*********************************************************************######
#### Worldclim #######
worldclim_ecotypesdata<-read.csv("data-external/1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv")
head(worldclim_ecotypesdata)

# Merge
mer<-merge(genesubset,worldclim_ecotypesdata, by.x='id',by.y="ecotypeid")
mer2<-merge(transcriptsubset,worldclim_ecotypesdata, by.x='id',by.y="ecotypeid")


# Plot cam5 and temp
ggplot(mer)+
  geom_point(aes(y=AT4G20370, x=bio1))+
  stat_smooth(aes(y=AT4G20370, x=bio1),method='glm')

#  splice forms
transcript1<-
  ggplot(mer2)+
  geom_point(aes(y=AT4G20370.1, x=bio1))+
  stat_smooth(aes(y=AT4G20370.1, x=bio1),method='glm', color= "#2166AC")+
  labs(x="Annual temperature of accession (C)")
transcript1

cor.test(mer$AT4G20370, mer$bio1, method='s')
cor.test(mer2$AT4G20370.1, mer$bio1, method='p')
cor.test(mer2$AT4G20370.1, mer$bio1, method='s')

#####*********************************************************************######
#### Phenotypes #######
pheno<-read.csv("data-external/atlas1001_phenotype_matrix_imputed_onlypheno_NEW.csv")
pheno<-read.csv("data-external/atlas_phenotype_matrix_withid.csv")
interstcolumns<-c("id","FT16","FT10")
phenosub<-pheno[,interstcolumns]

head(phenosub)

mm2<-merge(transcriptsubset,phenosub, by.x='id',by.y="id", all=T)

ggplot(mm2)+
  geom_point(aes(y=FT10,x=AT4G20370.1))+
  ylab("Flowering time (d)")+xlab("TSF")
ggplot(mm2)+
  geom_point(aes(y=FT10,x=AT1G65480.1))+
  ylab("Flowering time (d)")+xlab("FT")

# with TSF
cor.test(mm2$FT10,mm2$AT4G20370.1, method="s")
cor.test(mm2$FT10,mm2$AT4G20370.1, method="s")$p.value
# cor.test(mm2$FT10,mm2$AT4G20370.1, method="p")
# cor.test(mm2$FT10,mm2$AT4G20370.1, method="k")

# with FT
cor.test(mm2$FT10,mm2$AT1G65480.1, method="s")
cor.test(mm2$FT10,mm2$AT1G65480.1, method="s")$p.value
# cor.test(mm2$FT10,mm2$AT1G65480.1, method="p")
# cor.test(mm2$FT10,mm2$AT1G65480.1, method="k")


#####*********************************************************************######
##### SNP and expression ########

# Read VCF
library(vcfR)
vcf_file<-"data-intermediate/subset_TSF.recode.vcf"
vcf_file<-"data-external/arabidopsistranscriptomes/subset_TSF.recode.vcf"
vcf<-my_vcf <- read.vcfR(vcf_file)
vcf_subset <- vcf
# vcf_subset@gt <- vcf@gt[, c(1,which(colnames(vcf@gt) %in% individuals_to_keep))  ] # the 1 is to keep the first column format


gt <- extract.gt(vcf_subset, element = 'GT', as.numeric = TRUE) %>%
  t() %>%
  data.frame
dim(gt)
head(gt)
gt[1:5,1:5]

gt$X4_10999638 %>% table

# are the freq data overlap with vcf?
fc<-function (data.frame)
{
  as.character(as.matrix(data.frame))
}
table(fc(unique(subres$snp))  %in% fc(unique(vcf@fix[,3])))
table(fc(unique(vcf@fix[,3])) %in% fc(unique(subres$snp)) )

#
myinds<- row.names(gt)
myinds<-fc(sapply(myinds, function(i)strsplit(i,split = "_",fixed = T)[[1]][1]))
gt$id<-myinds

# get the top SNPs in freq environment correlation present in vcf
subsummary<-
  subres %>%
  dplyr::filter(snp %in% vcf@fix[,3]) %>%
  dplyr::group_by(snp) %>%
  dplyr::summarise(p=min(p)) %>%
  arrange(p) %>%
  head()
# Subset genotype table with top association snps
gtsub<-gt[,colnames(gt) %in% paste0("X",subsummary$snp)]
dim(gtsub)
gtsub$id<-myinds

# NO SUBSET OF ASSOCIATIONS
gtsub<-gt

# merge VCF with expression andphenotype data
phenosub$id<-as.numeric(phenosub$id)
gtsub$id<-as.numeric(gtsub$id)
gtm<-
  merge(gtsub, phenosub, by="id", all=T)

head(gtm)
dim(gtm)


# Plot this hit X4_10999638 is the one identified in them mini GWA above with correlation
gtm %>%
  dplyr::select(FT10, X4_10999638) %>%
  na.omit() %>%
ggplot(.)+
  # geom_violin(aes(y=FT10,x=factor(X4_10999638)))+
  geom_boxplot(aes(y=FT10,x=factor(X4_10999638)))+
  geom_jitter(aes(y=FT10,x=factor(X4_10999638), color=factor(X4_10999638)), alpha=0.5)+
  # colour_scheme=c("#d4b9da","#e7298a","#980043")
  scale_color_manual("",values = c(
                              "#d4b9da",
                              # "#e7298a",
                              "#980043"))+
  ylab("Flowering time in growth chamber (10C)")+
  xlab("SNP 4_10999638")+
  theme_minimal()->
    fig_tsf_floweringtime
fig_tsf_floweringtime

# from LFMM plot the SNP is X4_10999396!
gtm %>%
  dplyr::select(FT10, X4_10999396) %>%
  na.omit() %>%
  ggplot(.)+
  geom_boxplot(aes(y=FT10,x=factor(X4_10999396)))+
  geom_jitter(aes(y=FT10,x=factor(X4_10999396), color=factor(X4_10999396)), alpha=0.5)+
  # colour_scheme=c("#d4b9da","#e7298a","#980043")
  scale_color_manual("",values = c(
    "#d4b9da",
    # "#e7298a",
    "#980043"))+
  ylab("Flowering time in growth chamber (10C)")+
  xlab("SNP 4_10999638")+
  theme_minimal()->
  fig_tsf_floweringtime
fig_tsf_floweringtime

save_plot("figs/fig-TSF-snps-floweringtime.pdf",fig_tsf_floweringtime,
          base_height = 5,base_width = 6)
save_plot("figs/fig-TSF-snps-floweringtime.png",fig_tsf_floweringtime,
          base_height = 5,base_width = 6)


sink("tables/test-correlation-TSF-SNP-with-flowering-time-and-transcript.txt")

mean( gtm$FT10[gtm$X4_10999638==1] , na.rm = T) - mean( gtm$FT10[gtm$X4_10999638==0] , na.rm = T)

table(gtm$X4_10999638)
wilcox.test(gtm$FT10~gtm$X4_10999638)
wilcox.test(gtm$AT4G20370.1~gtm$X4_10999638)

sink()

wilcox.test(gtm$FT10~gtm$X4_10999396)
t.test(gtm$FT10~gtm$X4_10999396)

# merge VCF with geographic and climateic locations
gtmer2<-
  merge(gtsub, worldclim_ecotypesdata, by.x="id",by.y="ecotypeid", all=T)

# Plot to see where the alleles are
ggplot(gtmer2)+
  geom_point(aes(x = bio1, y=bio12, color=factor(X4_10999638) ))+
  scale_color_manual("",values = c(
    "#d4b9da",
    # "#e7298a",
    "#980043"))+
  theme_minimal()

#
#
# ggplot(gtm)+
#   geom_point(aes(y=FT16,x=AT4G20370.1,
#             size=X4_10999638,
#             color=factor(X4_10999638) ), # top SNP
#             # color=factor(X4_10998998)),
#             # color=factor(X4_10998766)),
#             # color=factor(X4_10998495)),
#             # color=factor(X4_10999638)),
#              # size=2
#              )+
#   scale_color_manual(values = c("#d4b9da",
#                                 # "#e7298a",
#                                 "#980043"))+
#   ylab("Flowering time (d)")+xlab("TSF")

