################################################################################
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
library(topr)
library(dplyr)
library(tidyr)

################################################################################

lfmm<-read.csv("data-intermediate/gea/lfmm_gen1/bio1/lfmm_bio1_k25_results.csv",header=T)
lfmm<-read.csv("data-intermediate/gea/lfmm_full//results/lfmm_bio1_results.csv",header=T) # old one

wa<-read.csv("data-intermediate/gea/lfmm_gen1/bio1/wza_results_lfmm_bio1.csv",header = T)
head(wa)
wa$pos<-1:nrow(wa)
wa<-
  wa %>%  mutate(chr = str_extract(gene, "^[^_]+"))

waplot<-
  ggplot(wa)+
  geom_point(aes(y=-log10(Z_pVal),x=pos, color=chr))

waplot+
  geom_point(data=dplyr::filter(wa,gene=="5_2922"),
             aes(y=-log10(Z_pVal),x=pos),
             color="green"
             )


waplot+
  geom_point(data=dplyr::filter(wa,gene=="2_970"), # this is svp
             aes(y=-log10(Z_pVal),x=pos),
             color="green"
  )

lfmm %>% head
lfmmtmp<-
lfmm %>% 
  dplyr::filter(block=="2_970") %>% 
  dplyr::mutate(snp=snp_id) %>% 
  separate(snp_id,into = c("chr","pos"),sep = "_") %>% 
  dplyr::mutate(chr=as.numeric(chr), pos=as.numeric(pos))

ggplot(lfmmtmp)+
  geom_point(
             aes(y=-log10(p_value),x=pos),
             color="green"
  )+
  geom_vline( xintercept = genstart)+
  geom_vline( xintercept = genend)+
  # xlim(c(genstart-5000,genend+5000))+
  theme_minimal()

lfmmtmp %>%
  dplyr::filter(-log10(p_value)>6.5) %>%
    merge(.,annotated_snps,by="snp") %>% 
  dplyr::filter(closest_gene=="AT2G22540")


"AGL17"
"AT2G22630"

tair %>% 
  dplyr::filter(Gene=="AT2G22630")
tair %>% 
  dplyr::filter(Gene=="AT2G22540")



ggplot(lfmmtmp)+
  geom_point(
    aes(y=-log10(p_value),x=pos),
    color="green"
  )+
  # geom_vline( xintercept = genstart)+
  # geom_vline( xintercept = genend)+
  
  geom_vline( xintercept = 9618372)+
  geom_vline( xintercept = 9621957)+
  theme_minimal()
  
# tmp2<-
#   lfmmtmp %>% 
#   dplyr::filter(-log10(p_value)>8) %>% 
#   merge(.,annotated_snps,by="snp")
# tmp2
# 
# myalleles <- annotated_snps %>% dplyr::filter(closest_gene==atname)
