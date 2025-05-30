################################################################################
### Goal
### Estimate allele climate associations simple GLM

################################################################################
#### If transformed to collab can use the stuff below
# from google.colab import drive
# drive.mount('/content/drive')
# %load_ext rpy2.ipython
# %%R

################################################################################
#### Libraries
library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(patchwork)

################################################################################
#### Location if run locally
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
# myfolder<-"~/grenephase1-analyses"

myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
# load("data-intermediate/allele_allyears_frequencies_long_raw_climate_flowers.rda")
# load("data-intermediate/allele_terminal_frequencies_long_raw_climate.rda")
load("data-intermediate/allele_terminal_frequencies_long_raw_climate_flowers.rda")
a<-allele_terminal_frequencies_long_raw_climate_flowers

# Parse it
a1<- a[grepl("1_", a$snp),]
a2<- a1 %>% 
  mutate(SNP=snp) %>% 
  separate(snp, into = c("chromosome", "position"), sep = "_") %>% 
  mutate(position=as.numeric(position))
  
################################################################################
#### Test a glm gaussian and a glm binomial
asub<-
  a2 %>% 
  dplyr::filter(SNP=="1_10001716")

#######################
# GLM functions
glmfunctionp<-function(testdelta,testyear,testbio){
  lmmod<-lm(testdelta ~ testyear + testbio)
  p=coefficients(summary(lmmod))[3,4]
  b=coefficients(summary(lmmod))[3,1]
  return(p)
}
glmfunctionb<-function(testdelta,testyear,testbio){
  lmmod<-lm(testdelta ~ testyear + testbio)
  p=coefficients(summary(lmmod))[3,4]
  b=coefficients(summary(lmmod))[3,1]
  return(b)
}
# Results
resglm<-
  a12%>%
  group_by(chr_pos, chromosome, position) %>%
  summarise(b_bio1= glmfunctionb(maxfreq, max_year,sites_bio1),
            p_bio1= glmfunctionp(maxfreq, max_year,sites_bio1)
  )
resglm<-
  resglm %>% 
  mutate(position=as.numeric(position))

library(ggman)
ggmanhattan(data = resglm ,chr = "chromosome",bp = "position",P = "p_bio1", SNP = "chr_pos",
            sparsify=FALSE,
            build="NA", # need to remove the default hg19 otherwise map incorrect
            expand.x = c(0.03, 0.1), expand.y = c(0.01,0.03),
            significance=0.05/nrow(resglm)
)
#######################
# GLM binomial functions
glmfunctionpbinomial<-function(testfreq, testcount,testyear,testbio){
  lmmod<-glm(cbind(testfreq,testcount) ~ testyear + testbio,
             family = quasibinomial(link = logit))
  p=drop1(lmmod, test="Chisq")
  p=as.data.frame(p)[3,4]
  return(p)
}
#Example on 1 snp
glmfunctionpbinomial(round(a1sub$maxfreq*a1sub$flowers), 
                     a1sub$flowers-round(a1sub$maxfreq*a1sub$flowers),
                     a1sub$max_year, 
                     a1sub$sites_bio1)
# Results
resglmbinom<-
  a2 %>%
  group_by(SNP, chromosome, position) %>%
  summarise(
            p_bio1= glmfunctionpbinomial(
                    round(maxfreq*flowers), 
                    flowers-round(maxfreq*flowers),
                    max_year,
                    sites_bio1)
           )

resglmbinom %>% head

library(ggman)
manhattanbinom<-
ggmanhattan(data = resglmbinom ,chr = "chromosome",bp = "position",P = "p_bio1", 
            SNP = "SNP",
            sparsify=FALSE,
            build="NA", # need to remove the default hg19 otherwise map incorrect
            expand.x = c(0.03, 0.1), expand.y = c(0.01,0.03),
            significance=0.05/nrow(resglm)
)
save_plot("fig-GEA-binomial-simple-bio1.pdf",manhattanbinom, base_height = 4, base_width = 9)
save_plot("fig-GEA-binomial-simple-bio1-qq.pdf",ggman:::gq(x = resglmbinom$p_bio1), base_height = 4, base_width = 4)