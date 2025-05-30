################################################################################
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
library(topr)
library(dplyr)
library(tidyr)
################################################################################

# Get LFMM subsetted
# path_lfmm = 'data-intermediate/gea/lfmm_full/results/lfmm_bio1_4topr.csv'
# path_lfmm = 'data-intermediate/gea/lfmm_full/results/l'
# lfmm <- read.csv(path_lfmm)
# lfmm$SNP<-paste0(lfmm$CHROM,"_",lfmm$POS)

# and all the results
# lfmmraw <- read.csv("data-intermediate/gea/lfmm_full/results/lfmm_bio1_results.csv")
lfmmraw <- read.csv("data-intermediate/gea/lfmm_last_gen/lfmm_bio1_k16_results.csv")
lfmmraw<-
  lfmmraw %>%
  mutate(snp=snp_id) %>%
  rename(P=p_value) %>%
  separate(snp_id, c("CHROM","POS")) %>%
  dplyr::select(P,CHROM,POS, snp)
head(lfmmraw)

# Get results gemma
path_gemma='data-intermediate/gea/climate_gwas/lmm_gemma/bio1/output/result_4topr.csv'
path_gemma='data-intermediate/gea/climate_gwas/lmm_gemma/bio1/output/results_lmm.csv'
gemma <- read.csv(path_gemma)
gemma<-
  gemma %>% 
  mutate(snp=rs) %>%
  rename(P=p_wald) %>%
  separate(rs, c("CHROM","POS")) %>%
  dplyr::select(P,CHROM,POS, snp)
head(gemma)

# Get SNPs annotated but subsetted
path_annotations = 'data-intermediate/TAIR10_GFF3_genes_transposons_formatted4topr.csv'
annotations <- read.csv(path_annotations)

annotationsraw <- read.csv('data-external//TAIR10_parsedgenes.csv')
annotationsraw<-
  annotationsraw %>%
  rename(chrom=Chromosome,gene_start=Position_start ,gene_end=Position_end, gene_symbol=Gene) %>%
  mutate(biotype="protein_coding_gene", gene_id=gene_symbol)

# Get common genes genes
commonnames<-read.csv("data-external//common-genes-and-identifiers.csv")
commonnames<-
  commonnames %>% 
  merge(.,annotationsraw,by.x="identifier",by.y="gene_id")

# Load annotated_snps
load("data-intermediate/alleles-TAIR10_parsedgenes.rda")


#####*********************************************************************######
#####* Plot region plot TSF
myatgene="AT4G20370"
mygene="TSF"

buffer=5000


myregion<-
  annotationsraw %>% 
  dplyr::filter(gene_id==myatgene)
 
myregion<-
  paste0(myregion$chrom,":",
         myregion$gene_start-buffer , "-",
         myregion$gene_end+buffer)
myregion



myfile=paste0("figs/fig-climate_gwas_lfmm_zoom-",
              paste0(mygene,"-",
                     myatgene),".png")

png(myfile,width = 9,height = 9,res = 600, units = "in")
regionplot(list(gemma, lfmmraw), 
           region=myregion, 
           build = annotations,
           title= paste(myatgene,mygene),
           annotate_with_vline = 1e-3, color=c("darkgrey","#41AB5D"),
           legend_labels = c("1001G GEA", "GrENE-net GEA"), show_gene_legend = F,
           gene_color="#646464")

dev.off()