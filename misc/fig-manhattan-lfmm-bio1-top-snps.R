################################################################################
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
library(topr)
library(dplyr)
library(tidyr)
################################################################################

# Get LFMM subsetted
path_lfmm = 'data-intermediate/gea/lfmm_full/results/lfmm_bio1_4topr.csv'
lfmm <- read.csv(path_lfmm)
lfmm$SNP<-paste0(lfmm$CHROM,"_",lfmm$POS)

# and all the results
lfmmraw <- read.csv("data-intermediate/gea/lfmm_full/results/lfmm_bio1_results.csv")
lfmmraw<-
  lfmmraw %>%
  mutate(snp=snp_id) %>%
  rename(P=p_value) %>%
  separate(snp_id, c("CHROM","POS")) %>%
  dplyr::select(P,CHROM,POS, snp)


# Get results gemma
path_gemma='data-intermediate/gea/climate_gwas/lmm_gemma/bio1/output/result_4topr.csv'
gemma <- read.csv(path_gemma)

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

################################################################################

toplfmm<-
  lfmm %>% 
  dplyr::filter(P<1e-10) %>% 
  merge(.,annotated_snps,by.x="SNP",by.y="snp")

for( i in 1:nrow(toplfmm)){

  # Region wuth the highest peak
buffer=5000
chr=2
chr=toplfmm[i,"CHROM"]
pos=toplfmm[i,"POS"]
myregion<-paste0(chr,":",pos-buffer , "-", pos+buffer)

png(paste0("figs/fig-climate_gwas_lfmm_zoom-",myregion,".png"),width = 9,height = 9,res = 600, units = "in")
regionplot(list(gemma, lfmm), region=myregion, build = annotations,
           annotate_with_vline = 1e-8, color=c("darkgrey","#41AB5D"),
           legend_labels = c("1001G GEA", "GrENE-net GEA"), show_gene_legend = F,
           gene_color="#646464")
dev.off()
pdf(paste0("figs/fig-climate_gwas_lfmm_zoom-",myregion,".pdf"),width = 9,height = 9)
regionplot(list(gemma, lfmm), region=myregion, build = annotations,
           annotate_with_vline = 1e-8, color=c("darkgrey","#41AB5D"),
           legend_labels = c("1001G GEA", "GrENE-net GEA"), show_gene_legend = F,
           gene_color="#646464")
dev.off()

}