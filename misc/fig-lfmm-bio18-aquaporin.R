################################################################################
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
library(topr)
library(dplyr)
library(tidyr)
scale_pvalues <- function(pvalues) {
  # Desired minimum p-value corresponds to -log10(p) = 10, hence p_min = 10^-10
  target_min_p <- 10^-10
  
  # Find the actual minimum p-value in the input vector
  actual_min_p <- min(pvalues)
  
  # Calculate the scaling factor
  scaling_factor <- actual_min_p / target_min_p
  
  # Scale all p-values
  scaled_pvalues <- pvalues / scaling_factor
  
  return(scaled_pvalues)
}


################################################################################

# Get LFMM subsetted
path_lfmm = 'data-intermediate/gea/lfmm_full/results/lfmm_bio18_results.csv'
lfmm <- read.csv(path_lfmm)
head(lfmm)
# Make new snp column
lfmm$SNP<-paste0(lfmm$CHROM,"_",lfmm$POS)
lfmm<-lfmm %>% 
  separate(col = snp_id, into = c("CHROM","POS")) %>% 
  rename(P=p_value)
# lfmm$P<-scale_pvalues(lfmm$P) ### REMOVE SCALING
# 
# 
# # Get results gemma - bio1
# # path_gemma='data-intermediate/gea/climate_gwas/lmm_gemma/bio1/output/result_4topr.csv'
# path_gemma='~/natvar/gwas/gemma/ABA_96h_low_water_potential/output/Verslues_Kalladan_PNAS_2017_PID_29073083.lmm.assoc.txt'
# 
# path_gemma='~/natvar/gwas/gemma/norm_base_perc/output/Schmitt_Martinez-Berdeja_PNAS_2020_PID_foreditMR.lmm.assoc.txt'
# gemma <- read.csv(path_gemma, sep="\t")
# head(gemma)
# gemma <- gemma %>% 
#   # separate(col = snp_id, into = c("CHROM","POS")) %>% 
#   rename(CHROM=chr, POS=ps,P=p_lrt)
# gemma$P<-scale_pvalues(gemma$P)


################################################################################
# ANNOTAIONS
# Get SNPs annotated but subsetted
annotations <- read.csv('data-intermediate/TAIR10_GFF3_genes_transposons_formatted4topr.csv')
annotations %>% head

# annotationsraw <- read.csv('data-external//TAIR10_parsedgenes.csv')
# annotationsraw<-
#   annotationsraw %>%
#   rename(chrom=Chromosome,gene_start=Position_start ,gene_end=Position_end, gene_symbol=Gene) %>%
#   mutate(biotype="protein_coding_gene", gene_id=gene_symbol)

# Load annotated_snps
load("data-intermediate/alleles-TAIR10_parsedgenes.rda")
annotated_snps %>% head

################################################################################
# Exploration of manhattans

savetrue=F
if(savetrue) pdf(file="figs/fig-manhattan-lfmm-bio18.pdf")
manhattan(lfmm, 
          build = annotations,
          # sign_thresh = 1,
          # size=c(1,3),
          color=c("#41AB5D"),
          legend_labels=c("GrENE-net GEA"), verbose=F)
if(savetrue) dev.off()

################################################################################

# Region plot of top gene
mygene="AT1G79780"
mynames<-annotations %>% 
  dplyr::filter(gene_id==mygene)
buffer=5000
myregion<-paste0(mynames[,"chrom"],":",
                 mynames[,"gene_start"]-buffer , "-",
                 mynames[,"gene_end"]+buffer)
myregion



################################################################################

# Region plot of top gene
mygene="AT1G52180"
mynames<-annotations %>% 
  dplyr::filter(gene_id==mygene)
buffer=1000
myregion<-paste0(mynames[,"chrom"],":",
                 mynames[,"gene_start"]-buffer , "-",
                 mynames[,"gene_end"]+buffer)
myregion


savetrue=FALSE
if(savetrue) pdf(file="figs/fig-manhattan-lfmm-bio18-zoomregion-aquaporinlike-AT1G52180.pdf")
regionplot(
  list(lfmm),
  region = myregion,
  build = annotations,
  title = mygene,
  annotate_with_vline = 1e-3, 
  color=c("#41AB5D"), 
  legend_labels=c("GrENE-net GEA") ,
  show_gene_legend = F,
  gene_color="#646464"
) 
if(savetrue) dev.off()
