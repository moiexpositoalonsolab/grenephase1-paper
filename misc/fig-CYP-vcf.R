################################################################################
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
library(topr)
library(dplyr)
library(tidyr)
library(vcfR)
library(cowplot)

myatgene<-myatname<-"AT4G19230"
mygene<-commonname<-"CYP707A1"


################################################################################
# Accessions info
load("data-intermediate/alleles-TAIR10_parsedgenes.rda")
head(annotated_snps)

tair<-read.csv("data-external/TAIR10_parsedgenes.csv")
tair<-read.table("data-external/TAIR10_GFF3_genes_transposons.gff")
head(tair)

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

################################################################################
# Extract Gene
# Define your start and end rows as variables

# Define start and end
myalleles <- annotated_snps %>% dplyr::filter(closest_gene==myatname)
chrom=myalleles$Chr %>% head(1)
buffer=2800
start_bp <- min(myalleles$Pos)
end_bp <- max(myalleles$Pos)
start_row <- which(annotated_snps$snp==head(myalleles$snp,1)) - buffer
end_row <- which(annotated_snps$snp==tail(myalleles$snp,1)) + buffer

# Construct the sed command dynamically
sed_command <- paste("sed -n '", start_row, ",", end_row, "p' data/merged_hapFIRE_allele_frequency.csv > data-intermediate/merged_hapFIRE_allele_frequency-",
                     myatname,"-",commonname,".csv", sep="")
sed_command_locations <- paste("sed -n '", start_row, ",", end_row, "p' data/greneNet_final_v1.1.bim > data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-",
                               myatname,"-",commonname,".csv", sep="")
cat(sed_command)
cat(sed_command_locations)

# Extract vcf
vcftools_command_locations <- paste0("vcftools --vcf greneNet_final_v1.1.recode.vcf --chr " ,chrom, " --from-bp ",start_bp," --to-bp ",end_bp ,
                                     " --recode --recode-INFO-all --out subset_" , commonname )
cat(vcftools_command_locations)



# buffer=5000 # The top hit of TSF is actually a few thousands bp away upstream 
# start_bp <- min(myalleles$Pos)
# end_bp <- max(myalleles$Pos)
# start_row <- which(annotated_snps$snp==head(myalleles$snp,1)) - buffer
# end_row <- which(annotated_snps$snp==tail(myalleles$snp,1)) + buffer
# 
# # Construct the sed command dynamically
# sed_command <- paste("sed -n '", start_row, ",", end_row, "p' data/merged_hapFIRE_allele_frequency.csv > data-intermediate/merged_hapFIRE_allele_frequency-",
#                      myatname,"-",commonname,".csv", sep="")
# sed_command_locations <- paste("sed -n '", start_row, ",", end_row, "p' data/greneNet_final_v1.1.bim > data-intermediate/merged_hapFIRE_allele_frequency-SNPlocations-",
#                                myatname,"-",commonname,".csv", sep="")
# cat(sed_command)
# cat(sed_command_locations)
# 
# # Extract vcf
# vcftools_command_locations <- paste0("vcftools --vcf plink.vcf --chr " ,chrom, " --from-bp ",start_bp," --to-bp ",end_bp ,
#                                      " --recode --recode-INFO-all --out subset_" , commonname )
# cat(vcftools_command_locations)
# 
# 
# 
# # Read VCF and visualize?
# library(GenotypePlot)
# 
# # Read VCF
# vcf_file<-paste0("data-external/arabidopsistranscriptomes/subset_",commonname,".recode.vcf")
# vcf<-my_vcf <- read.vcfR(vcf_file)
# vcf_subset <- vcf
# vcf_subset@gt <- vcf@gt[, c(1,which(colnames(vcf@gt) %in% individuals_to_keep))  ] # the 1 is to keep the first column format
# 
# 
# # gene model
# source("functions-gene-model-plot.R")
# genmod<-gene_model_plot(myatname)
# 
# # Sorted samples based on bio1
# our_popmap<-data.frame( ind=colnames(vcf_subset@gt)[-1],
#                         pop= NA) %>% 
#   # mutate(ecotypes = ind)%>%
#   separate(col = ind, into = c("id2","ecotypeid"),sep = "_") %>% 
#   merge(.,worldclim_ecotypesdata, by.x="ecotypeid",by.y="ecotypeid", all.x=T) %>% 
#   dplyr::mutate(pop=round(bio18_new/25)*25) %>% # dryness
#   dplyr::mutate(ind = paste0(id2,"_",id2)) %>% 
#   dplyr::select(ind, pop)
# our_popmap <- our_popmap[order(our_popmap$pop, our_popmap$ind), ]
# head(our_popmap)
# length(our_popmap$ind)
# our_popmap$pop %>% unique %>% length
# 
# # Plot VCF
# new_plot2 <- genotype_plot(vcf_object  =  vcf_subset,
#                            popmap = our_popmap,
#                            snp_label_size = 10000,
#                            plot_allele_frequency=F,
#                            cluster        = FALSE,
#                            colour_scheme=c("#d4b9da","#e7298a","#980043")
# )
# 
# 
# # Make a grid
# vcfplot<-
#   plot_grid(
#     genmod,
#     new_plot2$positions,
#     new_plot2$genotypes,
#     ncol=1,align = "v",
#     rel_heights = c(1,1,5)
#   )
# 
# vcfplot
# # 
# # save_plot(
# #   paste0("figs/fig-TSF-vcf-genemodel.pdf"),
# #   vcfplot,
# #   base_width = 7,base_height = 9 )
# # 
