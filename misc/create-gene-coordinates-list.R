################################################################################
### Goal
### Create a easily readable map of SNPs to gene names

################################################################################
#### If transformed to collab can use the stuff below
# from google.colab import drive
# drive.mount('/content/drive')
# %load_ext rpy2.ipython
# %%R

################################################################################
#### Libraries
library(tidyverse)

################################################################################
#### Location if run locally
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
# myfolder<-"~/grenephase1-analyses"

myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
#### Get the TAIR genes

tair<-read.table(file = "data-external/TAIR10_GFF3_genes_transposons.gff")

# Rename and select 
tair<-
  tair %>% 
  dplyr::filter(V3=="gene") %>% 
  dplyr::select(-V6,-V8)

# Extract just the gene name
tair$gene<- sapply( tair$V9, function(i){ gsub(tail(unlist(strsplit(i, ";",fixed = T) ) ,1 ), pattern="Name=",replacement ="" )  })

# Final filtering and renaming
tair<-
  tair %>% dplyr::select(-V9, -V2, -V3) %>% 
  dplyr::rename(Chromosome=V1,Position_start=V4, Position_end=V5, Strand=V7,Gene=gene) %>% 
  dplyr::mutate(Chromosome=gsub(Chromosome,pattern="Chr", replace="",fixed = T))


write.csv(file = "data-external/TAIR10_parsedgenes.csv", tair,row.names = F)
tair<-read.csv("data-external/TAIR10_parsedgenes.csv")

################################################################################
#### Annotate snps for the LD prunned data

# Load the SNPs from the pruned 13K dataset
p0<-read.table("data/average_seedmix_p0_LDpruned.txt", header=F)
p0$snp<-paste0(p0$V1,"_",p0$V2)

# Rename to merge with tair
snps<-
  p0 %>% 
  dplyr::select(V1, V2, snp) %>% 
  dplyr::rename(Chr=V1,Pos=V2) #%>% 

################################################################################
## Annotate SNPs with genes


# Function 
annotate_snps <- function(snps, tair) {
  annotated_snps <- snps %>%
    rowwise() %>%
    mutate(
      genic = if_else(any(
          Chr == tair$Chromosome & 
          Pos >= tair$Position_start & 
          Pos <= tair$Position_end), 
        "genic", 
        "intergenic"),
      closest_gene = if (genic == "genic") {
        tair %>%
          filter(Chromosome == Chr & 
                   Position_start <= Pos  & 
                   Position_end >= Pos) %>%
          pull(Gene) %>%
          paste(collapse = ", ")
      } else {
        tair %>%
          filter(Chromosome == Chr) %>%
          arrange(abs(Pos - Position_start)) %>%
          slice(1) %>%
          pull(Gene)
      }
    ) %>%
    ungroup()
  
  return(annotated_snps)
}

annotated_snps <- annotate_snps(snps, tair)
print(annotated_snps)

write.csv(file = "data-external/alleles-LDpruned-TAIR10_parsedgenes.csv", annotated_snps)
save(file = "data-intermediate/alleles-LDpruned-TAIR10_parsedgenes.rda", annotated_snps)

 

################################################################################
#### Annotate snps for all snps

# Load the SNPs from the pruned 13K dataset
p0<-read.table("data/average_seedmix_p0.txt", header=F)
p0$snp<-paste0(p0$V1,"_",p0$V2)

# Rename to merge with tair
snps<-
  p0 %>% 
  dplyr::select(V1, V2, snp) %>% 
  dplyr::rename(Chr=V1,Pos=V2) #%>% 

################################################################################
## Annotate SNPs with genes


# Function 
annotate_snps <- function(snps, tair) {
  annotated_snps <- snps %>%
    rowwise() %>%
    mutate(
      genic = if_else(any(
        Chr == tair$Chromosome & 
          Pos >= tair$Position_start & 
          Pos <= tair$Position_end), 
        "genic", 
        "intergenic"),
      closest_gene = if (genic == "genic") {
        tair %>%
          filter(Chromosome == Chr & 
                   Position_start <= Pos  & 
                   Position_end >= Pos) %>%
          pull(Gene) %>%
          paste(collapse = ", ")
      } else {
        tair %>%
          filter(Chromosome == Chr) %>%
          arrange(abs(Pos - Position_start)) %>%
          slice(1) %>%
          pull(Gene)
      }
    ) %>%
    ungroup()
  
  return(annotated_snps)
}

annotated_snps <- annotate_snps(snps, tair)
print(annotated_snps)

write.csv(file = "data-external/alleles-TAIR10_parsedgenes.csv", annotated_snps)
save(file = "data-intermediate/alleles-TAIR10_parsedgenes.rda", annotated_snps)

################################################################################
# Install BiocManager if you haven't already
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install biomaRt from Bioconductor
BiocManager::install("biomaRt")

# Load the biomaRt package
library(biomaRt)

# Connect to the Ensembl Plants BioMart database
ensembl <- useMart("plants_mart", host = "https://plants.ensembl.org")
ensembl <- useDataset("athaliana_eg_gene", mart = ensembl)

# List of gene IDs
genes <- c("AT1G52050", "AT1G80500", "AT1G23060", "AT1G03710", "AT2G41380",
           "AT2G02160", "AT2G23270", "AT3G29185", "AT3G62800", "AT3G01860",
           "AT4G30860", "AT4G08910", "AT4G01440", "AT5G65207", "AT5G02030",
           "AT5G26010")

# Retrieve the annotations, including descriptions
annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = ensembl
)

# Print the annotations
print(annotations)
