################################################################################
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
library(topr)
################################################################################

path_lfmm = 'data-intermediate/gea/lfmm_full/results/lfmm_bio1_4topr.csv'
lfmm <- read.csv(path_lfmm)
head(lfmm)
lfmm$SNP<-paste0(lfmm$CHROM,"_",lfmm$POS)
lfmm$ID<-lfmm$SNP

lfmmraw <- read.csv("data-intermediate/gea/lfmm_full/results/lfmm_bio1_results.csv")
lfmmraw<-
  lfmmraw %>% 
    mutate(snp=snp_id) %>% 
    rename(P=p_value) %>% 
    separate(snp_id, c("CHROM","POS")) %>% 
    select(P,CHROM,POS, snp)

path_annotations = 'data-intermediate/TAIR10_GFF3_genes_transposons_formatted4topr.csv'
annotations <- read.csv(path_annotations)

path_gemma='data-intermediate/gea/climate_gwas/lmm_gemma/bio1/output/result_4topr.csv'
gemma <- read.csv(path_gemma)


load("data-intermediate/alleles-TAIR10_parsedgenes.rda")

################################################################################
# Manhattan

lead_snps <- get_lead_snps(lfmm, thresh = 1e-8, region_size = 10000000)
lead_snps <-lead_snps %>% 
  merge(., annotated_snps,
        by.x="ID", by.y="snp",all.x=T) %>% 
  dplyr::select(P,CHROM,POS, ID, genic, closest_gene)

# manhattan(lfmmraw)
# manhattan(lfmm)

pdf("figs/fig-climate_gwas_lfmm_bio1.pdf",width = 10,height = 3)
# lfmmsubset <- lfmm %>% dplyr::filter(P<1e-5)
manhattan(list(lfmm, lead_snps), 
          build = annotations,
          sign_thresh = 1,
          size=c(1,3),
          color=c("darkgrey","#41AB5D"), 
          legend_labels=c("GrENE-net GEA", "Lead snp's"), verbose=F) 
dev.off()

png("figs/fig-climate_gwas_lfmm_bio1.png",width = 12,height = 3, units="in",res=1000)
# lfmmsubset <- lfmm %>% dplyr::filter(P<1e-5)
manhattan(list(lfmm, lead_snps), 
          build = annotations,
          sign_thresh = 1,
          size=c(1,3),
          color=c("darkgrey","#41AB5D"), 
          legend_labels=c("GrENE-net GEA", "Lead snp's"), verbose=F) 
dev.off()

cat(lead_snps$closest_gene, collapse="\t")
################################################################################

## plot lfmm and gemma results together
## start and end of the intersting block 2_973 2_9718274 2_9888341
## make it smaller 9750000
# climate_gwas_lfmm_block_2_973.pdf

### Block middle zoom out
pdf("figs/fig-climate_gwas_lfmm_block_2_973.pdf")
regionplot(list(gemma, lfmm), region="2:9750000-9888341", build = annotations,
           annotate_with_vline = 1e-8, color=c("darkgrey","#41AB5D"),
           legend_labels = c("1001G GEA", "GrENE-net GEA"), show_gene_legend = F,
           gene_color="#646464")
dev.off()
png("figs/fig-climate_gwas_lfmm_block_2_973.png",width = 9,height = 9,res = 600, units = "in")
regionplot(list(gemma, lfmm), region="2:9750000-9888341", build = annotations,
           annotate_with_vline = 1e-8, color=c("darkgrey","#41AB5D"),
           legend_labels = c("1001G GEA", "GrENE-net GEA"), show_gene_legend = F,
           gene_color="#646464")
dev.off()

### Block middle zoom in
# climate_gwas_lfmm_block_2_973.pdf ZOOM
pdf("figs/fig-climate_gwas_lfmm_block_2_973-zoom.pdf")
regionplot(list(gemma, lfmm), region="2:9800000-9820000", build = annotations,
           annotate_with_vline = 1e-8, color=c("darkgrey","#41AB5D"),
           legend_labels = c("1001G GEA", "GrENE-net GEA"), show_gene_legend = F,
           gene_color="#646464")
dev.off()

png("figs/fig-climate_gwas_lfmm_block_2_973-zoom-2_9809464.png",width = 9,height = 9,res = 600, units = "in")
regionplot(list(gemma, lfmm), region="2:9800000-9820000", build = annotations,
           annotate_with_vline = 1e-8, color=c("darkgrey","#41AB5D"),
           legend_labels = c("1001G GEA", "GrENE-net GEA"), show_gene_legend = F,
           gene_color="#646464")
dev.off()

png("figs/fig-climate_gwas_lfmm_block_2_973-zoom-2_9882323.png",width = 9,height = 9,res = 600, units = "in")
buffer=5000
chr=2
pos=9882323
myregion<-paste0(chr,":",pos-buffer , "-", pos+buffer)
regionplot(list(gemma, lfmm), region=myregion, build = annotations,
           annotate_with_vline = 1e-8, color=c("darkgrey","#41AB5D"),
           legend_labels = c("1001G GEA", "GrENE-net GEA"), show_gene_legend = F,
           gene_color="#646464")
dev.off()


# Region wuth the highest peak
buffer=5000
chr=2
pos=9902238
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


# CAM5
buffer=5000
chr=2
pos=9902238
myregion<-paste0(chr,":",pos-buffer , "-", pos+buffer)
regionplot(list(gemma, lfmm), region=myregion, build = annotations,
           annotate_with_vline = 1e-8, color=c("#41AB5D"),
           legend_labels = c("1001G GEA", "GrENE-net GEA"), show_gene_legend = F,
           gene_color="#646464")