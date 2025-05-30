#### Quantify phenotypic evolutionary response from frequency and phenotype
#### correlations. Create supplemental table

#####*********************************************************************######

library(tidyverse)
library(dplyr)

library(ggplot2)
library(RColorBrewer)
library(tidyr)


# %%R
#### Location
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

#####*********************************************************************######
##### Load datasets #####

# Load the long format with ecotype frequencies and climates
# load("data-intermediate/ecotype_terminal_frequencies_long_raw.rda")
# load("data-intermediate/ecotype_terminal_frequencies_long_raw_climate_popsize.rda")
# load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize.rda")
load("data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize_flowers.rda")
load("data-external//atlas1001_phenotype_matrix_imputed_onlypheno_NEW.rda")

# Get the ecotype frequency long, rename
# e = ecotype_terminal_frequencies_long_raw
# e = ecotype_terminal_frequencies_long_raw_climate_popsize %>%
#   dplyr::rename(freq=maxfreq)
# e = ecotype_allyears_frequencies_long_raw_climate_popsize
e = ecotype_allyears_frequencies_long_raw_climate_popsize_flowers
head(e)

##### Merge phenotypes #####
# Combine with phenotype data, subset first
phenosub <-
  pheno %>%
  dplyr::select(id,
                Delta_13C,
                Leaf_Area,
                leafsize,
                rosette_DM,
                stomatasize,
                stomata_density,
                FT10,
                FT16,
                RGR,
                DSDS50,
                Root_horizontal_index_avg,
                Root_growth_rate_avg,
                Relative_root_length
  )
# Root names colnames(pheno)[grep("root", colnames(pheno))]

# merge frequency and phenotypes and group to correlate per site freq and phenotype
freqpheno<-mmmtmp<-merge(e,phenosub, by="id")

# Merge e and environment
mmerge<-merge(e,phenosub, by="id") %>%
  dplyr::select(-starts_with("sites_vap")) %>%
  dplyr::select(-starts_with("ecotypes_vap"))

#####*********************************************************************######
##### GROUP BY SITE #####
# Do the correlation between frequency and phenotypes
mmm<-mmerge %>%
  # group correlations by site
  group_by(
    across(c(site, starts_with("sites_"))),
    # rep, # if you want one per replicate
    year, # if you want one per year
    latitude, longitude, altitude) %>%
  # Scale variables,
  dplyr::mutate(
    FT16=FT16/sd(FT16,na.rm=T),
    FT10=FT10/sd(FT10,na.rm=T),
    Delta_13C=Delta_13C/sd(Delta_13C,na.rm=T),
    RGR=RGR/sd(RGR,na.rm=T),
    rosette_DM=rosette_DM/sd(rosette_DM,na.rm=T),
    DSDS50=DSDS50/sd(DSDS50,na.rm=T),
    Root_growth_rate_avg=Root_growth_rate_avg/sd(Root_growth_rate_avg,na.rm=T),
    Root_horizontal_index_avg=Root_horizontal_index_avg/sd(Root_horizontal_index_avg,na.rm=T),
    stomatasize=stomatasize/sd(stomatasize,na.rm=T),
    stomata_density=stomata_density/sd(stomata_density,na.rm=T)
  ) %>%

  # do the correlations
  dplyr::summarise(
    # Correlation
    r_c13= cor(Delta_13C,freq),
    n_c13= length(freq),
    r_p_c13= cor.test(Delta_13C,freq)$p.value,
    r_ft10= cor(FT10,freq),
    n_ft10= length(freq),
    r_p_ft10= cor.test(FT10,freq)$p.value,
    r_ft16= cor(FT16,freq),
    r_p_ft16= cor.test(FT16,freq)$p.value,
    r_la= cor(Leaf_Area,freq),
    r_p_la= cor.test(Leaf_Area,freq)$p.value,
    r_rosettearea= cor(rosette_DM,freq),
    r_p_rosettearea= cor.test(rosette_DM,freq)$p.value,
    r_p_roothorizontality= cor.test(Root_horizontal_index_avg,freq)$p.value,
    r_roothorizontality= cor(Root_horizontal_index_avg,freq),
    r_rootlength= cor(Relative_root_length,freq),
    r_p_rootlength= cor.test(Relative_root_length,freq)$p.value,
    r_rootgrowth= cor(Root_growth_rate_avg,freq),
    r_dormancy= cor(DSDS50,freq),
    r_p_dormancy= cor.test(DSDS50,freq)$p.value,
    r_stomata_density=cor(stomata_density,freq),
    r_stomata_size=cor(stomatasize,freq),
    r_p_stomata= cor.test(stomata_density,freq)$p.value,
    # Gradients
    # ft
    b_FT= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[2,1],
    p_FT= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[2,4],
    # c13
    b_c13= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[3,1],
    p_c13= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[3,4],
    # rgr
    b_c13= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[4,1],
    p_c13= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[4,4],
    # rosette
    b_rosette= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[3,1],
    p_rosette= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[3,4],
    # dormancy
    b_c13= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[4,1],
    p_c13= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[4,4],
    # root growth
    b_rootgrowth= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[5,1],
    p_rootgrowth= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[5,4],
    # rootangle
    b_roothorizontality= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[6,1],
    p_roothorizontality= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[6,4],
    # stomatasize
    b_stomatasize= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[7,1],
    p_stomatasize= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[7,4],
    # stomata_density
    b_stomata_density= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[8,1],
    p_stomata_density= coefficients(
      summary(lm(freq ~ FT16 + Delta_13C + RGR + rosette_DM + DSDS50 +Root_growth_rate_avg + Root_horizontal_index_avg + stomatasize + stomata_density)) # correcting for many things
    )[8,4],
  )

# Save table
sites_simple_names<-read.csv("grene/data/sites_simple_names.csv")

table_selection_coefficients<-
  df<-
  mmm %>%
  ungroup %>%
  dplyr::select( -starts_with("sites_"), starts_with("sites_bio")) %>%
  merge(.,sites_simple_names, by="site") %>%
  dplyr::select(names(sites_simple_names), everything()) # Reorder columns

df

write_tsv(file = "tables/table_selection_coefficients.tsv",
          x=table_selection_coefficients
)
