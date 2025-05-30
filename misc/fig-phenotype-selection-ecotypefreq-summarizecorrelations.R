#### Summarize phenotype evolution in paper

#####*********************************************************************######

library(tidyverse)
library(dplyr)

library(ggplot2)
library(RColorBrewer)
library(tidyr)


#### Location
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)


#####*********************************************************************######
#### Read files #####

mmm<-df<-
    read_tsv(file = "tables/table_selection_coefficients.tsv")
head(mmm)


#####*#### Example locations correlations ####

df %>%
  dplyr::filter(r_p_dormancy <0.01) %>%
  dplyr::filter(year ==3 ) %>%
  ungroup %>%
  dplyr::select(sites_bio1, sites_bio18,city ,r_p_dormancy, r_dormancy) %>%
  arrange(r_dormancy)

df %>%
  dplyr::filter(r_p_ft10 <0.01) %>%
  dplyr::filter(year ==3 ) %>%
  # dplyr::filter(year ==1 ) %>%
  ungroup %>%
  dplyr::select(sites_bio1,city ,r_p_ft10, r_ft10) %>%
  arrange(r_ft10)

df %>%
  dplyr::filter(r_p_rosettearea <0.01) %>%
  dplyr::filter(year ==3 ) %>%
  # dplyr::filter(year ==1 ) %>%
  ungroup %>%
  dplyr::select(sites_bio1,city ,r_p_ft10, r_ft10,r_rosettearea,r_p_rosettearea) %>%
  arrange(r_rosettearea)

df %>%
  # dplyr::filter(r_p_ft10 < exp(-100)) %>%
  # dplyr::arrange(r_p_ft10) %>%
  # dplyr::filter(r_p_rootlength<0.01) %>%
  dplyr::filter(sites_bio18 <100) %>%
  dplyr::filter(any(r_p_c13<0.05, r_p_stomata<0.05,r_p_roothorizontality<0.05 ,r_p_rootlength<0.05)) %>%
  dplyr::filter(year ==3 ) %>%
  data.frame %>%
  dplyr::select(site, city,
                # starts_with("r_"),
                r_rootlength, r_la, r_c13, r_ft10
                ) %>%
  arrange(r_rootlength, r_la, r_c13)


df %>%
  dplyr::filter(year ==3 ) %>%
  # dplyr::filter(sites_bio1 <10 ) %>%
  data.frame %>%
  dplyr::select(city, starts_with("r_")) %>%
  # dplyr::select(site, city,
  #               # starts_with("r_"),
  #               r_rootlength, r_la, r_c13, r_ft10
  # ) %>%
  dplyr::select(city, contains("stoma")) %>%
  dplyr::filter(r_p_stomata<0.05) %>%
  arrange(r_stomata_density) %>%
  ht()



df %>%
  dplyr::filter(year ==3 ) %>%
  data.frame %>%
  dplyr::select(city, contains("c13"), starts_with("sites_bio")) %>%
  # dplyr::filter(r_p_stomata<0.05) %>%
  arrange(r_c13) %>%
  ht()


df %>%
  dplyr::filter(year ==3 ) %>%
  dplyr::filter(sites_bio1 <10 ) %>%
  data.frame %>%
  dplyr::select(city, starts_with("r_")) %>%
  # dplyr::select(site, city,
  #               # starts_with("r_"),
  #               r_rootlength, r_la, r_c13, r_ft10
  # ) %>%
  arrange(r_c13)


df %>%
  # dplyr::filter(r_p_ft10 < exp(-100)) %>%
  # dplyr::arrange(r_p_ft10) %>%
  dplyr::filter(r_p_rootlength<0.01) %>%
  # dplyr::filter(sites_bio18 <100) %>%
  # dplyr::filter(r_c13>0) %>%
  # dplyr::filter(any(r_p_c13<0.05)) %>%
  dplyr::filter(year ==3 ) %>%
  data.frame %>%
  dplyr::select(site, city,
                # starts_with("r_"),
                r_rootlength, r_la, r_c13, r_ft10
  ) %>%
  arrange(r_rootlength)

df %>%
  # dplyr::filter(r_p_ft10 < exp(-100)) %>%
  # dplyr::arrange(r_p_ft10) %>%
  # dplyr::filter(r_p_rootlength<0.01) %>%
  # dplyr::filter(sites_bio18 <100) %>%
  dplyr::filter(r_c13>0) %>%
  # dplyr::filter(any(r_p_c13<0.05)) %>%
  dplyr::filter(year ==3 ) %>%
  data.frame %>%
  dplyr::select(site, city,
                # starts_with("r_"),
                r_rootlength, r_la, r_c13, r_ft10
  ) %>%
  arrange(r_rootlength, r_la, r_c13)


df %>%
  # dplyr::filter(year ==3 ) %>%
  data.frame %>%
  dplyr::select(site, city,
                # starts_with("r_"),
                r_rootlength, r_la, r_c13, r_ft10,r_stomata_density
  ) %>%
  arrange(r_rootlength, r_la, r_c13,r_stomata_density)


df %>%
  dplyr::filter(year ==3 ) %>%
  dplyr::filter(city=="Tartu")
df %>%
  dplyr::filter(year ==3 ) %>%
  dplyr::filter(city=="Warsaw")


r_stomata_density


# sink("tables/test-phenotype-evolution-by-climate.txt")
cor.test(df$r_dormancy[df$year==3], df$sites_bio18[df$year==3])
cor.test(df$r_ft10[df$year==3], df$sites_bio1[df$year==3])
cor.test(df$r_la[df$year==3], df$sites_bio1[df$year==3])
cor.test(df$r_rosettearea[df$year==3], df$sites_bio18[df$year==3])
cor.test(df$r_rosettearea[df$year==3], df$sites_bio1[df$year==3])
cor.test(df$r_roothorizontality[df$year==3], df$sites_bio1[df$year==3])
cor.test(df$r_stomata_density[df$year==3], df$sites_bio1[df$year==3])
cor.test(df$r_stomata_density[df$year==3], df$sites_bio18[df$year==3])

cor.test(df$r_rootlength[df$year==3], df$sites_bio16[df$year==3])
cor.test(df$r_c13[df$year==3], df$sites_bio1[df$year==3])
# sink()


df %>%
  dplyr::filter(year ==3 ) %>%
  dplyr::select(r_ft10, starts_with("sites_bio")) %>%
  lm(data=.,  r_ft10 ~ . ) %>% summary

df %>%
  dplyr::filter(year ==3 ) %>%
  dplyr::select(r_rootlength, starts_with("sites_bio1")) %>%
  lm(data=.,  r_rootlength ~ . ) %>% summary

cor.test(df$r_rootlength[df$year==3], df$sites_bio16[df$year==3])
cor.test(df$r_rootlength[df$year==3], df$sites_bio18[df$year==3])
cor.test(df$r_rootlength[df$year==3], df$sites_bio15[df$year==3])


df %>%
  dplyr::filter(year ==1 ) %>%
  # dplyr::select(r_c13, starts_with("sites_bio")) %>%
  dplyr::select(r_c13, starts_with("sites_bio1")) %>%
  lm(data=.,  r_c13 ~ . ) %>% summary


# df %>%
#   dplyr::select(r_rootlength, starts_with("sites_bio")) %>%
#   dplyr::select(r_rootlength, sites_bio1,
#                 sites_bio2,
#                 sites_bio3,
#                 sites_bio4,
#                 sites_bio5,
#                 sites_bio6,
#                 sites_bio7,
#                 sites_bio8,
#                 sites_bio9,
#                 sites_bio10,
#                 sites_bio11
#                 ) %>%
#   lm(data=.,  r_rootlength ~ . ) %>% summary


df %>%
  dplyr::filter(year ==3 ) %>%
  dplyr::select(r_ft10, sites_bio1) %>%
  # dplyr::mutate(r_ft10 =r_ft10/mean(r_ft10) ) %>%
  lm(data=.,  r_ft10 ~ . ) %>% summary

df %>%
  dplyr::filter(year ==3 ) %>%
  dplyr::select(r_dormancy, sites_bio1) %>%
  # dplyr::mutate(r_ft10 =r_ft10/mean(r_ft10) ) %>%
  lm(data=.,  r_dormancy ~ . ) %>% summary



df %>%
  # dplyr::filter(year ==3 ) %>%
  dplyr::select(r_dormancy, starts_with("sites_bio")) %>%
  # dplyr::mutate(r_ft10 =r_ft10/mean(r_ft10) ) %>%
  lm(data=.,  r_dormancy ~ . ) %>% summary
