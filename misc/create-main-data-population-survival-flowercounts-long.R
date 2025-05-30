################################################################################
### Goal
### Parse the survival.csv file from Tati
### Create a long format dataset with summary survival and flower counts
################################################################################

myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives/MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses" # GOOGLE DRVIE

setwd(myfolder)

################################################################################
library(tidyr)

################################################################################
survival<-read.csv("data/survival.csv")

survival<-
  survival %>% 
  dplyr::filter(plot %in% c(1:12)) %>% 
  dplyr::select(site, plot, contains("flower"), contains("surv")) %>% 
    pivot_longer(
      cols = starts_with("X"),
      names_to = c("year", "variable"),
      names_pattern = "X(\\d+)_(.*)"
    ) %>%
    pivot_wider(
      names_from = variable,
      values_from = value
    ) %>%
    mutate(year = as.integer(year))

flowers_survival_long = survival
write.csv( file= paste0(myfolder,"/data/flowers_survival_long.csv"), 
           flowers_survival_long)
save(file=paste0(myfolder,"/data/flowers_survival_long.rda"),
     flowers_survival_long)