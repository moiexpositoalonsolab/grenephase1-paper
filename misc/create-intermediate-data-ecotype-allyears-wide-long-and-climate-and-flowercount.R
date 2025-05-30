################################################################################
### Goal
### Make several convenient intermediate files for faster plot production later

################################################################################
library(tidyverse)
library(dplyr)

################################################################################

# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives/MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses" # GOOGLE DRVIE


# Specify the file path
file_edelta = '/data/delta_ecotype_freq.txt'
file_efreq = '/data/merged_ecotype_frequency.txt'
file_e0 = '/data/founder_ecotype_frequency.txt'
file_eid = '/data/greneNet_final_v1.1.fam'
file_climate_sites="/grene/data/worldclim_sitesdata.rda"
file_sites_info="/grene/data/sites.clim.rda"

file_path_edelta <- paste0(myfolder, file_edelta)
file_path_efreq <- paste0(myfolder, file_efreq)
file_path_e0 <- paste0(myfolder, file_e0)
file_path_eid <- paste0(myfolder, file_eid)
file_path_climate_sites=paste0(myfolder, file_climate_sites)
file_path_sites_info=paste0(myfolder, file_sites_info)


################################################################################
# Load ecotypes and their ID
eco0<-read.table(file_path_e0, header=F) %>% 
  dplyr::rename(id=V1,e0=V2)
ecofreq<-read.table(file_path_efreq, header=T) 
ecodelta<-read.table(file_path_edelta, header=T)


ecofreq<-data.frame(ecofreq)
ecodelta<-data.frame(ecodelta)
fam<-read.table(file_path_eid)
fam<-data.frame(fam)
# eco<-apply(ecodelta,2,function(x) x+eco0)  %>% do.call(cbind, .) # comment out

################################################################################
# Make long to parse the column names
ecodelta$id <- fam$V1
ecolong<-
  ecodelta %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to=c("site","year","rep"),
    values_to = c("freq"),
    names_pattern = "X?(.*)_(.)_(.*)"
  ) %>%
  dplyr::select(site, year, rep, freq, id)

df<- ecolong
df_delta<-df

################################################################################
ecotype_allyears_frequencies_long = df
write.csv( file= paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long.csv"), 
           ecotype_allyears_frequencies_long)
save(file=paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long.rda"),
     ecotype_allyears_frequencies_long)
# ################################################################################
# # Spread the data to wide format
# df_wide <-
#   df %>%
#   mutate(site_rep=paste0(site,"_",rep,"_",)) %>% 
#   pivot_wider(names_from = id, values_from = freq, id_cols = site_rep)
# head(df_wide)
# 
# ################################################################################
# ecotype_allyears_frequencies_wide = df_wide 
# write.csv( file= paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_wide.csv"),
#            ecotype_allyears_frequencies_wide)
# save(file=paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_wide.rda"),
#      ecotype_allyears_frequencies_wide)

################################################################################
################################################################################
# Now the raw frequencies or ratio
################################################################################
################################################################################

# Get the raw frequencies and parse

ecofreq$id <- fam$V1
# Make long to parse the column names
ecolong<-
  ecofreq %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to=c("site","year","rep"),
    values_to = c("freq"),
    names_pattern = "X?(.*)_(.)_(.*)"
  ) %>%
  dplyr::select(site, year, rep, freq, id)

df<- ecolong

# Add the starting frequencies
eco0<-eco0 %>% apply(.,2,as.numeric) %>% data.frame 
df <- df %>% 
  merge(.,eco0,by='id') %>% 
  dplyr::rename(startfreq=e0)

df_raw<-df


message("saving ecotype all years frequency, long format, raw")
ecotype_allyears_frequencies_long_raw = df 
write.csv( file= paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long_raw.csv"), 
           ecotype_allyears_frequencies_long_raw)
save(file=paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long_raw.rda"),
     ecotype_allyears_frequencies_long_raw)


################################################################################
# Add climate
###############
# Load climates
file_path_climate_sites=paste0(myfolder, "/grene/data/worldclim_sitesdata.rda")
file_path_sites_info=paste0(myfolder, "/grene/data/locations_data.rda")
file_path_climate_ecotypes=paste0(myfolder, "/grene/data/worldclim_ecotypesdata.rda")

load(file_path_climate_sites)
load(file_path_sites_info)
load(file_path_climate_ecotypes)

colnames(worldclim_sitesdata)[-1]<-paste0("sites_",colnames(worldclim_sitesdata)[-1])
colnames(worldclim_ecotypesdata)[-1] <- paste0("ecotypes_",colnames(worldclim_ecotypesdata)[-1])
locations_data<-locations_data %>% 
  dplyr::select(site,longitude,latitude)

### Merge with delta
e<-df_delta
e<-
  e %>% 
  merge(
    .,
    locations_data, by='site'
  ) %>% 
  merge(
    .,
    worldclim_sitesdata, by="site"
  ) %>% 
  merge(
    .,
    worldclim_ecotypesdata, by.x="id", by.y="ecotypeid"
  )

## Save
message("saving ecotype all years frequency, long format, delta, with climate")
ecotype_allyears_frequencies_long_delta_climate = e
e[1:5,1:5]
write.csv( file= paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long_delta_climate.csv"), 
           ecotype_allyears_frequencies_long_delta_climate)
save(file=paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long_delta_climate.rda"),
     ecotype_allyears_frequencies_long_delta_climate)
################################################################################

### Merge with raw
e<-df_raw
e<-
  e %>% 
  merge(
    .,
    locations_data, by='site'
  ) %>% 
  merge(
    .,
    worldclim_sitesdata, by="site"
  ) %>% 
  merge(
    .,
    worldclim_ecotypesdata, by.x="id", by.y="ecotypeid"
  )

## Save
message("saving ecotype all years frequency, long format, raw, with climate")
ecotype_allyears_frequencies_long_raw_climate = e
e[1:5,1:5]
write.csv( file= paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long_raw_climate.csv"), 
           ecotype_allyears_frequencies_long_raw_climate)
save(file=paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long_raw_climate.rda"),
     ecotype_allyears_frequencies_long_raw_climate)

################################################################################
# Add population size
###############

# This is the population size created by Tati
popsize=read.csv(paste0(myfolder, "/data-intermediate//pop_size_estimations.csv"))
popsize<-data.frame(popsize)
popsize<-popsize %>% dplyr::select(-X,-longitude,-latitude, -predicted)
popsize[1:5,1:8]

# Merge with ecotype frequency
e<-ecotype_allyears_frequencies_long_raw_climate
e[1:5,1:5]
epop<-
  merge(
  x=e, by.x=c("site","rep","year"),
  y=popsize, by.y=c("site","plot","generation")
)

head(epop)
epop[1:5,1:10]
colnames(epop)
dim(epop)

## Save
message("saving ecotype all years frequency, long format, raw, with climate, with population size")
ecotype_allyears_frequencies_long_raw_climate_popsize = epop
epop[1:5,1:10]
ecotype_allyears_frequencies_long_raw_climate_popsize[1:5,1:10]
write.csv( file= paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize.csv"), 
           ecotype_allyears_frequencies_long_raw_climate_popsize)
save(file=paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize.rda"),
     ecotype_allyears_frequencies_long_raw_climate_popsize)


################################################################################
# Add actual flower counts
###############

load("grene/data/samples_data.rda")
samples_tmp<-samples_data %>% 
  dplyr::filter(year<2021) %>% 
  dplyr::mutate(year= year-2017) %>% 
  dplyr::rename(flowers=flowerscollected) %>% 
  dplyr::select(site,plot,year, date,month,day,flowers) %>% 
  dplyr::mutate(julian_day = yday(ymd(date)))

# Reduce the size to only bio1-19
dat<-ecotype_allyears_frequencies_long_raw_climate_popsize

dat<-
  dat %>% 
  merge(
    .,
    samples_tmp, 
    by.x=c("site","year","rep"), by.y=c("site","year","plot")
  )

ecotype_allyears_frequencies_long_raw_climate_popsize_flowers <- dat
write.csv( file= paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize_flowers.csv"), 
           ecotype_allyears_frequencies_long_raw_climate_popsize_flowers)
save(file=paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_long_raw_climate_popsize_flowers.rda"),
     ecotype_allyears_frequencies_long_raw_climate_popsize_flowers)

