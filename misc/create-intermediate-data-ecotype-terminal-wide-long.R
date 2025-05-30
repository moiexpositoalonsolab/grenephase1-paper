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


# Extract the terminal year of frequency change
df2<-
  df %>%
  group_by(site, rep, id) %>%
  mutate(year=as.numeric(year)) %>%
  reframe(max_year = if(all(is.na(freq))) NA else max(year[!is.na(freq)], na.rm = TRUE),
          maxfreq = freq[year==max_year])

################################################################################
ecotype_terminal_frequencies_long = df2
write.csv( file= paste0(myfolder,"/data-intermediate/ecotype_terminal_frequencies_long_delta.csv"), 
           ecotype_terminal_frequencies_long)
save(file=paste0(myfolder,"/data-intermediate/ecotype_terminal_frequencies_long_delta.rda"),
            ecotype_terminal_frequencies_long)
################################################################################
# Spread the data to wide format
df_wide <-
  df2 %>%
  mutate(site_rep=paste0(site,"_",rep)) %>% 
  pivot_wider(names_from = id, values_from = maxfreq, id_cols = site_rep)
head(df_wide)

################################################################################
ecotype_terminal_frequencies_wide = df_wide 
write.csv( file= paste0(myfolder,"/data-intermediate/ecotype_terminal_frequencies_wide_delta.csv"),
           ecotype_terminal_frequencies_wide)
save(file=paste0(myfolder,"/data-intermediate/ecotype_terminal_frequencies_wide.rda_delta"),
     ecotype_terminal_frequencies_wide)

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

# Extract the terminal year of frequency change
df3<-
  df %>%
  group_by(site, rep, id) %>%
  mutate(year=as.numeric(year)) %>%
  reframe(max_year = if(all(is.na(freq))) NA else max(year[!is.na(freq)], na.rm = TRUE),
          maxfreq = freq[year==max_year])

head(df3)
# Add the starting frequencies
eco0<-eco0 %>% apply(.,2,as.numeric) %>% data.frame 
df3 <- df3 %>% 
  merge(.,eco0,by='id') %>% 
  dplyr::rename(startfreq=e0)

################################################################################
ecotype_terminal_frequencies_long_raw = df3 
write.csv( file= paste0(myfolder,"/data-intermediate/ecotype_terminal_frequencies_long_raw.csv"), 
           ecotype_terminal_frequencies_long_raw)
save(file=paste0(myfolder,"/data-intermediate/ecotype_terminal_frequencies_long_raw.rda"),
     ecotype_terminal_frequencies_long_raw)

################################################################################
# Now wide
df_wide <-
  df3 %>%
  dplyr::select(-startfreq) %>% 
  mutate(site_rep=paste0(site,"_",rep)) %>% 
  pivot_wider(names_from = id, values_from = maxfreq, id_cols = site_rep)
head(df_wide)

################################################################################
ecotype_terminal_frequencies_wide_raw = df_wide 
write.csv( file= paste0(myfolder,"/data-intermediate/ecotype_allyears_frequencies_wide_raw.csv"), 
           ecotype_terminal_frequencies_wide_raw)
save(file=paste0(myfolder,"/data-intermediate/ecotype_terminal_frequencies_wide_raw.rda"),
     ecotype_terminal_frequencies_wide_raw)
