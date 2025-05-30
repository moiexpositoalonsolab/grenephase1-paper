library(dplyr)
library(purrr)
library(tidyr)

myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses"
myfolder<-"~/grenephase1-analyses"

# Specify the file path
file_edelta = '/data/delta_ecotype_freq.txt'
file_e0 = '/data/merged_ecotype_frequency.txt'
file_eid = '/data/greneNet_final_v1.0.fam'
file_climate_sites="/grene/data/worldclim_sitesdata.rda"

file_path_edelta <- paste0(myfolder, file_edelta)
file_path_e0 <- paste0(myfolder, file_e0)
file_path_eid <- paste0(myfolder, file_eid)
file_path_climate_sites=paste0(myfolder, file_climate_sites)

file_plot_correlation="/fig-intrasite-correlation-ecotypes-bio1.pdf"
################################################################################
# library(grene)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
# library(cowplot)
# theme_set(theme_cowplot())
# theme_set(theme_classic)

# Load ecotypes and their ID
eco<-read.table(file_path_e0, header=T)
ecodelta<-read.table(file_path_edelta, header=T)
fam<-read.table(file_path_eid)
ecodelta$id<- fam$V1

# Check dataset
head(ecodelta[,1:5])

# head(fam)
# print(eco[,1:5])

# df<-eco
# Make long to parse the column names
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
df<-
df %>%
  group_by(site, rep, id) %>%
  mutate(year=as.numeric(year)) %>%
  reframe(max_year = if(all(is.na(freq))) NA else max(year[!is.na(freq)], na.rm = TRUE),
         maxfreq = freq[year==max_year]) 
         
# Spread the data to wide format
df_wide <- 
df %>%
  pivot_wider(names_from = rep, values_from = maxfreq, id_cols = c(site, id))

# # Apply the function to each site
# df_wide %>%
#   group_by(site) %>%
#   summarise(avg_cor = cor_fun(across(-site)))


# Function to calculate correlation for each pair of columns
cor_fun <- function(df) {
  cor_values <- combn(df[, -1], 2, function(x) {
    if(sum(!is.na(x[1])) < 2 || sum(!is.na(x[2])) < 2) {
      return(NA)
    } else {
      return(cor(x[1], x[2], use = "pairwise.complete.obs"))
    }
  }, simplify = FALSE)
  mean(unlist(cor_values), na.rm = TRUE)
}

# Split the data by site
df_split <- split(df_wide, df_wide$site)

# Apply the function to each subset
avg_cor <- sapply(df_split, cor_fun)

# Create a data frame to summarize the results
result <- data.frame(site = names(avg_cor), avg_cor = avg_cor)
result<-apply(result,2,as.numeric)

###############
# Load climates
load(file_path_climate_sites)
worldclim_sitesdata<-apply(worldclim_sitesdata,2,as.numeric)

resultmerge=merge(result,worldclim_sitesdata,by="site")

# Plot
corplot<-
ggplot(resultmerge, aes(x=bio1, y=avg_cor, size=avg_cor)) +
  geom_point()+
  stat_smooth(method = "glm")+
  scale_size_continuous(guide="none")+
  ylab("Evolutionary parallelism (r)")+
  xlab("Annual temperature (C)")+
  theme_minimal()
corplot
save_plot(corplot,
          filename=paste0(myfolder,file_plot_correlation),
          base_height = 4.5,base_width = 5)
