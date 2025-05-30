rm(list=ls())

library(dplyr)

## local adaptation

## This script generates the stabilizing selection data frame which has the meta information of all ecotypes,
## experimental sites, plot information, and biolcimate variable distance for the TERMINAL GENERATION

prefix <- "LOCAL_ADAPTATION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)



## load the meta data and frequency files
meta <- read.delim("data/merged_sample_table.csv",sep=",",stringsAsFactors = F)
ecotype_frequency <- read.delim("data/merged_ecotype_frequency.txt")
terminal_index <- read.delim("data/site_plot_terminal_index.txt")

## Due to technical or biological reasons, some site-plots dont have all 3 years data after quality filtering. 
##It is hard to conclude the whether it is due to a bad year or bad sampling, therefore, we only use the site-plot
## that have year 1, year 1&2, year 1&2&3 to build models to see how well the model predicts across years

## generation 1
table(terminal_index$generation1==1)
## generation 2
table(terminal_index$generation1==1 &terminal_index$generation2==2)
## generation 3
table(terminal_index$generation1==1 &terminal_index$generation2==2 & terminal_index$generation3==3)

passed_index <- which(terminal_index$generation1==1 | terminal_index$generation1==1 &terminal_index$generation2==2 | terminal_index$generation1==1 &terminal_index$generation2==2 & terminal_index$generation3==3)
terminal_index_passed <- terminal_index[passed_index,]


terminal_meta <- meta[terminal_index_passed$index3,]
terminal_ecotype_frequency <- ecotype_frequency[,terminal_index_passed$index3]
experimental_sites <- unique(terminal_meta$site)

## load the ecotype initial frequency 
ecotype_p0 <- read.delim("data/founder_ecotype_frequency.txt",head=F)

## define some global parameters
num_ecotype <- nrow(ecotype_p0)
num_pop <- nrow(terminal_meta)
num_experiment_site <- length(experimental_sites)
founder_ecotype <- ecotype_p0$V1
founder_ecotype_frequency <- ecotype_p0$V2


## read the site and ecotype climate data
site_climate <- read.delim("data-external/bioclimvars_experimental_sites_era5.csv",sep=",")
colnames(site_climate) <- c("site",paste0("bio",1:19,"_site"))

## add longitude and latitude info to the site climate
site_info <- read.delim("data/sites_clim.csv",sep=",")
colnames(site_info)[c(3,14,15)] <- c("site","longitude_site","latitude_site")


site_climate <- site_climate %>%
  left_join(site_info[,c(3,14,15)],by="site") %>%
  filter(site %in% experimental_sites)


## now use the era5 data for ecotype climate of origin
ecotype_climate <- read.delim("data-external/bioclimvars_ecotypes_era5.csv",sep=",")
colnames(ecotype_climate) <- c("ecotype",paste0("bio",1:19,"_ecotype"))
ecotype_info <- read.delim("data/1001g_regmap_grenet_ecotype_info_corrected_2024May16.csv",sep = ",")
colnames(ecotype_info)[c(1,5,6)] <- c("ecotype","latitude_ecotype","longitude_ecotype")

## filter ecotype_info down to GrENE.net ecotypes only
ecotype_info <- ecotype_info %>%
  filter(GrENE.net == TRUE) %>%
  arrange(match(ecotype,founder_ecotype))

ecotype_climate <- ecotype_climate %>%
  filter(ecotype %in% founder_ecotype) %>%
  arrange(match(ecotype,founder_ecotype)) %>%
  left_join(ecotype_info[,c(1,6,5)],by="ecotype")

table(ecotype_climate$longitude==ecotype_info$longitude)
table(ecotype_climate$latitude==ecotype_info$latitude)
table(ecotype_climate$ecotype==ecotype_info$ecotype)





## now lets generate the pt0 and pt1 for every generation (p0 stands for the ecotype frequency of the founder generation, and pt stands for the ecotype frequency of the terminal generation)
pt <- c()
p0 <- c()

for(i in 1:num_ecotype){
  pt <- c(pt,as.numeric(terminal_ecotype_frequency[i,]))
  p0 <- c(p0,rep(founder_ecotype_frequency[i],num_pop))
}

M <- c()
for(j in 1:num_ecotype){
    M <- c(M,rep(founder_ecotype[j],num_pop))
}

table(unique(M) == founder_ecotype)


## generate N (site)
N <- rep(terminal_meta$site,231)

## generate R (plot)
R <- rep(terminal_meta$plot,231)

## generate G (generation)
G <- rep(terminal_meta$generation,231)

length(R)
length(M)
length(N)
length(G)


mydata <- as.data.frame(matrix(nrow = length(pt),ncol = 6))
colnames(mydata) <- c("pt","p0","ecotype","site","plot","generation")
mydata$pt <- pt
mydata$p0 <- p0
mydata$ecotype <- M
mydata$site <- N
mydata$plot <- R
mydata$generation <- G

mydata <- mydata %>%
  left_join(ecotype_climate,by="ecotype") %>%
  left_join(site_climate,by="site")

write.table(mydata,"data-intermediate/stabilizing_selection_data_2024Jun18_collectionsite_era5_terminal_generation_pt_p0.txt",append = F,quote = F,sep = "\t",row.names = F)

