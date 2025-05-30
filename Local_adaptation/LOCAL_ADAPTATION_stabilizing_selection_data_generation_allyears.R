rm(list=ls())

library(dplyr)

## local adaptation

## This script generates the stabilizing selection data frame which has the meta information of all ecotypes,
## experimental sites, plot information, and biolcimate variable distance by all three generations

## To derive the  model log(p_t+1/p_t) = log(Wmax) - log(Wavg) - 1/V (Z-Z0)^2, we need to rely on the equation 
## p_t+1 = p_t * W / W_avg assuming the population has passed one generation and has infinite population size, p_t+1 and p_t are the ecotype 
## frequency of the current and previous generation

## following the logic, it is necessary to generate the dataframe per generation and update the p_t with true
## ecotype frequency from the previous generation

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

site_plot_index <- c(terminal_index$index1[terminal_index$generation1==1],
                     terminal_index$index2[terminal_index$generation1==1 &terminal_index$generation2==2],
                     terminal_index$index3[terminal_index$generation1==1 &terminal_index$generation2==2 & terminal_index$generation3==3])
experimental_sites <- unique(meta$site[site_plot_index])

experimental_meta <- meta[site_plot_index,]

## load the ecotype initial frequency 
ecotype_p0 <- read.delim("data/founder_ecotype_frequency.txt",head=F)

## define some global parameters
num_ecotype <- nrow(ecotype_p0)
num_pop <- nrow(experimental_meta)
num_experiment_site <- length(experimental_sites)
founder_ecotype <- ecotype_p0$V1
founder_ecotype_frequency <- ecotype_p0$V2


## read the site and ecotype climate data
## load the era5 site climate data
site_climate_2018 <- read.delim("data-external/bioclimvars_sites_era5_year_2018.csv",sep=",")
colnames(site_climate_2018) <- c("site",paste0(colnames(site_climate_2018)[2:20],"_site"))
site_climate_2019 <- read.delim("data-external/bioclimvars_sites_era5_year_2019.csv",sep=",")
colnames(site_climate_2019) <- c("site",paste0(colnames(site_climate_2019)[2:20],"_site"))
site_climate_2020 <- read.delim("data-external/bioclimvars_sites_era5_year_2020.csv",sep=",")
colnames(site_climate_2020) <- c("site",paste0(colnames(site_climate_2020)[2:20],"_site"))

## add longitude and latitude info to the site climate
site_info <- read.delim("data/sites_clim.csv",sep=",")
colnames(site_info)[c(3,14,15)] <- c("site","longitude_site","latitude_site")

site_climate_2018 <- site_climate_2018 %>%
  left_join(site_info[,c(3,14,15)],by="site") %>%
  filter(site %in% experimental_sites)

site_climate_2019 <- site_climate_2019 %>%
  left_join(site_info[,c(3,14,15)],by="site") %>%
  filter(site %in% experimental_sites)

site_climate_2020 <- site_climate_2020 %>%
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





## now lets generate the pt0 and pt1 for every generation (pt0 stands for the ecotype frequency of the previous generation, and pt1 stands
## for the ecotype frequency of the current generation)
pt1 <- c()
pt0 <- c()

## from generation 0 to generation 1
pt1_1 <- c()
pt0_1 <- c()
gen1_index <- terminal_index$index1[terminal_index$generation1==1]
gen1_ecotype_frequency <- ecotype_frequency[,gen1_index]

for(i in 1:num_ecotype){
  pt1_1 <- c(pt1_1,as.numeric(gen1_ecotype_frequency[i,]))
  pt0_1 <- c(pt0_1,rep(founder_ecotype_frequency[i],length(gen1_index)))
}
hist(log(pt1_1 / pt0_1),breaks = 100)

# ## from generation 0 to generation 1
# log_p1_p0 <- c()
# gen1_index <- terminal_index$index1[terminal_index$generation1==1]
# gen1_ecotype_frequency <- ecotype_frequency[,gen1_index]
# 
# for(i in 1:num_ecotype){
#   log_p1_p0 <- c(log_p1_p0,log(as.numeric(gen1_ecotype_frequency[i,]) / founder_ecotype_frequency[i]))
# }
# hist(log_p1_p0,breaks = 100)
# 
# table(log(pt1_1 / pt0_1) == log_p1_p0)


## from generation 1 to generation 2
pt1_2 <- c()
pt0_2 <- c()
gen2_index <- terminal_index$index2[terminal_index$generation1==1 & terminal_index$generation2==2 ]
gen1_index <- terminal_index$index1[terminal_index$generation1==1 & terminal_index$generation2==2 ]
gen2_ecotype_frequency <- ecotype_frequency[,gen2_index]
gen1_ecotype_frequency <- ecotype_frequency[,gen1_index]

for(i in 1:num_ecotype){
  pt1_2 <- c(pt1_2,as.numeric(gen2_ecotype_frequency[i,]))
  pt0_2 <- c(pt0_2,as.numeric(gen1_ecotype_frequency[i,]))
}
hist(log(pt1_2 / pt0_2),breaks = 100)


# ## from generation 1 to generation 2
# log_p2_p1 <- c()
# gen2_index <- terminal_index$index2[terminal_index$generation1==1 & terminal_index$generation2==2 ]
# gen1_index <- terminal_index$index1[terminal_index$generation1==1 & terminal_index$generation2==2 ]
# gen2_ecotype_frequency <- ecotype_frequency[,gen2_index]
# gen1_ecotype_frequency <- ecotype_frequency[,gen1_index]
# 
# for(i in 1:num_ecotype){
#   log_p2_p1 <- c(log_p2_p1,log(as.numeric(gen2_ecotype_frequency[i,]) / as.numeric(gen1_ecotype_frequency[i,])))
# }
# hist(log_p2_p1,breaks = 100)
# 
# table(log(pt1_2 / pt0_2) == log_p2_p1)


## from generation 2 to generation 3
pt1_3 <- c()
pt0_3 <- c()
gen3_index <- terminal_index$index3[terminal_index$generation1==1 & terminal_index$generation2==2 & terminal_index$generation3==3]
gen2_index <- terminal_index$index2[terminal_index$generation1==1 & terminal_index$generation2==2 & terminal_index$generation3==3]
gen3_ecotype_frequency <- ecotype_frequency[,gen3_index]
gen2_ecotype_frequency <- ecotype_frequency[,gen2_index]

for(i in 1:num_ecotype){
  pt1_3 <- c(pt1_3,as.numeric(gen3_ecotype_frequency[i,]))
  pt0_3 <- c(pt0_3,as.numeric(gen2_ecotype_frequency[i,]))
}
hist(log(pt1_3 / pt0_3),breaks = 100)

# ## from generation 2 to generation 3
# log_p3_p2 <- c()
# gen3_index <- terminal_index$index3[terminal_index$generation1==1 & terminal_index$generation2==2 & terminal_index$generation3==3]
# gen2_index <- terminal_index$index2[terminal_index$generation1==1 & terminal_index$generation2==2 & terminal_index$generation3==3]
# gen3_ecotype_frequency <- ecotype_frequency[,gen3_index]
# gen2_ecotype_frequency <- ecotype_frequency[,gen2_index]
# 
# log_p3_p2 <- c()
# for(i in 1:num_ecotype){
#   log_p3_p2 <- c(log_p3_p2,log(as.numeric(gen3_ecotype_frequency[i,]) / as.numeric(gen2_ecotype_frequency[i,])))
# }
# hist(log_p3_p2,breaks = 100)

#table(log(pt1_3 / pt0_3) ==log_p3_p2 )



#log_p_ratio <- c(log_p1_p0,log_p2_p1,log_p3_p2)

pt1<- c(pt1_1,pt1_2,pt1_3)
pt0<- c(pt0_1,pt0_2,pt0_3)

## We can then construct the list M (ecotypes) 
## and list N (experimental sites)
## and list R (plots for each site)
## and list G (generation)

M <- c()
for(i in c(326,203,142)){
  for(j in 1:num_ecotype){
    M <- c(M,rep(founder_ecotype[j],i))
  }
}

table(unique(M) == founder_ecotype)



## generate N (site)
N <- c(rep(experimental_meta$site[experimental_meta$generation==1],231),
       rep(experimental_meta$site[experimental_meta$generation==2],231),
       rep(experimental_meta$site[experimental_meta$generation==3],231))


## generate R (plot)
R <- c(rep(experimental_meta$plot[experimental_meta$generation==1],231),
       rep(experimental_meta$plot[experimental_meta$generation==2],231),
       rep(experimental_meta$plot[experimental_meta$generation==3],231))

## generate G (generation)
G <- c(rep(1,length(pt1_1)),
       rep(2,length(pt1_2)),
       rep(3,length(pt1_3)))

length(R)
length(M)
length(N)
length(G)


mydata <- as.data.frame(matrix(nrow = length(pt1),ncol = 6))
colnames(mydata) <- c("pt1","pt0","ecotype","site","plot","generation")
mydata$pt1 <- pt1
mydata$pt0 <- pt0
mydata$ecotype <- M
mydata$site <- N
mydata$plot <- R
mydata$generation <- G

mydata <- mydata %>%
  left_join(ecotype_climate,by="ecotype")

mydata_1 <- mydata[mydata$generation==1,]
mydata_2 <- mydata[mydata$generation==2,]
mydata_3 <- mydata[mydata$generation==3,]


mydata_1 <- mydata_1 %>%
  left_join(site_climate_2018,by = "site")
mydata_2 <- mydata_2 %>%
  left_join(site_climate_2019,by = "site")
mydata_3 <- mydata_3 %>%
  left_join(site_climate_2020,by = "site")

mydata_site_merged <- bind_rows(mydata_1,mydata_2,mydata_3,)


write.table(mydata_site_merged,"data-intermediate/stabilizing_selection_data_2024Jun18_collectionsite_era5_allyears_pt1_pt0.txt",append = F,quote = F,sep = "\t",row.names = F)


#write.table(mydata_site_merged,"data-intermediate/stabilizing_selection_data_2024Jun18_collectionsite_era5_allyears.txt",append = F,quote = F,sep = "\t",row.names = F)



#plot(mydata_site_merged$bio1_ecotype - mydata_site_merged$bio1_site,mydata_site_merged$pt1)
