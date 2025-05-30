rm(list=ls())

library(dplyr)
library(geodata)
require(raster)
require(sp)


## This script use stabilizing selection model to predict Moi's 2019 nature paper

prefix <- "FUTURE_PREDICTION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

nature_data <- read.delim("data-external/Moi_2019_nature_data.txt")

stabilizing_selection <- read.delim("data-intermediate/FUTURE_PREDICTION_1001g_grenet_predicted_Wmax_Vs.txt",sep=" ")

stabilizing_selection_nature <- stabilizing_selection %>%
  filter(ecotype %in% nature_data$id) %>%
  filter(Vs_prd > 0)

nature_data <- nature_data %>%
  filter(id %in% stabilizing_selection_nature$ecotype)

## 513 accessions

## get the bio1 for these 513 accessions

ecotype_climate <- read.delim("data-external/bioclimvars_ecotypes_era5.csv",sep = ",")

ecotype_climate <- ecotype_climate %>%
  filter(ecotype %in% stabilizing_selection_nature$ecotype)

ecotype_info <- stabilizing_selection_nature %>%
  left_join(.,ecotype_climate[,c(1,2)],by="ecotype")

## extract madrid and tuebingen bio1

site <- c("Madrid","tuebingen")
lat <- c(40.40805,48.545809)
long <- c(-3.83535,9.042449)
coords <- data.frame(x=long, y = lat)
mycrs <- "+proj=longlat +datum=WGS84"

wc <- raster::stack(x = "data-external/wc2.1_2.5/wc2.1_2.5m_bio_1.tif")

pt <- sp::SpatialPoints(coords = coords,
                        proj4string = CRS(mycrs))
# extract from the given worldclim raster
wc_value <- raster::extract(x = wc, y = pt, method = 'simple', buffer=NULL)
wc_value
## madrid 14.117833 tuebingen 8.591833

## now calculate the fitness for spain and germany w/o water

stabilizing_selection_fitness <- function(origin,location,wmax,vs){
  w <- wmax * exp(-(origin - location)^2 / vs)
  return(w)
}

## spain

spain_fitness <- c()

for(i in 1:nrow(stabilizing_selection_nature)){
  spain_fitness <- c(spain_fitness,stabilizing_selection_fitness(ecotype_info$bio1_new[i],14.117833,ecotype_info$Wmax_prd[i],ecotype_info$Vs_prd[i]))
}

germany_fitness <- c()

for(i in 1:nrow(stabilizing_selection_nature)){
  germany_fitness <- c(germany_fitness,stabilizing_selection_fitness(ecotype_info$bio1_new[i],8.591833,ecotype_info$Wmax_prd[i],ecotype_info$Vs_prd[i]))
}

ecotype_info$spain_prd_fitness <- spain_fitness
ecotype_info$germany_prd_fitness <- germany_fitness

spain_true <- nature_data %>%
  filter(site=="madrid") %>%
  filter(water=="h") %>%
  filter(indpop=="i") %>%
  arrange(match(id,ecotype_info$ecotype))

#summary(lm(spain_true$Fitness~ecotype_info$spain_prd_fitness[ecotype_info$ecotype %in% spain_true$id]))
#cor.test(spain_true$Fitness,ecotype_info$spain_prd_fitness,method = "spearman",alternative = "greater")
#cor.test(spain_true$Fruits,ecotype_info$spain_prd_fitness,method = "spearman",alternative = "greater")
#cor.test(spain_true$Seeds,ecotype_info$spain_prd_fitness,method = "spearman",alternative = "greater")
cor.test(spain_true$Green,ecotype_info$spain_prd_fitness[ecotype_info$ecotype %in% spain_true$id],method = "spearman",alternative = "greater")
cor.test(spain_true$Fitness,ecotype_info$spain_prd_fitness[ecotype_info$ecotype %in% spain_true$id],method = "spearman",alternative = "greater")



summary(lm(rank(ecotype_info$spain_prd_fitness)~rank(spain_true$Green)))

germany_true <- nature_data %>%
  filter(site=="tuebingen") %>%
  filter(water=="l") %>%
  filter(indpop=="i") %>%
  arrange(match(id,ecotype_info$ecotype))

#summary(lm(germany_true$Green~ecotype_info$germany_prd_fitness[ecotype_info$ecotype %in% germany_true$id]))
cor.test(germany_true$Green,ecotype_info$germany_prd_fitness,method = "spearman",,alternative = "greater")


cor.test(germany_true$Green,spain_true$Green,method="spearman")

stabilizing_selection_fitness_ratio <- function(origin,location1,location2,vs){
  ratio <- exp((-(origin - location1 )^2 + (origin - location2)^2 ) / vs)
  return(ratio)
}

prd_ratio <- c()

for(i in 1:nrow(stabilizing_selection_nature)){
  prd_ratio <- c(prd_ratio,stabilizing_selection_fitness_ratio(ecotype_info$bio1_new[i],8.591833,14.117833,ecotype_info$Vs_prd[i]))
}

true_ratio <- germany_true$Green / spain_true$Green

#summary(lm(true_ratio[!is.infinite(true_ratio)]~prd_ratio[!is.infinite(true_ratio)]))
cor.test(true_ratio[!is.infinite(true_ratio)],prd_ratio[!is.infinite(true_ratio)],method = "spearman",alternative = "greater")

mean(spain_true$Fitness)
mean(germany_true$Fitness)
