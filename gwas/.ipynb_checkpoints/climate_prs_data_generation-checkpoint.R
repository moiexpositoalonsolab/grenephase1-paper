## this script generates the climate data for climate prs

setwd("~/scratch/grenet/climate_prs/")

climate <- read.delim("../metadata/1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv",sep=",")
info <- read.delim("../metadata/1001g_regmap_grenet_ecotype_info_corrected_2024May16.csv",sep=",")

grenet_accession <- read.delim("../ecotypes_names.txt",header=F)$V1
accession_1001g <- read.delim("1001g/1001g_grenet_accessions.txt",header=F)$V1
accession_1001g <- accession_1001g[accession_1001g %in% climate$ecotypeid]
write.table(accession_1001g,"1001g/1001g_grenet_accessions_climate.txt",append = F,quote = F,row.names = F,col.names = F)

accession_1001g_regmap <- read.delim("1001g_regmap/1001g_regmap_grenet_accessions.txt",header=F)$V1
accession_1001g_regmap <- accession_1001g_regmap[accession_1001g_regmap %in% climate$ecotypeid ]
write.table(accession_1001g_regmap,"1001g_regmap/1001g_regmap_grenet_accessions_climate.txt",append = F,quote = F,row.names = F,col.names = F)

USA_samples <- info$ecotype_id[info$country=="USA"]

table(grenet_accession%in%accession_1001g_regmap)

## 1001g training (exclude grenet and USA samples) for hapfm
train_1001g <- accession_1001g[!accession_1001g %in% grenet_accession & !accession_1001g%in%USA_samples]
train_1001g_climate <- climate[climate$ecotypeid %in% train_1001g,]
train_1001g_climate <- train_1001g_climate[match(train_1001g,train_1001g_climate$ecotypeid),]

for(i in 1:19){
  name = paste0("1001g/hapfm/pheno_training/bio",i,".txt")
  write.table(train_1001g_climate[,i+2],file = name,append = F,quote = F,row.names = F,col.names = F)
}
write.table(train_1001g,"1001g/hapfm/training_accessions.txt",append = F,quote = F,row.names = F,col.names = F)

## 1001g pheno data for bslmm
climate_1001g <- climate[climate$ecotypeid %in% accession_1001g,]
climate_1001g <- climate_1001g[match(accession_1001g,climate_1001g$ecotypeid),]
NA_1001g <- accession_1001g[accession_1001g %in% grenet_accession | accession_1001g%in%USA_samples]
climate_1001g$bio1[climate_1001g$ecotypeid%in%NA_1001g] <- "NA"

for(i in 1:19){
  name = paste0("1001g/bslmm/bioclim/bio",i,".txt")
  climate_1001g[climate_1001g$ecotypeid%in%NA_1001g,i+2] <- "NA"
  write.table(climate_1001g[,i+2],file = name,append = F,quote = F,row.names = F,col.names = F)
}

## 1001g_regmap training (exclude grenet and USA samples) for hapfm
train_1001g_regmap <- accession_1001g_regmap[!accession_1001g_regmap %in% grenet_accession & !accession_1001g_regmap%in%USA_samples]
train_1001g_regmap_climate <- climate[climate$ecotypeid %in% train_1001g_regmap,]
train_1001g_regmap_climate <- train_1001g_regmap_climate[match(train_1001g_regmap,train_1001g_regmap_climate$ecotypeid),]

for(i in 1:19){
  name = paste0("1001g_regmap/hapfm/pheno_training/bio",i,".txt")
  write.table(train_1001g_regmap_climate[,i+2],file = name,append = F,quote = F,row.names = F,col.names = F)
}
write.table(train_1001g_regmap,"1001g_regmap/hapfm/training_accessions.txt",append = F,quote = F,row.names = F,col.names = F)

## 1001g_regmap pheno data for bslmm
climate_1001g_regmap <- climate[climate$ecotypeid %in% accession_1001g_regmap,]
climate_1001g_regmap <- climate_1001g_regmap[match(accession_1001g_regmap,climate_1001g_regmap$ecotypeid),]
NA_1001g_regmap <- accession_1001g_regmap[accession_1001g_regmap %in% grenet_accession | accession_1001g_regmap%in%USA_samples]

for(i in 1:19){
  name = paste0("1001g_regmap/bslmm/bioclim/bio",i,".txt")
  climate_1001g_regmap[climate_1001g_regmap$ecotypeid%in%NA_1001g_regmap,i+2] <- "NA"
  write.table(climate_1001g_regmap[,i+2],file = name,append = F,quote = F,row.names = F,col.names = F)
}

