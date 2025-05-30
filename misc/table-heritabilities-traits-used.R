# Get the heritability


st_all <- readRDS(file="~/Google Drive/natvar/data/all_gwa_stats_table.rda")
statsSS <- st_all[!is.na(st_all$stressstrategy),]
p1 <- ggplot(statsSS, aes(x = PVE_median, fill = stressstrategy,)) +
  geom_density(alpha=0.25) + xlim(0,1) + xlab("h2") +
  theme(legend.position = c(0.35, .7), legend.key.size =unit(.7, 'cm'),
        legend.text = element_text(size=7), legend.title = element_blank()) +
  scale_fill_manual("",values = c("Avoidance"="red","Escape"='green4', "Tolerance"="navy"))

#statsPC <- st_all[!is.na(st_all$phenotypecategory),]
p2 <- ggplot(st_all, aes(x = pve.x), show.legend=F) +
  geom_density(alpha=0.25, fill="dodgerblue1") + xlim(0,1) +xlab("h2") + ylab("")
p3 <- ggplot(st_all, aes(x = PVE_median), show.legend=F) +
  geom_density(alpha=0.25, fill="dodgerblue1") + xlim(0,1) +xlab("h2") + ylab("")

p_combined <- ggplot(st_all) +
  geom_density(aes(x = pve.x, fill = "pve - lmm"), alpha = 0.25, fill = "dodgerblue1") +
  geom_density(aes(x = PVE_median, fill = "pve - bslmm"), alpha = 0.25, fill = "green4") +
  xlim(0, 1) +
  theme(legend.position = c(0.35, .7), legend.key.size =unit(.7, 'cm'),
        legend.text = element_text(size=7), legend.title = element_blank()) +
  labs(title = "green - bslmm, blue - lmm",
       x = "h2",
       y = "")
p_combined


st_all %>% head


st_all %>%  dplyr::filter(phenotype=="FT10") %>% summarize(hsummary= paste0(format(h_mean,digits=3), "[",format(h_2.5.,digits=3), "-",format(h_97.5.,,digits=3),"]"))

st_all %>%  dplyr::filter(phenotype=="Delta_13C") %>% summarize(hsummary= paste0(format(h_mean,digits=3), "[",format(h_2.5.,digits=3), "-",format(h_97.5.,,digits=3),"]"))

st_all %>%  dplyr::filter(phenotype=="DSDS50") %>% summarize(hsummary= paste0(format(h_mean,digits=3), "[",format(h_2.5.,digits=3), "-",format(h_97.5.,,digits=3),"]"))

st_all %>%  dplyr::filter(phenotype=="stomata_density") %>% summarize(hsummary= paste0(format(h_mean,digits=3), "[",format(h_2.5.,digits=3), "-",format(h_97.5.,,digits=3),"]"))

# st_all %>%  dplyr::filter(phenotype=="rosette_DM") %>% summarize(hsummary= paste0(format(h_mean,digits=3), "[",format(h_2.5.,digits=3), "-",format(h_97.5.,,digits=3),"]"))


st_all %>%  dplyr::filter(phenotype=="Leaf_Area") %>% summarize(hsummary= paste0(format(h_mean,digits=3), "[",format(h_2.5.,digits=3), "-",format(h_97.5.,,digits=3),"]"))

st_all %>%  dplyr::filter(phenotype=="Root_horizontal_index_avg") %>% summarize(hsummary= paste0(format(h_mean,digits=3), "[",format(h_2.5.,digits=3), "-",format(h_97.5.,,digits=3),"]"))
st_all %>%  dplyr::filter(phenotype=="Relative_root_length") %>% summarize(hsummary= paste0(format(h_mean,digits=3), "[",format(h_2.5.,digits=3), "-",format(h_97.5.,,digits=3),"]"))


allphenotypenames<-st_all$phenotype %>% unique
allphenotypenames[grep("Delta", allphenotypenames)]
length(allphenotypenames)

st_all %>%  dplyr::filter(phenotype %in% c("FT10","Delta_13C","DSDS50","stomata_density","Root_horizontal_index_avg", "Leaf_Area") ) %>% summarize(hsummary= mean(h_mean))

