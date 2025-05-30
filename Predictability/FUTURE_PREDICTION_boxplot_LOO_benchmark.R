rm(list=ls())

library(dplyr)

## future prediction

## This script generates the comparision figure for different prediction methods

prefix <- "FUTURE_PREDICTION_"
PATH <-  "~/Google Drive/My Drive/Research_Projects/GRENE-net_PHASE1/grenephase1-analyses/"
setwd(PATH)

climate <- read.delim("data-intermediate/Leave_one_out_prediction/climate_distance/Genomic_offset_site_climate_distance.txt")
go <- read.delim("data-intermediate/Leave_one_out_prediction/binomial_regression/genomic_offset_firstgen_binom_reg_results.csv",sep=",") %>%
  mutate(source = "genomic offset",
         r2 = r_squared,
         pearson_r = pearsonr,
         sp_r = sp_correlation) %>%
  dplyr::select(site,plot,pearson_r,sp_r,r2,source)
replicates <- read.delim("data-intermediate/Leave_one_out_prediction/replicates/Genomic_offset_plot_replication_prediction_summary.txt")
stabilizing_era5 <- read.delim("data-intermediate/Leave_one_out_prediction/stabilizing_selection_era5/Genomic_offset_stabilizing_selection_loo_era5_prediction_summary.txt")
stabilizing_prs <- read.delim("data-intermediate/Leave_one_out_prediction/stabilizing_selection_prs/Genomic_offset_stabilizing_selection_loo_prs_prediction_summary.txt")

parallelism <- read.delim("data-intermediate/generation_1_parallelism.txt")

df <- bind_rows(climate,go,replicates,stabilizing_era5,stabilizing_prs) %>%
  left_join(.,parallelism[parallelism$source=="snp",c(1,2)],by="site") %>%
  mutate(high_parallelism = mean > 0.3)
write.table(df,"data-intermediate/FUTURE_PREDICTION_LOO_accuracy.csv",append = F,quote = F,row.names = F,sep = ",")

unique(df$source)

df$source <- factor(df$source,levels = c("climate_distance","genomic offset","stabilizing_prs","stabilizing_era5","plot_mean"))

df_summary <- df %>%
  group_by(source)%>%
  summarise(mean_sp = mean(sp_r),
            mean_r = mean(pearson_r),
            mean_r2 = mean(r2))




ggplot(df,aes(x=source,y=r2,fill=high_parallelism))+
  geom_boxplot(width=0.5)+
  scale_fill_manual(values = c("grey80","grey30")) +
  coord_flip()+
  theme_minimal(base_size = 18)

ggsave("figs/FUTURE_PREDICTION_LOO_comparison_r2.pdf",device = "pdf",height = 5,width = 8,units = "in")


ggplot(df,aes(x=source,y=sp_r,fill=high_parallelism))+
  geom_boxplot(width=0.5)+
  scale_fill_manual(values = c("grey80","grey30")) +
  coord_flip()+
  theme_minimal(base_size = 18)
ggsave("figs/FUTURE_PREDICTION_LOO_comparison_sp_r.pdf",device = "pdf",height = 5,width = 8,units = "in")

ggplot(df,aes(x=source,y=pearson_r,fill=high_parallelism))+
  geom_boxplot(width=0.5)+
  scale_fill_manual(values = c("grey80","grey30")) +
  coord_flip()+
  theme_minimal(base_size = 18)
ggsave("figs/FUTURE_PREDICTION_LOO_comparison_pearson_r.pdf",device = "pdf",height = 5,width = 8,units = "in")

