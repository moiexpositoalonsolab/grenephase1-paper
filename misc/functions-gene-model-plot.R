
gene_model_plot<-function(atname="AT2G27030"){

if (exists("tairraw")) {
  print("The object 'tairraw' exists.")
} else {
  print("Reading tairraw")
  tairraw<-read.table("data-external/TAIR10_GFF3_genes_transposons.gff")
}

tairsub<-
  tairraw %>% 
  dplyr::filter(grepl(atname,V9))
# tairsub

exon_data <- tairsub %>%
  filter(V3 == "exon") %>%
  mutate(start = as.numeric(V4),
         end = as.numeric(V5),
         isoform = V9)  # Extract exon information
ggplot()+
  # geom_point(aes(y=-log10(p), x=pos))+
  # geom_vline( xintercept = 11531967)+
  # geom_vline( xintercept = 11534358)+
  # xlim(c(11531967-500,11534358+500))+
  geom_rect(data=exon_data,  # adding the genes
            aes(xmin = start, xmax = end,
                ymin = -as.numeric(as.factor(isoform)) - 0.4, ymax = -as.numeric(as.factor(isoform)) -0.1
            ), fill = "grey") +
  theme_classic()+
  theme(axis.line.y = element_blank())+
  theme(axis.text.y = element_blank())+
  theme(axis.ticks.y = element_blank())

}