
extract_genotypes_vcf<-function(myvcf){
  gt <- extract.gt(myvcf, element = 'GT', as.numeric = TRUE) %>% 
    t() %>% 
    data.frame
  return(gt)
}
