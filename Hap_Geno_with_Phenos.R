hapgeno_with_pheno<-function(hap_out4) {
  library(tidyverse)
  load("/Volumes/MyFolder$/yunya418/Genes_Journal/Back_Elimination_on_sweeps/phenotype_with_sex_nobadsmpls.RData")
  source("/Volumes/MyFolder$/yunya418/Genes_Journal/Back_Elimination_on_sweeps/Sample_Geno.R")
  hap_geno_df_stch<-Sample_Geno(hap_out4)
  #hap_geno_df_stch<-hap_geno_df
  hap_geno_df_stch$Sample<-gsub("X","",hap_geno_df_stch$Sample)
  hap_geno_df_stch<-hap_geno_df_stch[hap_geno_df_stch$Sample %in% phenotype_with_sex_nobadsmpls_in_stch$ID,]
  hap_geno_df_stch<-hap_geno_df_stch[order(match(hap_geno_df_stch$Sample,phenotype_with_sex_nobadsmpls_in_stch$ID)),]
  hapgeno_with_pheno<-cbind(hap_geno_df_stch,phenotype_with_sex_nobadsmpls_in_stch)
}