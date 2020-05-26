getPhe4STCH_geno_DF<-function(genoDF) {
  library(tidyverse)
  load("/Volumes/MyFolder$/yunya418/Genes_Journal/Back_Elimination_on_sweeps/phenotype_with_sex_nobadsmpls.RData")
  #genoDF<-"STITCH_on_all28_Chrs_SCOREHWE_BGL4_2667spls_4334Snps_with_24ZYJ_fixedSNPs_newname_HWS_genoDF.txt"
  genos<-read.table(genoDF,header = T,stringsAsFactors = F,sep = "\t") %>% select(-c(3,6:9)) %>% 
    mutate(Chr=CHROM) %>% select(1,Chr,everything())

  for (m in 6:ncol(genos)) {
    genos[,m]<-str_split_fixed(genos[,m],":",3)[,1]
  }
  
  chrs_GG5<-read.table("/Volumes/MyFolder$/yunya418/Genes_Journal/Back_Elimination_on_sweeps/chr_id.match_modified_by_yyz.txt",header = T,stringsAsFactors = F,sep = "\t") %>% 
    select(c(2,4)) %>% filter(Name %in% (1:28))
  STCH_chrs<-levels(as.factor(genos$Chr))
  for (p in 1:length(STCH_chrs)) {
    #p=1
    pos11<-which(chrs_GG5$INSDC==STCH_chrs[p])
    genos$Chr<-gsub(STCH_chrs[p],chrs_GG5$Name[pos11],genos$Chr)
  }

  smples<-colnames(genos)[6:2672];smples<-gsub("X","",smples);colnames(genos)[6:2672]<-smples
  genos<-genos %>% mutate(chr_pos=paste0(genos$Chr,"_",genos$POS)) %>% select(chr_pos,everything()) %>% select(-c(2:6))

  for (m in 2:2668) {
    #print(m)
    #m=2
    genos[,m]<-as.character(genos[,m])
    genos[,m]<-gsub("0\\|0","0",genos[,m])
    genos[,m]<-gsub("0\\|1","1",genos[,m])
    genos[,m]<-gsub("1\\|0","1",genos[,m])
    genos[,m]<-gsub("1\\|1","2",genos[,m])
    genos[,m]<-gsub("\\.|\\.",NA,genos[,m])
  }
  ## missing is not allowed in BE analysis
  genos_nomiss_t<-genos %>% na.omit() %>% t()
  
  colnames(genos_nomiss_t)<-as.character(unlist(genos_nomiss_t[1,]))
  genos_nomiss_t<-genos_nomiss_t[-1,]
  genos_nomiss_t<-as.data.frame(genos_nomiss_t)
  
  sm<-row.names(genos_nomiss_t)
  genos_nomiss_t$sm<-sm
  genos_nomiss_t<-genos_nomiss_t[,c(ncol(genos_nomiss_t),1:(ncol(genos_nomiss_t)-1))]
  
  phenotype_with_sex_nobadsmpls_in_stch<-
    phenotype_with_sex_nobadsmpls_in_stch[order(match(phenotype_with_sex_nobadsmpls_in_stch$ID,genos_nomiss_t$sm)),]
  phe_geno_2667<-cbind(phenotype_with_sex_nobadsmpls_in_stch,genos_nomiss_t)
  return(phe_geno_2667)
}

