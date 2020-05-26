get_HapGenos4for012copy_withPheno<-function(HapGenowithPheno,hapltoype) {
  #HapGenowithPheno<-hws21SNPschr9_hapgenoPheno
  #hapltoype<-"Hap63"
  copy0<-paste0("no",hapltoype)
  copy1<-paste0(hapltoype,"/-")
  copy2<-paste0(hapltoype,"/",hapltoype)
  hws11SNPschr11_hapgenoPheno_2Copy<-HapGenowithPheno[HapGenowithPheno$Genotypes==copy2,]
  hws11SNPschr11_hapgenoPheno_2Copy$Genotypes<-copy2
  hws11SNPschr11_hapgenoPheno_1Copy<-HapGenowithPheno[grep(hapltoype,HapGenowithPheno$Genotypes),]
  hws11SNPschr11_hapgenoPheno_1Copy<-hws11SNPschr11_hapgenoPheno_1Copy[hws11SNPschr11_hapgenoPheno_1Copy$Genotypes!=copy2,]
  hws11SNPschr11_hapgenoPheno_1Copy$Genotypes<-copy1
  hws11SNPschr11_hapgenoPheno_0Copy<-HapGenowithPheno[grep(hapltoype,HapGenowithPheno$Genotypes),]
  hws11SNPschr11_hapgenoPheno_0Copy<-HapGenowithPheno[!(HapGenowithPheno$Genotypes%in%hws11SNPschr11_hapgenoPheno_0Copy$Genotypes),]
  hws11SNPschr11_hapgenoPheno_0Copy$Genotypes<-copy0
  hws11SNPschr11_hapgenoPheno_3Genos<-rbind(hws11SNPschr11_hapgenoPheno_2Copy,hws11SNPschr11_hapgenoPheno_0Copy,hws11SNPschr11_hapgenoPheno_1Copy)
  return(hws11SNPschr11_hapgenoPheno_3Genos)
}
