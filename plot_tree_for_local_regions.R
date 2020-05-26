library(qqman);library(tidyverse);library(ape)
k10AFD0maf005_hws_chr4<-k10AFD0maf005_hws %>% filter(chr == 4, log >7.6)
k10AFD0maf005_hws_chr2<-k10AFD0maf005_hws %>% filter(chr == 2, pos %in% (1.07e8:1.09e8))

plot(k10AFD0maf005_hws_chr2$pos,k10AFD0maf005_hws_chr2$log)
with(subset(k10AFD0maf005_hws_chr2_4,21570503<pos&pos<21667322),points(pos,log,pch=16,col="red"))
text(21570000,12,labels = "GRIA2",col = "red",cex=1.2)

with(subset(k10AFD0maf005_hws_chr2_4,pos<21509410&pos<21558077),points(pos,log,pch=16,col="green"))
text(21509410,3,labels = "GLRB",col = "green",pos = 2,cex =1.2)

with(subset(k10AFD0maf005_hws_chr2_4,pos<21933586&pos>21915869),points(pos,log,pch=16,col="orange"))
text(21870000,6,labels = "ABCG2",col = "orange",pos = 4,cex =1.2)

with(subset(k10AFD0maf005_hws_chr2_4,pos<21952314&pos>21936057),points(pos,log,pch=16,col="blue"))
text(21900000,5,labels = "TREM144",col = "blue",pos = 4.3,cex =1.2)


k10AFD0maf005_hws$chr_pos<-paste0(k10AFD0maf005_hws$chr,"_",k10AFD0maf005_hws$pos)
k10AFD0maf005_hws<-k10AFD0maf005_hws[k10AFD0maf005_hws$chr_pos %in% af_hws53,]
geno_pheno_hw_chr9_3.6Mb_logQ2<-
  getPhe4STCH_geno_DF("/Volumes/MyFolder$/yunya418/Genes_Journal/Back_Elimination_on_sweeps/Hap_Fre_on_HWS_by_hapFLK/Chr4/STITCH_on_all28_Chrs_SCOREHWE_BGL4_2667spls_newname_Chr9_3.6Mb_159SNPs_genoDF.txt")
par(mfrow=c(3,2))
hws.tree=read.tree('chip25546m_ail5524212m_togather_imputed_BEAGLE4_gp09_smratio09_noneAFD0_noWS40_noP255_AIL2WS40_4_HW_pops_MAF005_plk_KING_tree.txt')
plot(hws.tree,align=T);axis(1,line=1);title(xlab='WGS data in HW lines')
fig_label("A",region = "plot",pos = "topleft",cex = 1.5)
## Chr1 82 Mb
hws_chr1_82Mb_tree<-read.tree("hws_tree_chr1_from_82232828_to_82482203_KING_tree.txt")
plot(hws_chr1_82Mb_tree,align=T);axis(1,line=1);title(xlab='Chr1 82Mb in HW')
fig_label("B",region = "plot",pos = "topleft",cex = 1.5)
## Chr1 165 Mb
hws_chr1_165Mb_tree<-read.tree("hws_tree_chr1_from_165806484_to_166109377_KING_tree.txt")
plot(hws_chr1_165Mb_tree,align=T);axis(1,line=1);title(xlab='Chr1 165Mb in HW')
fig_label("C",region = "plot",pos = "topleft",cex = 1.5)
## Chr1
hws_chr1_194Mb_tree<-read.tree("hws_tree_chr1_from_194544809_to_194735127_KING_tree.txt")
plot(hws_chr1_194Mb_tree,align=T);axis(1,line=1);title(xlab='Chr1 194Mb in HW')
fig_label("D",region = "plot",pos = "topleft",cex = 1.5)
## Chr2 16Mb
hws_chr2_16Mb_tree<-read.tree("hws_tree_chr2_16018852_to_17431874_KING_tree.txt")
plot(hws_chr2_16Mb_tree,align=T);axis(1,line=1);title(xlab='Chr2 16 ~ 17.4 Mb in HW lines')
fig_label("E",region = "plot",pos = "topleft",cex = 1.5)
## Chr2 107Mb
hws_chr2_107Mb_tree<-read.tree("hws_tree_chr2_from_107765315_to_108418074_KING_tree.txt")
plot(hws_chr2_107Mb_tree,align=T);axis(1,line=1);title(xlab='Chr2 107Mb in HW')
fig_label("F",region = "plot",pos = "topleft",cex = 1.5)
## Chr3 56Mb
hws_chr3_56Mb_tree<-read.tree("hws_tree_chr3_from_56150130_to_56366315_KING_tree.txt")
plot(hws_chr3_56Mb_tree,align=T);axis(1,line=1);title(xlab='Chr3 56Mb in HW')
fig_label("G",region = "plot",pos = "topleft",cex = 1.5)
## Chr3 71Mb
hws_chr3_71Mb_tree<-read.tree("hws_tree_chr3_from_71958103_to_72267239_KING_tree.txt")
plot(hws_chr3_71Mb_tree,align=T);axis(1,line=1);title(xlab='Chr3 71Mb in HW')
fig_label("H",region = "plot",pos = "topleft",cex = 1.5)
## Chr3 107Mb
hws_chr3_107Mb_tree<-read.tree("hws_tree_chr3_from_107354197_to_107434021_KING_tree.txt")
plot(hws_chr3_107Mb_tree,align=T);axis(1,line=1);title(xlab='Chr3 107Mb in HW')
fig_label("B",region = "plot",pos = "topleft",cex = 1.5)
## Chr4 21Mb
hws_chr4_21Mb_tree<-read.tree("hws_tree_chr4_from_21562649_to_21978762_KING_tree.txt")
plot(hws_chr4_21Mb_tree,align=T);axis(1,line=1);title(xlab='Chr4 21Mb in HW')
fig_label("C",region = "plot",pos = "topleft",cex = 1.5)
## Chr4 21Mb
hws_chr4_24Mb_tree<-read.tree("hws_tree_chr4_from_24141974_to_24817854_KING_tree.txt")
plot(hws_chr4_24Mb_tree,align=T);axis(1,line=1);title(xlab='Chr4 24Mb in HW')
fig_label("D",region = "plot",pos = "topleft",cex = 1.5)
## Chr4 35Mb
hws_chr4_35Mb_tree<-read.tree("hws_tree_chr4_from_35195793_to_35226658_KING_tree.txt")
plot(hws_chr4_35Mb_tree,align=T);axis(1,line=1);title(xlab='Chr4 35Mb in HW')
fig_label("E",region = "plot",pos = "topleft",cex = 1.5)
## Chr4 41Mb
hws_chr4_41Mb_tree<-read.tree("hws_tree_chr4_from_40873479_to_43113646_KING_tree.txt")
plot(hws_chr4_41Mb_tree,align=T);axis(1,line=1);title(xlab='Chr4 41Mb in HW')
fig_label("F",region = "plot",pos = "topleft",cex = 1.5)
## Chr4 51Mb
hws_chr4_51Mb_tree<-read.tree("hws_tree_chr4_from_51569335_to_51581439_KING_tree.txt")
plot(hws_chr4_51Mb_tree,align=T);axis(1,line=1);title(xlab='Chr4 51Mb in HW')
fig_label("G",region = "plot",pos = "topleft",cex = 1.5)
## Chr5 22Mb
hws_chr5_22Mb_tree<-read.tree("hws_tree_chr5_from_22672866_to_22679104_KING_tree.txt")
plot(hws_chr5_22Mb_tree,align=T);axis(1,line=1);title(xlab='Chr5 22Mb in HW')
fig_label("H",region = "plot",pos = "topleft",cex = 1.5)
## Chr6 1.3Mb
hws_chr6_1.3Mb_tree<-read.tree("hws_tree_chr6_from_1299200_to_1334638_KING_tree.txt")
plot(hws_chr6_1.3Mb_tree,align=T);axis(1,line=1);title(xlab='Chr6 1.3Mb in HW')
fig_label("B",region = "plot",pos = "topleft",cex = 1.5)
## Chr7 5Mb
hws_chr7_5Mb_tree<-read.tree("hws_tree_chr7_from_5113635_to_5124381_KING_tree.txt")
plot(hws_chr7_5Mb_tree,align=T);axis(1,line=1);title(xlab='Chr7 5Mb in HW')
fig_label("C",region = "plot",pos = "topleft",cex = 1.5)
## Chr8 5Mb
hws_chr8_2Mb_tree<-read.tree("hws_tree_chr8_from_2352132_to_2684657_KING_tree.txt")
plot(hws_chr8_2Mb_tree,align=T);axis(1,line=1);title(xlab='Chr8 2Mb in HW')
fig_label("D",region = "plot",pos = "topleft",cex = 1.5)
## Chr8 19Mb
hws_chr8_19Mb_tree<-read.tree("hws_tree_chr8_from_19580353_to_19660822_KING_tree.txt")
plot(hws_chr8_19Mb_tree,align=T);axis(1,line=1);title(xlab='Chr8 19Mb in HW')
fig_label("E",region = "plot",pos = "topleft",cex = 1.5)
## Chr9 13Mb
hws_chr9_13Mb_tree<-read.tree("hws_tree_chr9_from_11962587_to_13663676_KING_tree.txt")
plot(hws_chr9_13Mb_tree,align=T);axis(1,line=1);title(xlab='Chr9 13Mb in HW')
fig_label("F",region = "plot",pos = "topleft",cex = 1.5)
## Chr9 16Mb
hws_chr9_16Mb_tree<-read.tree("hws_tree_chr9_from_16684306_to_17498371_KING_tree.txt")
plot(hws_chr9_16Mb_tree,align=T);axis(1,line=1);title(xlab='Chr9 16Mb in HW')
fig_label("G",region = "plot",pos = "topleft",cex = 1.5)
## Chr10 1.8
hws_chr10_1.8Mb_tree<-read.tree("hws_tree_chr10_from_1864672_to_1867346_KING_tree.txt")
plot(hws_chr10_1.8Mb_tree,align=T);axis(1,line=1);title(xlab='Chr10 1.8Mb in HW')
fig_label("H",region = "plot",pos = "topleft",cex = 1.5)
## Chr11 10Mb
hws_chr11_10Mb_tree<-read.tree("hws_tree_chr11_from_10647635_to_11278128_KING_tree.txt")
plot(hws_chr11_10Mb_tree,align=T);axis(1,line=1);title(xlab='Chr11 10Mb in HW')
fig_label("B",region = "plot",pos = "topleft",cex = 1.5)
## Chr12 14Mb
hws_chr12_14Mb_tree<-read.tree("hws_tree_chr12_from_14412479_to_14520447_KING_tree.txt")
plot(hws_chr12_14Mb_tree,align=T);axis(1,line=1);title(xlab='Chr12 14Mb in HW')
fig_label("C",region = "plot",pos = "topleft",cex = 1.5)
## Chr18 1.8Mb
hws_chr18_1.8Mb_tree<-read.tree("hws_tree_chr18_from_1814832_to_1954119_KING_tree.txt")
plot(hws_chr18_1.8Mb_tree,align=T);axis(1,line=1);title(xlab='Chr18 1.8Mb in HW')
fig_label("D",region = "plot",pos = "topleft",cex = 1.5)
## Chr20 10Mb
hws_chr20_10Mb_tree<-read.tree("hws_tree_chr20_from_10692416_to_10737192_KING_tree.txt")
plot(hws_chr20_10Mb_tree,align=T);axis(1,line=1);title(xlab='Chr20 10.6 ~ 10.7 Mb in HW lines')
fig_label("E",region = "plot",pos = "topleft",cex = 1.5)
## Chr24 4Mb
hws_chr24_4Mb_tree<-read.tree("hws_tree_chr24_from_4046230_to_4106152_KING_tree.txt")
plot(hws_chr24_4Mb_tree,align=T);axis(1,line=1);title(xlab='Chr24 4Mb in HW')
fig_label("F",region = "plot",pos = "topleft",cex = 1.5)
################ Trees in LWS logq 2 ####################
par(mfrow=c(4,2))
## Whole Genome
lws.tree<-read.tree("chip25546m_ail5524212m_togather_imputed_BEAGLE4_gp09_smratio09_noneAFD0_LW_maf005_tree.txt")
plot(lws.tree,align=T);axis(1,line=1);title(xlab='WGS data in LW lines')
fig_label("A",region = "plot",pos = "topleft",cex = 1.5)
## Chr4 29Mb
hws_chr4_29Mb_tree<-read.tree("lws_tree_chr4_from_29007140_to_29428159_KING_tree.txt")
plot(hws_chr4_29Mb_tree,align=T);axis(1,line=1);title(xlab='Chr4 29Mb in LW')
fig_label("D",region = "plot",pos = "topleft",cex = 2)
## Chr6 16Mb
hws_chr6_16Mb_tree<-read.tree("lws_tree_chr6_from_16942912_to_17731291_KING_tree.txt")
plot(hws_chr6_16Mb_tree,align=T);axis(1,line=1);title(xlab='Chr6 16.9 ~ 17.7 Mb in LW lines')
## Chr11 7Mb
hws_chr11_7Mb_tree<-read.tree("lws_tree_chr11_from_7252453_to_7578662_KING_tree.txt")
plot(hws_chr11_7Mb_tree,align=T);axis(1,line=1);title(xlab='Chr11 7Mb in LW')
## Chr13 7Mb
hws_chr13_3Mb_tree<-read.tree("lws_tree_chr13_from_3213046_to_3493989_KING_tree.txt")
plot(hws_chr13_3Mb_tree,align=T);axis(1,line=1);title(xlab='Chr13 3.2 ~ 3.4 Mb in LW lines')
## Chr18 7Mb
hws_chr18_9Mb_tree<-read.tree("lws_tree_chr18_from_9541522_to_9597536_KING_tree.txt")
plot(hws_chr18_9Mb_tree,align=T);axis(1,line=1);title(xlab='Chr18 9Mb in LW')
## Chr21 4Mb
hws_chr21_4Mb_tree<-read.tree("lws_tree_chr21_from_4911726_to_4962699_KING_tree.txt")
plot(hws_chr21_4Mb_tree,align=T);axis(1,line=1);title(xlab='Chr21 4Mb in LW')


################ Trees in LWS logq 1.3####################
par(mfrow=c(1,2))
## Whole Genome
lws.tree<-read.tree("chip25546m_ail5524212m_togather_imputed_BEAGLE4_gp09_smratio09_noneAFD0_LW_maf005_tree.txt")
plot(lws.tree,align=T);axis(1,line=1);title(xlab='Whole Genome in LW')
fig_label("A",region = "plot",pos = "topleft",cex = 1.5)
## Chr1 37 Mb
hws_chr1_37Mb_tree<-read.tree("lws_tree_chr1_from_35202489_to_37589459_KING_tree.txt")
plot(hws_chr1_37Mb_tree,align=T);axis(1,line=1);title(xlab='Chr1 37Mb in LW')
fig_label("B",region = "plot",pos = "topleft",cex = 1.5)
## Chr1 65 Mb
hws_chr1_65Mb_tree<-read.tree("lws_tree_chr1_from_65331948_to_65677701_KING_tree.txt")
plot(hws_chr1_65Mb_tree,align=T);axis(1,line=1);title(xlab='Chr1 65Mb in LW')
fig_label("C",region = "plot",pos = "topleft",cex = 1.5)
## Chr1 177 Mb
hws_chr1_177Mb_tree<-read.tree("lws_tree_chr1_from_177438078_to_177438879_KING_tree.txt")
plot(hws_chr1_177Mb_tree,align=T);axis(1,line=1);title(xlab='Chr1 177Mb in LW')
fig_label("D",region = "plot",pos = "topleft",cex = 1.5)
## Chr2 55 Mb
hws_chr2_55Mb_tree<-read.tree("lws_tree_chr2_from_55068124_to_56076475_KING_tree.txt")
plot(hws_chr2_55Mb_tree,align=T);axis(1,line=1);title(xlab='Chr2 55Mb in LW')
fig_label("E",region = "plot",pos = "topleft",cex = 1.5)
## Chr2 106 Mb
hws_chr2_106Mb_tree<-read.tree("lws_tree_chr2_from_106406529_to_106502538_KING_tree.txt")
plot(hws_chr2_106Mb_tree,align=T);axis(1,line=1);title(xlab='Chr2 106Mb in LW')
fig_label("F",region = "plot",pos = "topleft",cex = 1.5)
## Chr2 134 Mb
hws_chr2_134Mb_tree<-read.tree("lws_tree_chr2_from_134160794_to_146134161_KING_tree.txt")
plot(hws_chr2_134Mb_tree,align=T);axis(1,line=1);title(xlab='Chr2 134Mb in LW')
fig_label("G",region = "plot",pos = "topleft",cex = 1.5)
## Chr3 Mb
hws_chr3_66Mb_tree<-read.tree("lws_tree_chr3_from_66629672_to_67225643_KING_tree.txt")
plot(hws_chr3_66Mb_tree,align=T);axis(1,line=1);title(xlab='Chr3 66Mb in LW')
fig_label("H",region = "plot",pos = "topleft",cex = 1.5)
## Chr4 29Mb
hws_chr4_29Mb_tree<-read.tree("lws_tree_chr4_from_29007140_to_29428159_KING_tree.txt")
plot(hws_chr4_29Mb_tree,align=T);axis(1,line=1);title(xlab='Chr4 29Mb in LW')
fig_label("B",region = "plot",pos = "topleft",cex = 1.5)
## Chr6 16Mb
hws_chr6_16Mb_tree<-read.tree("lws_tree_chr6_from_16934374_to_17739115_KING_tree.txt")
plot(hws_chr6_16Mb_tree,align=T);axis(1,line=1);title(xlab='Chr6 16Mb in LW')
fig_label("C",region = "plot",pos = "topleft",cex = 1.5)
## Chr6 24Mb
hws_chr6_24Mb_tree<-read.tree("lws_tree_chr6_from_24711581_to_25439205_KING_tree.txt")
plot(hws_chr6_24Mb_tree,align=T);axis(1,line=1);title(xlab='Chr6 24Mb in LW')
fig_label("D",region = "plot",pos = "topleft",cex = 1.5)
## Chr7 0.5Mb
hws_chr7_0.5Mb_tree<-read.tree("lws_tree_chr7_from_594783_to_1591031_KING_tree.txt")
plot(hws_chr7_0.5Mb_tree,align=T);axis(1,line=1);title(xlab='Chr7 0.5Mb in LW')
fig_label("E",region = "plot",pos = "topleft",cex = 1.5)
## Chr7 2Mb
hws_chr7_2Mb_tree<-read.tree("lws_tree_chr7_from_2834822_to_2976468_KING_tree.txt")
plot(hws_chr7_2Mb_tree,align=T);axis(1,line=1);title(xlab='Chr7 2Mb in LW')
fig_label("F",region = "plot",pos = "topleft",cex = 1.5)
## Chr8 26Mb
hws_chr8_26Mb_tree<-read.tree("lws_tree_chr8_from_26903258_to_27630695_tree.txt")
plot(hws_chr8_26Mb_tree,align=T);axis(1,line=1);title(xlab='Chr8 26Mb in LW')
fig_label("G",region = "plot",pos = "topleft",cex = 1.5)
## Chr11 7Mb
hws_chr7_11Mb_tree<-read.tree("lws_tree_chr11_from_7199883_to_7607333_KING_tree.txt")
plot(hws_chr7_11Mb_tree,align=T);axis(1,line=1);title(xlab='Chr11 7Mb in LW')
fig_label("H",region = "plot",pos = "topleft",cex = 1.5)
## Chr13 3Mb
hws_chr13_3Mb_tree<-read.tree("lws_tree_chr13_from_3211917_to_3503070_KING_tree.txt")
plot(hws_chr13_3Mb_tree,align=T);axis(1,line=1);title(xlab='Chr13 3Mb in LW')
fig_label("B",region = "plot",pos = "topleft",cex = 1.5)
## Chr14 12Mb
hws_chr14_12Mb_tree<-read.tree("lws_tree_chr14_from_12627083_to_12644395_KING_tree.txt")
plot(hws_chr14_12Mb_tree,align=T);axis(1,line=1);title(xlab='Chr14 12Mb in LW')
fig_label("C",region = "plot",pos = "topleft",cex = 1.5)
## Chr15 6Mb
hws_chr14_6Mb_tree<-read.tree("lws_tree_chr15_from_6670393_to_6780570_KING_tree.txt")
plot(hws_chr14_6Mb_tree,align=T);axis(1,line=1);title(xlab='Chr15 6Mb in LW')
fig_label("D",region = "plot",pos = "topleft",cex = 1.5)
## Chr18 2Mb
hws_chr18_2Mb_tree<-read.tree("lws_tree_chr18_from_2844293_to_2872417_KING_tree.txt")
plot(hws_chr18_2Mb_tree,align=T);axis(1,line=1);title(xlab='Chr18 2Mb in LW')
fig_label("E",region = "plot",pos = "topleft",cex = 1.5)
## Chr18 2Mb
hws_chr18_9Mb_tree<-read.tree("lws_tree_chr18_from_9528753_to_9702003_KING_tree.txt")
plot(hws_chr18_9Mb_tree,align=T);axis(1,line=1);title(xlab='Chr18 9Mb in LW')
fig_label("F",region = "plot",pos = "topleft",cex = 1.5)
## Chr21 2Mb
hws_chr21_2Mb_tree<-read.tree("lws_tree_chr21_from_2139490_to_2184206_KING_tree.txt")
plot(hws_chr21_2Mb_tree,align=T);axis(1,line=1);title(xlab='Chr21 2Mb in LW')
fig_label("G",region = "plot",pos = "topleft",cex = 1.5)
## Chr21 4.3Mb
hws_chr21_4.3Mb_tree<-read.tree("lws_tree_chr21_from_4305097_to_5856178_KING_tree.txt")
plot(hws_chr21_4.3Mb_tree,align=T);axis(1,line=1);title(xlab='Chr21 4.3Mb in LW')
fig_label("H",region = "plot",pos = "topleft",cex = 1.5)
## Chr26 3.6Mb
hws_chr26_3.6Mb_tree<-read.tree("lws_tree_chr26_from_3600714_to_3624219_KING_tree.txt")
plot(hws_chr26_3.6Mb_tree,align=T);axis(1,line=1);title(xlab='Chr26 3.6Mb in LW')
fig_label("B",region = "plot",pos = "topleft",cex = 1.5)


#### Plot HWS and LWS together: just relaxed regions
par(mfrow=c(2,3),las=1,mar=c(4.5,3,2,2))
## Whole Genome in LWS
lws.tree<-read.tree("chip25546m_ail5524212m_togather_imputed_BEAGLE4_gp09_smratio09_noneAFD0_LW_maf005_tree.txt")
plot(lws.tree,align=T,cex = 1);axis(1,line=1,cex.axis=1.2)
title(xlab='Whole Genome Data in LW lines',cex.lab=1.2)
fig_label("A",region = "plot",pos = "topleft",cex = 2)
## Chr6 16Mb
lws_chr6_16Mb_tree<-read.tree("lws_tree_chr6_from_16934374_to_17739115_KING_tree.txt")
plot(lws_chr6_16Mb_tree,align=T,cex = 1);axis(1,line=1,cex.axis=1.2)
title(xlab='Chr6 16.9 ~ 17.7 Mb in LW lines',cex.lab=1.2)
fig_label("B",region = "plot",pos = "topleft",cex = 2)
## Chr6 24Mb
plot(hws_chr6_24Mb_tree,align=T,cex = 1);axis(1,line=1,cex.axis=1.2)
title(xlab='Chr6 24.7 ~ 25.4 Mb in LW lines',cex.lab=1.2)
fig_label("C",region = "plot",pos = "topleft",cex = 2)
### in HWS
hws.tree=read.tree('chip25546m_ail5524212m_togather_imputed_BEAGLE4_gp09_smratio09_noneAFD0_noWS40_noP255_AIL2WS40_4_HW_pops_MAF005_plk_KING_tree.txt')
plot(hws.tree,align=T,cex = 1);axis(1,line=1)
title(xlab='Whole Genome Data in HW lines',cex.lab=1.2)
fig_label("D",region = "plot",pos = "topleft",cex = 2)
## Chr8 5Mb
hws_chr8_2Mb_tree<-read.tree("hws_tree_chr8_from_2352132_to_2684657_KING_tree.txt")
plot(hws_chr8_2Mb_tree,align=T,cex=1);axis(1,line=1)
title(xlab='Chr8 2.3 ~ 2.7 Mb in HW lines',cex.lab=1.2)
fig_label("E",region = "plot",pos = "topleft",cex = 2)
hws_chr20_10Mb_tree<-read.tree("hws_tree_chr20_from_10692416_to_10737192_KING_tree.txt")
plot(hws_chr20_10Mb_tree,align=T,cex=1);axis(1,line=1);
title(xlab='Chr20 10.6 ~ 10.7 Mb in HW lines',cex.lab=1.2)
fig_label("F",region = "plot",pos = "topleft",cex = 2)
## Chr9 13Mb
hws_chr9_13Mb_tree<-read.tree("hws_tree_chr9_from_11962587_to_13663676_KING_tree.txt")
plot(hws_chr9_13Mb_tree,align=T,cex=1);axis(1,line=1)
title(xlab='Chr9 11.9 ~ 13.7 Mb in HW lines',cex.lab=1.2)
fig_label("F",region = "plot",pos = "topleft",cex = 2)




