Hap_Fre_Ori_TagSNP<-function(GT_file_AIL,GT_file_STCH,region,parts,left,right,out_dir,filterpositons) {
  library(ggplot2);library(ggpubr);library(reshape2)
  #GT_file_AIL<-"chip25546m_ail5524212m_togather_imputed_BEAGLE4_gp09_smratio09_noneAFD0_noWS40_noP255_AIL2WS40_Chr4_21Mb_152SNPs_genoDF.txt"
  #region<-"Chr4";parts="Peak1_21Mb";left=21618187;right=21675123;out_dir="./"
  hapchr6_ail<-read.csv2(GT_file_AIL,header = T,sep = "\t",stringsAsFactors = F)
  hapchr6_ail<-hapchr6_ail[,-c(1,3,6:9)]
  #print("Extract geotypes from AIL founders files......")
  for (r in 1:nrow(hapchr6_ail)) {
    for (c in 4:ncol(hapchr6_ail)) {
      hapchr6_ail[r,c]<-unlist(strsplit(as.character(hapchr6_ail[r,c]),":"))[1]
    }
  }
   
  #print("replace number genos with AGCT......")
  for (r in 1:nrow(hapchr6_ail)) {
    #print(r)
    for (c in 4:ncol(hapchr6_ail)) {
      if (hapchr6_ail[r,c]=="0|1" | hapchr6_ail[r,c]=="1|0") {
        hapchr6_ail[r,c]<-paste0(hapchr6_ail[r,2]," ",hapchr6_ail[r,3])
      } else if (hapchr6_ail[r,c]=="0|0") {
        hapchr6_ail[r,c]<-paste0(hapchr6_ail[r,2]," ",hapchr6_ail[r,2])
      } else if (hapchr6_ail[r,c]=="1|1") {
        hapchr6_ail[r,c]<-paste0(hapchr6_ail[r,3]," ",hapchr6_ail[r,3])
      }
    }
  }
  hapchr6_ail<-hapchr6_ail[,-c(2,3)]
  ############# Read GT files from STITCH pipeline for 2667 Samples
  #GT_file_STCH<-"STITCH_on_all28_Chrs_SCOREHWE_BGL4_2667spls_newname_hws_chr4_Bonf_21Mb_152Snps_genoDF.txt"
  hapchr6_stch<-read.csv2(GT_file_STCH,header = T,sep = "\t",stringsAsFactors = F)
  hapchr6_stch<-hapchr6_stch[,-c(1,3,6:9)]
  #print("Extract geotypes from STITCH files for 2667 Samples......")
  for (r in 1:nrow(hapchr6_stch)) {
    for (c in 4:ncol(hapchr6_stch)) {
      hapchr6_stch[r,c]<-unlist(strsplit(as.character(hapchr6_stch[r,c]),":"))[1]
    }
  }
  
  #print("replace number genos with AGCT......")
  for (r in 1:nrow(hapchr6_stch)) {
    #print(r)
    for (c in 4:ncol(hapchr6_stch)) {
      if (hapchr6_stch[r,c]=="0|1" | hapchr6_stch[r,c]=="1|0") {
        hapchr6_stch[r,c]<-paste0(hapchr6_stch[r,2]," ",hapchr6_stch[r,3])
      } else if (hapchr6_stch[r,c]=="0|0") {
        hapchr6_stch[r,c]<-paste0(hapchr6_stch[r,2]," ",hapchr6_stch[r,2])
      } else if (hapchr6_stch[r,c]=="1|1") {
        hapchr6_stch[r,c]<-paste0(hapchr6_stch[r,3]," ",hapchr6_stch[r,3])
      }
    }
  }
  hapchr6_stch<-hapchr6_stch[,-c(2,3)]
  ############## Mrege files from STITCH and AIL datasets ########
  hapchr6_stch<-hapchr6_stch[hapchr6_stch$POS %in% hapchr6_ail$POS,]
  hapchr6_ail<-hapchr6_ail[hapchr6_ail$POS %in% hapchr6_stch$POS,]
  hapchr6<-merge(hapchr6_stch,hapchr6_ail,by.x = "POS",by.y = "POS")
  ################################################################
  ################# Filter Positions if you want delete some SNPs#
  hapchr6<-hapchr6[hapchr6$POS %in% filterpositons,]
  ################################################################
  #print("remove 2/3 columns")
  rownames(hapchr6)<-hapchr6$POS;hapchr6<-hapchr6[,-1]
  ###########
  #hapchr6<-cbind(hapchr6,hapchr6_chip197s)
  ###########
  #
  hapchr6<-hapchr6 %>% t() %>% as.matrix() 
  
  ## freq
  genomatrix<-hapchr6
  hapL<-c()
  hapL_out<-c()
  for(i in 1:nrow(genomatrix)){
    geno<-unlist(strsplit(as.character(genomatrix[i,])," "))
    # get two hap together
    hapL<-c(hapL,paste(geno[seq(1,length(geno),2)],collapse=""))
    hapL<-c(hapL,paste(geno[seq(2,length(geno),2)],collapse=""))
    
    hapL_out<-rbind(hapL_out,t(paste(geno[seq(1,length(geno),2)],sep="\t")))
    hapL_out<-rbind(hapL_out,t(paste(geno[seq(2,length(geno),2)],sep="\t")))
  }
  ##name
  newrownames<-c()
  colnames(hapL_out)<-colnames(genomatrix)
  for ( r in 1:nrow(genomatrix)) {
    newrownames<-c(newrownames,paste0(rownames(genomatrix)[r],".1"),paste0(rownames(genomatrix)[r],".2"))
  }
  rownames(hapL_out)<-newrownames
  ##calculate haplotype frequecies in each pop
  hapL_out2<-as.data.frame(hapL_out)
  hapL_out3<-apply(format(hapL_out2), 1,paste, collapse="")
  hapL_out2$hap1<-hapL_out3
  # LWS40 pop
  hapL_out2$pop<-rownames(hapL_out2)
  dim(hapL_out2)
  a=dim(hapL_out2)[2]
  b=dim(hapL_out2)[2]-1
  c=dim(hapL_out2)[2]-2
  hapL_out2<-hapL_out2[,c(a,b,1:c)]
  ## Get pop names for STITCH samples
  for (r in 1:(2667*2)) {
    hapL_out2$pop[r]<-unlist(strsplit(hapL_out2$pop[r],"\\."))[1]
    ln<-length(unlist(strsplit(hapL_out2$pop[r],'')))
    x1<-paste(unlist(strsplit(hapL_out2$pop[r],''))[(ln-1):ln],collapse = "")
    hapL_out2$pop[r]<-paste0("F_",x1)
  }
  ## Get pop names for AIL samples
  for (r in (2667*2+1):dim(hapL_out2)[1]) {
    #r=1
    hapL_out2$pop[r]<-unlist(strsplit(hapL_out2$pop[r],"_"))[2]
    hapL_out2$pop[r]<-unlist(strsplit(hapL_out2$pop[r],'\\.'))[1]
  }
  
  
  ## use hapL_out5 to create haplotype fasta file to do cluster and Plot Heatmap
  hapL_out5<-hapL_out2
  hapL_out5$sample_names<-row.names(hapL_out5)
  hapL_out5$hap_new_name<-NA
  ## rename those haplotypes in HW_LW
  hapL_out4<-hapL_out2[,1:2]
  table_hap<-table(hapL_out4$hap1)
  for (m in 1:length(names(table(hapL_out4$hap1)))) {
    new_hap_name<-paste0("Hap",m)
    #print(new_hap_name)
    pos<-which(hapL_out4$hap1==names(table_hap)[m])
    #print(pos)
    for (rr in pos) {
      #print(rr)
      hapL_out4$hap1[rr]<-new_hap_name
      hapL_out5$hap_new_name[rr]<-new_hap_name
    }
  }
  d=dim(hapL_out5)[2]-1
  e=dim(hapL_out5)[2]
  f=dim(hapL_out5)[2]-2
  hapL_out5<-hapL_out5[,c(1,2,d,e,3:f)]
  ##save hapl_out5 for heatmap(a little tricky to save new variable created by Paste0 function)
  assign(paste0("Hap_Seq_on_",region,"_",parts,"_from_",left,"_to_",right,"_hapl_out5"),hapL_out5)
  xx<-paste0("Hap_Seq_on_",region,"_",parts,"_from_",left,"_to_",right,"_hapl_out5")
  save(list = xx,file = paste(out_dir,xx,".RData",sep = ""))
  ##calculate frequencies
  pop_hap_fre<-melt(t(table(hapL_out4)));pop_hap_fre$Hap_total<-NA
  pop1<-levels(as.factor(pop_hap_fre$pop))
  for (i in 1:length(pop1)) {
    #i=1
    #print(pop1[i])
    pos<-which(pop_hap_fre$pop==pop1[i])
    pop_hap_fre$Hap_total[pos]<-dim(hapL_out4[hapL_out4$pop==pop1[i],])[1]
  }
  ##Frequency of Haplotypes
  pop_hap_fre$Hap_fre=pop_hap_fre$value/pop_hap_fre$Hap_total
  LWS_pops<-c("LWS40","LWS50","LWS53","LWR9")
  HWS_pops<-c("HWS40","HWS50","HWS53","HWR9")
  #pop_hap_fre_197<-pop_hap_fre %>% filter(pop %in% c(LWS_pops,HWS_pops))
  #
  
  ##Add group info
  for (r in 1:nrow(pop_hap_fre)) {
    pop_hap_fre$Group[r]<-paste0(unlist(strsplit(as.character(pop_hap_fre$pop[r]),NULL))[1:2],collapse = "")
  }
  ####################### Save pop_hap_fre #####################
  print("Save pop_hap_fre ......")
  assign(paste0("Hap_Fre_on_197_HWSR_LWSR_and_2667_AIL_samples",region,"_from_",left,"_to_",right),pop_hap_fre)
  xz<-paste0("Hap_Fre_on_197_HWSR_LWSR_and_2667_AIL_samples",region,"_from_",left,"_to_",right)
  save(list = xz,file = paste(out_dir,xz,".RData",sep = ""))
  ###########################################################################
  
  ####################### Save Hap_Genotypes(Hap_Genos) #####################
  ##save which samples have which haplotypes,like 2132_LWS40.1 has two haplotypes and they are:
  #####2132_LWS40.1 Hap24
  #####2132_LWS40.2 Hap26
  hapL_out4$sample_names<-rownames(hapL_out4)
  print("Save samples with Hapltoype names ......")
  assign(paste0("Hap_Genos_on_197_HWSR_LWSR_and_2667_AIL_samples_",region,"_from_",left,"_to_",right,"hapL_out4"),hapL_out4)
  xy<-paste0("Hap_Genos_on_197_HWSR_LWSR_and_2667_AIL_samples_",region,"_from_",left,"_to_",right,"hapL_out4")
  save(list = xy,file = paste(out_dir,xy,".RData",sep = ""))
  ###########################################################################
  LW_pos<-pop_hap_fre[pop_hap_fre$pop %in% LWS_pops,]
  HW_pos<-pop_hap_fre[pop_hap_fre$pop %in% HWS_pops,]
  pop_hap_fre_AIL<-pop_hap_fre %>% filter(!(pop %in% c(LWS_pops,HWS_pops)))
  lw<-ggplot(LW_pos,aes(x=pop,y=Hap_fre,group=factor(LW_pos$hap1)))+scale_y_continuous(limits = c(0,1))+
     geom_line(aes(col=factor(LW_pos$hap1)))+geom_text(aes(label=factor(LW_pos$hap1)),hjust=0, vjust=0)+
     geom_point(aes(col=factor(LW_pos$hap1)))+
     labs(x=NULL,y="Haplotype frequecies in each population")+
     theme(plot.title = element_text(hjust = 0.5,size = 10),
           axis.title.x = element_text(size = 10),
           axis.title.y = element_text(size = 10),
           axis.text = element_text(size = 10))+scale_color_discrete(name="Haplotypes")
  
  #Hight weight pops
   hw<-ggplot(HW_pos,aes(x=pop,y=Hap_fre,group=factor(HW_pos$hap1)))+scale_y_continuous(limits = c(0,1))+
    geom_line(aes(col=factor(HW_pos$hap1)))+
    geom_point(aes(col=factor(HW_pos$hap1)))+geom_text(aes(label=factor(HW_pos$hap1)),hjust=0, vjust=0)+
    labs(x=NULL,y="Haplotype frequecies in each population")+
    theme(plot.title = element_text(hjust = 0.5,size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text = element_text(size = 10))+scale_color_discrete(name="Haplotypes")
  

    figure<-ggarrange(lw, hw, common.legend = TRUE, legend = "none", align = "h")
    fig_title<-paste0(region,"_from_",left,"_to_",right,"_",parts,"_in_LW_HW")
    final_fig<-annotate_figure(figure,top = text_grob(fig_title))
    final_fig
  
  
  print("Plot Haplotype Frequency......")
    hapfre_in_ail<-ggplot(pop_hap_fre_AIL,aes(x=pop,y=Hap_fre,group=factor(pop_hap_fre_AIL$hap1)))+scale_y_continuous(limits = c(0,1))+
      geom_line(aes(col=factor(pop_hap_fre_AIL$hap1)))+geom_text(aes(label=factor(pop_hap_fre_AIL$hap1)),hjust=0, vjust=0)+
      geom_point(aes(col=factor(pop_hap_fre_AIL$hap1)))+
      labs(x=NULL,y="Haplotype frequecies in each population")+
      theme(plot.title = element_text(hjust = 0.5,size = 10),legend.position = "none",
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            axis.text = element_text(size = 10))+scale_color_discrete(name="Haplotypes")
  #Hight weight pops
  ## if frequency of haplotypes in all pops is 0, then remove them.
  Hap_fre_in_AIL197_and_STCH2667<-ggarrange(final_fig,hapfre_in_ail,ncol = 1,nrow = 2)
  return(Hap_fre_in_AIL197_and_STCH2667)
}
