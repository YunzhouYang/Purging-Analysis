# Purging
# Introduction
Virginia chicken lines were an well-designed experimental populations which were used to study complicated tratis, such as body weight, fitness. Bi-directional selection on 56-day body weight was used to produce high body weight lines (HWS) and low body weight lines (LWS). Meanwhile, relaxed lines (HWR/LWR) were also created as a contemporary lines to HWS and LWS. By comparing genomes among XWS and XWR populations, we detected purged regions when selection was stopped. We created this repository to publish our data which were used in our recent work on [*GENES*](https://www.mdpi.com/2073-4425/11/6/630) to see full text.

## Imputation

### Impute chip data to Whole genome sequencing density
Genotypes by SNP-chip were already published by [Matts, *et al*, 2013](https://www.g3journal.org/content/3/12/2305.short). Within these samples, 59 individals which were founders for advanced intercross lines (AIL) were sequenced in high depth (~30X) by [Yanjun and Thibaut, 2019](https://link.springer.com/article/10.1186/s12711-019-0487-1). Here, we imputed genotypes for all lines (HWS40/50/53,HWR9 and LWS40/50/53, LWR9) to a whole-genome wide density by [*BEAGLE4*](https://faculty.washington.edu/browning/beagle/b4_0.html).

  - Download and install [*BEAGLE4*](https://faculty.washington.edu/browning/beagle/b4_0.html) software.
  - Impuation command: `java -jar beagle.r1399.jar gt=file_to_be_imputed out=imputed_file_name`, you will get an imputed file named as: **imputed_file_name.vcf.gz**.
  - Filtration: only SNPs with GP values above 0.9 in more than 90% samples will be kept. We wrote a python script. Python script used: `extract_sites_with_gp095sm1_197samples.py`.
  - Separate dataset to LW and HW sets, thereafter, calculate allele frequency for each generation in LW and HW, respectively.
  - Calculate allele frequency difference between XWS40/50/53 and XWR9 (X=L or H).
  - Extract sites with AFD >=0.5. If AFD (>0.5) between each of XWS40/50/53 and XWR9 on sites is found, use `intersect` function in R to choose sites.
  - Ather stpes above, you get a imputed and phased vcf file for 197 XWS and XWR samaples, like **all_197_Imputed.vcf.gz**.

### Impute low coverage sequencing data by *STITCH*
Variants infos are required for STITCH and  [Yanjun and Thibaut, 2019](https://link.springer.com/article/10.1186/s12711-019-0487-1) have provided these infomations. You can choose to use a reference or not and we didn't use it. Steps are: 
  
  - Select K values: K was considered as number of ancestral haplotypes in the population you are working on. We already have ~30X depth sequencing data (for 60 F0 samples) which means we have known genotypes with high confidence.
      - Choose 10 F0 and downsample them to ~0.4X by [*seqtk*](https://github.com/lh3/seqtk). Commands are: `seqtk sample -s100 read1.fq 10000 > sub1.fq;seqtk sample -s100 read2.fq 10000 > sub2.fq`. Change 10000(reads number) as you want.
      - Extract chromosome 28 for these 10 F0 and other 2667 Samples.
      - Produce bam files for these fastq files.
      - Give different K values to `STITCH_Rcpp.R` and run the script to get the imputation results.
        
        - Prepare large storage folder for computation as temparay files are also big.
        - Use **Fat Node** if you: (1)have very big chromosomes like chr1 in chicken (~200Mb) and (2) many samples (> 1000).
        
      - Compare imputed and 30x sequenced results for those 10 F0 samples.
      - Plot results and select K.
As you known, K is always good if it si bigger than 5. Theoretically, bigger K will give you better results but it will take very long time to run and you have compromise between feasibility and accurarcy. Anyway, K around 10 is good enough for current study and it took us 7 days to impute Chr1 for 3511 samples before.

Once you know which K is good, you can run the STITCH for other 27 chromosomes followed by: 
  
  - filtration by `bcftools filter -i 'INFO/INFO_SCORE > 0.4 & INFO/HWE > 1e-6' -Oz -o xx_SCOREHWE.vcf.gz`
  - filtration by `vcftools --gzvcf xx_SCOREHWE.vcf.gz --remove badsamples.txt --recode --stdout | bgzip > xx_SCOREHWE_noBadSamples.vcf.gz` (STITCH will record some samples when there are something wrong with it).
  - Phase it by BEAGLE as `java -jar beagle.r1399.jar gt= out= xx_SCOREHWE_noBadSamples.vcf.gz out=xx_SCOREHWE_noBadSamples_BGL4`.
  - Extract good sites by `extract_sites_with_gp095sm1_197samples.py` (change sample size: 197 to 3511).


## Selective sweep analysis
As seletions (directional and relaxed) were carrried out in XWS/XWR (X=H or L), [*hapFLK*](https://forge-dga.jouy.inra.fr/projects/hapflk) was used to detect sweeps owing to these selections.
  
  - Install hapFLK and python2.7.
  - transform the vcf files to fam/ped or binary formats by PLINK. As chicken genome has more than 23 chromosome pairs, **--aec --chr-set 28** should be added to your commands.
  - Run hapFLK for each chromosome respectviely and merge results (Set K=10).
  - Download *scale_chi2_hapflk.py* to calculated pvalues based on hapFLK results.
  - Calculate FDR threshold to select selected and relaxed regions. You cannot tell which regions is under ongoing selection or relaxed selection currently unless you continue the following steps.
  - Extract significant regions and plot local trees. 
    
      - Extract regions: get regions by vcftools and transform it to binary format for hapFLK.
        
        - `vcftools --gzvcf xxx.vcf.gz --from-bp 1111 --to-bp 2222 --stdout --recode | > bgzip 11112222.vcf.gz`
        - `vcftools --gzvcf 11112222.vcf.gz --plink-tped --out 11112222_plk`
        - `plink --tfile 11112222_plk --make-bed --aec --chr-set 28 --out 11112222_plk_bin`
      
      - Use hapFLK to produce tree files: `hapFLK --bfile 11112222_plk_bin -p 11112222_plk_binKIN`, you will get *11112222_plk_binKIN_tree.txt*
  - Plot phelogeny trees by  `plot_tree_for_local_regions.R`. 
    
    - Edit folder path to adapt it to your own environment.
    - Download and install required packages.
    - Use [*Figure_Lable.R*](https://www.r-bloggers.com/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/) to add lables for figures. I revised this function: `figure_label_topleft.R`.
    
Based on trees, you can find regions which are under ongoing selection  and which are under relexed selection.

## Select Tag SNPs within purging regions.

We have steps to do this: 
  - Select variants: (1) those located within purging regions from local tree, (2) having AFD >=0.5.
  - Extract these SNPs and store them in a vcf-like file. vcf-like file is a vcf file that has no annotation lines which start with ## and this file start from #CHROM.
  - Use `getPhe4STCH_geno_DF.R` function to produce a genotype file with phenotypes. Phenotype was **F2_to_F15_2667Samples_phetypes.txt**.
  - Use Back-Elimination method to choose SNPs which are independently related to body weight. This R package was published by [Yanjun *et al*, 2017](https://academic.oup.com/mbe/article/34/10/2678/3952745). You can download it from https://github.com/yanjunzan/BE. Now you get the Tag SNPs.
  - Extract these Tag SNPs from **xx_SCOREHWE_noBadSamples.vcf.gz** and **all_197_Imputed.vcf.gz**. Store them in  vcf-like file for hapltoype analysis, like **xxchr6_16Mb_TagSNP_STITCH_genoDF.txt** and **xxchr6_16Mb_TagSNP_AIL_genoDF.txt**. A recommended scritp is like this: `vcftools --gzvcf xx_SCOREHWE_noBadSamples.vcf.gz --positions TagSNP_Chr_Pos.txt --recode -stdout |grep -v -e "##.*" > xxchr6_16Mb_TagSNP_genoDF.txt`. Then,use `sed -i s/^#// xxchr6_16Mb_TagSNP_genoDF.txt` to remove *#* in the file.

## Haplotype analysis
This part we will use some function to 
  - extract haplotypes form **xxchr6_16Mb_TagSNP_genoDF.txt** file.
  - calculate haplotype frequency in XWS40/50/53 and XWR9 (X=L or H).
  - Find out purging haplotypes according to frequenc in selected and relaxed linex.
  - associations between puging haplotypes and body weight.
  
### Extract haplotypes and calculate hapltoype frequency
  We have **xxchr6_16Mb_TagSNP_STITCH_genoDF.txt** and **xxchr6_16Mb_TagSNP_AIL_genoDF.txt** file and they contain phased genotypes for all samples. What we need to do is to extract haplotypes from these phased files. The R function used here is called `Hap_Fre_in_HWLWSR_withTagSNPs.R`. It will extract hapltoypes, calculate frequencies and plot haplotype frequencies in XWS/XWR and AIL. `Hap_Fre_Ori_TagSNP(GT_file_AIL, GT_file_STCH,region,parts,filterpositons,left,right,out_dir)` and this function will produce three dataframe:
  - Hap_Fre_on_xxx_pop_hap_fre: haplotype frequency for all hapltypes in each generation.
  - Hap_Genos_on_xxx_hapL_out4: copy number of hapltoypes for each individual
  - Hap_Seq_on_xxx_hapl_out5: halotype sequences.

Once you get the haplotype frequency, you will know which is selected or purged haplotype.

### Associations between purged hapltypes and 56-day body weight

We have several steps to do association analysis:
  - Assign pheotypes to hapltpye genotypes which contain different copies of hapltoype of interest. `Hap_Geno_with_Phenos.R` was used. Just give **Hap_Genos_on_xxx_hapL_out4** file to it.
  - Do associations on purging or selected hapltoyeps by using `get_HapGenos4for012copy_withPheno.R` function. You just give the name of hapltype you're interested in.
  - Use a `aov(bw56~haplotype_genotypes+sex+generation)` to test the difference among samples with different copies of purged or selected hapltoyeps.
  - Visualizations: Plot in the way you like. I used `ggplot+geom_boxplot`.

  
  
  








