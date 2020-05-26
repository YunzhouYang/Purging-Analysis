#!/bin/bash -l
  
#SBATCH -A  snic2018-8-148
#SBATCH -p node
#SBATCH -C fat
#SBATCH -t 7-00:00:00
#SBATCH --mail-user=yang.yunzhou@imbim.uu.se
#SBATCH --mail-type=ALL

## $1 >>> tmp_dir 
## $2 >>> chr in CM000xxx.4 format
## $3 >>> posfile
## $4 >>> outputdir in 
## $5 >>> nCores
## $6 >>> Numeric Chr,eg: Chr1 or Chr23 
cd $7
Rscript /proj/uppstore2018093/private/F2_fq2bam/Impute_Chr1toChr28_with_F1/STITCH_Rcpp.R $1 $2 $3 $4 $5 $6 $7

bcftools filter -i 'INFO/INFO_SCORE > 0.4 & INFO/HWE > 1e-6' ./output_dir/"STITCH_"$2"_"$3"_K10nGen19.vcf.gz" -Oz -o ./output_dir/"STITCH_"$2"_"$3"_K10nGen19_SCOREHWE.vcf.gz"
