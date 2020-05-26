#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
require(data.table,lib.loc="/home/yunzhou/yyz_stitchR/")
require(rrbgen,lib.loc="/home/yunzhou/yyz_stitchR/")
require(STITCH,lib.loc="/home/yunzhou/yyz_stitchR/")
setwd(args[7])
file_name=paste0("STITCH_",args[2],"_",args[6],"_K10nGen19.vcf.gz")
core=as.numeric(as.character(args[5]))
bamlist="/proj/uppstore2018093/private/F2_fq2bam/Low_coverage_fastq_F0_to_F19/all_F0toF18_3610Samples_20200129_ThibautCheckSamples_bamlist.txt"
STITCH(tempdir = args[1], chr = args[2], bamlist = bamlist, posfile = args[3], outputdir = args[4], K = 10, nGen = 19, nCores = core, output_filename = file_name)
