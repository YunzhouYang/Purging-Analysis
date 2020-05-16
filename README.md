# Purging
# Introduction
Virginia chicken lines were an well-designed experimental populations which were used to study complicated tratis, such as body weight, fitness. Bi-directional selection on 56-day body weight was used to produce high body weight lines (HWS) and low body weight lines (LWS). Meanwhile, relaxed lines (HWR/LWR) were also created as a contemporary lines to HWS and LWS. By comparing genomes among XWS and XWR populations, we detected purged regions when selection was stopped. We created this repository to publish our data which were used in our recent work on *GENES* [Click here] to see full text.

## Imputation
Genotypes by SNP-chip were already published by [Matts, *et al*, 2013](https://www.g3journal.org/content/3/12/2305.short). Within these samples, 59 individals which were founders for advanced intercross lines (AIL) were sequenced in high depth (30X) by [Yanjun and Thibaut, 2019](https://link.springer.com/article/10.1186/s12711-019-0487-1). Here, we imputed genotypes for all lines (HWS40/50/53,HWR9 and LWS40/50/53, LWR9) to a whole-genome wide density by [*BEAGL4*](https://faculty.washington.edu/browning/beagle/b4_0.html).

  - Download and install [BEAGLE4](https://faculty.washington.edu/browning/beagle/b4_0.html) software.
  - Impuation command: `java -jar beagle.r1399.jar gt=file_to_be_imputed out=imputed_file_name`
  - Filtration: only SNPs with GP values above 0.9 in more than 90% samples will be kept. We wrote a python script.