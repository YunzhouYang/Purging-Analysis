## Input 4 arguments here:
	## sys.argv[1]: original imputed vcf files by BEAGLE4;
 	## sys.argv[2]: new txt file containing the failed SNP sites info;

import re
import pandas as pd
import sys
import gzip

all_txt=open(sys.argv[2],"w")
all_txt.close()
all_txt=open(sys.argv[2],"a")
x0=0
with gzip.open(sys.argv[1],"rt") as all_cohort:
	for line in all_cohort:
		#print(line)
		if re.match("^##", line):
			continue
		elif re.match("^#CHROM", line):
			continue
		elif re.match("\n", line):
			continue
		else:
			line=line.rstrip("\n")
			con=re.split("\t",line)
			con_length=len(con)
			x=0
			for n in range(9,con_length):
				gp_vector=re.split(":",con[n])[2]
				gp_max=max(re.split(",",gp_vector))
				if float(gp_max) < 0.95:
					x=x+1
				else:
					continue
			ratio=x/197
			if ratio>0:
				chr_pos=con[0]+"\t"+con[1]
				print(chr_pos,file=all_txt)
		x0+=1
		print(str(x0))
all_cohort.close()

all_txt.close()
