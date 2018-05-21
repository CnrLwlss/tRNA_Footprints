ERRctrl=('ERR2208504' 'ERR2208505') # 	control osteosarcoma 143B, control normal HEK293T
ERRpatient=('ERR2208506' 'ERR2208507' 'ERR2208508') # mitochondrial disease
SRRs=('SRR935452' 'SRR935453') # Cybrid control RP rep1, Cybrid control RP rep2
roots=("${ERRctrl[@]}" "${ERRpatient[@]}" "${SRRs[@]}")

for root in "${roots[@]}"
do
	echo ${root}
	samtools view -Sb  ${root}_polyA_header.sam  >  ${root}_polyA_header.bam
	samtools depth ${root}_polyA_header.bam > ${root}_polyA_header.coverage
	#samtools flagstat ${root}_polyA_header.bam
	samtools view -F 0x04 -c ${root}_polyA_header.bam
done


	#|\
	#awk 'BEGIN { prev_chr="NC_012920.1";prev_pos=0;} { if($1==prev_chr && prev_pos+1!=int($2)) {for(i=prev_pos+1;i<int($2);++i) {printf("%s\t%d\t0\n",$1,i);}} print; prev_chr=$1;prev_pos=int($2);}' > ${root}_polyA_header.coverage
