# Install bowtie2 and samtools:
# sudo apt-get install bowtie2 samtools
# Install cutadapt:
# sudo pip3 install cutadapt
# Install SRA toolkit:
# https://github.com/ncbi/sra-tools/wiki/HowTo:-Binary-Installation
# Don't be tempted to apt-get install sra-tools.  That version is broken.

# Data from Gao et al. (2018)
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6284/
# Data from Rooijers et al. (2013)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48933

ERRctrl=('ERR2208504' 'ERR2208505') # 	control osteosarcoma 143B, control normal HEK293T
ERRpatient=('ERR2208506' 'ERR2208507' 'ERR2208508') # mitochondrial disease
SRRs=('SRR935452' 'SRR935453') # Cybrid control RP rep1, Cybrid control RP rep2
#roots=("${ERRctrl[@]}" "${SRRs[@]}")
roots=("${ERRctrl[@]}" "${ERRpatient[@]}" "${SRRs[@]}")

#wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR220/006/ERR2208506/ERR2208506.fastq.gz &
#wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR220/007/ERR2208507/ERR2208507.fastq.gz &
#wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR220/008/ERR2208508/ERR2208508.fastq.gz &

adapt_Short='TCGTATGCCGTCTTCTGCTTG'
adapt_Long='TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG'

for root in "${roots[@]}"
do
	echo ${root}

	# Use appropriate adapter sequence for each sample
	if [ "$root" == "SRR935452" ]; then
	  adapt="${adapt_Short}"
	else
	  adapt="${adapt_Long}"
	fi

	echo ${adapt}

	# Download and uncompress annotated human genome  
	#wget ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
	#gunzip Homo_sapiens.GRCh38.ncrna.fa.gz
	# Filter RNA species
	#awk '/^>/ {P=($0~"gene_biotype:rRNA")||($0~"gene_biotype:rRNA_pseudogene")||($0~"gene_biotype:Mt_rRNA")} {if(P) print} ' Homo_sapiens.GRCh38.ncrna.fa > Homo_sapiens.GRCh38.rrna.fa
	# Download tRNA sequences
	#wget http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/hg19-tRNAs.fa
	# Combine rRNA and tRNA sequences into one fasta file and index
	#cat Homo_sapiens.GRCh38.rrna.fa hg19-tRNAs.fa > hgRNA.fa
	#mkdir hgRNA
	#bowtie2-build hgRNA.fa hgRNA/hgRNA

	# Download sequence and build index for human mtDNA
	#curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text" >mtDNA.fa
	#mkdir mtDNA
	#bowtie2-build mtDNA.fa mtDNA/mtDNA

	# Download reads in fastq format
	if [[ " ${SRRs[@]} " =~ " ${root} " ]]; then
	  echo ${root}
	  #fastq-dump ${root}
	fi
	
	if [[ " ${ERRs[@]} " =~ " ${root} " ]]; then
	  echo ${root}
	  #wget -c 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR220/005/${root}/${root}.fastq.gz' # URL building slightly wrong (not all 005)
	  #gunzip ${root}.fastq.gz
	fi
	
	# Cut approriate adapter sequences
	cutadapt -a ${adapt} -O 12 -m 20 -j 22 ${root}.fastq -o ${root}_trimmed.fastq
	
	# Drop all reads that do not end in AAAAA
	awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if(seq ~ /[A]{5,15}$/){tailstart = match(seq, /[A]{5,15}$/); print header, substr(seq,0,tailstart-1), qheader, substr(qseq,0,tailstart-1)}}' < ${root}_trimmed.fastq > ${root}_filtered.fastq

	# Align reads to decoy Trna/Rrna from nuclear and MT genome
	bowtie2 -p 22 -D20 -R 10 -N 1 -L 20 -i C,1 --un ${root}_screened.fastq -x hgRNA/hgRNA -U ${root}_filtered.fastq -S ${root}_screened_out.sam

	# Align to mitochondrial genome
	bowtie2 -p 22 -D20 -R 10 -N 1 -L 20 -i C,1 -x mtDNA/mtDNA -U ${root}_screened.fastq -S ${root}_aligned_mito.sam

	# Filter out unaligned reads, secondary alignments, aligned reads with 0 quality and sort.
	samtools view -Sb ${root}_aligned_mito.sam -u| samtools view -f 0 -q 1 - -u|samtools sort - -f ${root}_polyA_header.bam
	samtools index ${root}_polyA_header.bam ${root}_polyA_header.bai
	samtools view ${root}_polyA_header.bam > ${root}_polyA_header.sam
	
done