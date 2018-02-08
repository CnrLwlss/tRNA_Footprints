# Install bowtie2 and samtools:
# apt-get install bowtie2 samtools
# Install cutadapt:
# pip3 install cutadapt
# Install SRA toolkit:
# https://github.com/ncbi/sra-tools/wiki/HowTo:-Binary-Installation
# Don't be tempted to apt-get install sra-tools.  That version is broken.

root='SRR935452'
#root='SRR935453'

adapt_Rooijers='TCGTATGCCGTCTTCTGCTTG'
adapt_Gao='TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG'

# Collect sequences and build index of human rRNA and tRNA 
#wget ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
#gunzip Homo_sapiens.GRCh38.ncrna.fa.gz
#awk '/^>/ {P=($0~"gene_biotype:rRNA")||($0~"gene_biotype:rRNA_pseudogene")||($0~"gene_biotype:Mt_rRNA")} {if(P) print} ' Homo_sapiens.GRCh38.ncrna.fa > Homo_sapiens.GRCh38.rrna.fa
#wget http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/hg19-tRNAs.fa
#cat Homo_sapiens.GRCh38.rrna.fa hg19-tRNAs.fa > hgRNA.fa
#mkdir hgRNA
#bowtie2-build hgRNA.fa hgRNA/hgRNA

# Collect sequences and build index of human mtDNA
#curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text" >mtDNA.fa
#mkdir mtDNA
#bowtie2-build mtDNA.fa mtDNA/mtDNA

#fastq-dump ${root}

# For SRR935453: "One or more adapter sequences may be incomplete: 91% of the bases preceding the sequence are C.  Add C at start of adapter?
# Should there be two adapters?  One 3' and one 5'?
cutadapt -a ${adapt_Rooijers} -O 12 -m 20 -j 5 ${root}.fastq -o ${root}_trimmed.fastq

# Align reads to decoy Trna/Rrna from nuclear and MT genome (4% alignment?)
bowtie2 -p 5 -D20 -R 10 -N 1 -L 20 -i C,1 --un screened.fastq -x hgRNA/hgRNA -U ${root}_trimmed.fastq -S screened_out.sam

# Align to mitochondrial genome
bowtie2 -p 5 -D20 -R 10 -N 1 -L 20 -i C,1 -x mtDNA/mtDNA -U screened.fastq -S ALIGNMITO.sam

# Filter out multi mapping reads, convert to a sorted BAM file
#samtools view -q 1 -bS screened_out.sam -u|samtools sort - -f screened.sorted.bam
#samtools view -q 1 -S screened_out.sam -u|samtools sort - -f screened.sorted.sam
#samtools view -q 1 -S screened_out.sam