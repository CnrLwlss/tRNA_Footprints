export PATH=$PATH

root='GSE48933'
chrlen=${#root}
subr=${root:0:(chrlen-3)}

# Download...
#wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/${subr}nnn/${root}/suppl/${root}_RAW.tar

# Untar...
mkdir -p ${root}
tar -xvf ${root}_RAW.tar --directory ${root}
# Remove adapter
#cutadapt -a ????? -O 12 -m ????.fastq -o ????_trimmed.fastq

# Align reads to decoy Trna/Rrna from nuclear and MT genome
#bowtie2 -p 4 -D20 -R 10 -N 1 -L 20 -i C,1 --un screened.fastq -x prescreen -U ????_trimmed.fastq -S screened_out.sam

# Align to mitochondrial genome
#bowtie2 -p 10 -D20 -R 10 -N 1 -L 20 -i C,1 -x HumanMT -U screened.fastq -Q ALIGNMITO.sam

# Filter out multi mapping reads, convert to a sorted BAM file
#samtools view -q 1 -bS screened_out.sam -u|samtools sort - -f screened.sorted.bam

