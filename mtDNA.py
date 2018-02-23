from Bio import SeqIO
from Bio import SeqFeature
fasta = SeqIO.parse("mtDNA.fa","fasta")

mito_stop=["TAA","TAG","AGA","AGG"]

def nextstop(seq,stops):
    triplets = [seq[i:min(i+3,len(seq))] for i in range(0, len(seq), 3)]
    ind = 0
    for i,trip in enumerate(triplets):
        if trip in mito_stop:
            return(3*i)
    return(0)

for seq_record in fasta:
    print(seq_record.id)
    mtDNA = seq_record.seq
    print(len(seq_record))

triplets = [str(mtDNA)[i:i+3] for i in range(0, len(mtDNA), 3)]
is_stop = ''.join(['111' if t in mito_stop else '000' for t in triplets])

starts = [i for i in range(0,len(str(mtDNA))-3) if str(mtDNA)[i:(i+3)] == "ATG"]
stops = [i+nextstop(str(mtDNA)[i:len(mtDNA)],mito_stop) for i in starts]
stops = list(set(stops))

fout = open("mtDNA_stop.gff3","w")
with open("mtDNA.gff3","r") as fin:
    for i,line in enumerate(fin):
        if i<2:
            fout.write(line)
        else:
            if (("region" not in line) and ("sequence_alteration" not in line)):
                fout.write(line)
ind = 1
for i,s in enumerate(stops):
    i0,iN = s+1, s+3
    fout.write('NC_012920.1	BioPython_script	stop_codon	{}	{}	.	+	.	ID=sc{};Name=SC\n'.format(str(i0),str(iN),str(i)))
fout.close()
