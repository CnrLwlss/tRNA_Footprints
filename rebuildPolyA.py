import subprocess
from itertools import islice

def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return(reverse_complement)

for root in ['ERR2208504','ERR2208505','ERR2208506','ERR2208507','ERR2208508','SRR935452','SRR935453']:
    print(root)
    sam = root + "_polyA_header.sam"
    fastqinfname = root + ".fastq"

    readidfname = root +".ids"
    fastqoutfname = root+"_OUT.fa"
    samoutfname = sam.split(".sam")[0]+"_fullreads.sam"

    readsout = open(readidfname,mode="w")

    with open(sam,mode="r") as samfile:
        for line in samfile:
            if line[0] != "@":
                fields = line.rstrip().split("\t")
                readsout.write(fields[0]+"\n")
                
    readsout.close()

    command = "seqtk subseq {} {} > {}".format(fastqinfname,readidfname,fastqoutfname)
    print(command)
    process = subprocess.Popen(command,shell=True)
    process.wait()

    readict = {}
    counter = 0
    with open(fastqoutfname,"r") as reads:
        for line in reads:
            if counter == 0:
                readid = line.rstrip()[1:]
            if counter == 1:
                readval = line.rstrip()
            if counter == 3:
                qualval = line.rstrip()
                res = {"seq":readval, "qual":qualval}
                readict[readid] = res
                counter = -1
            counter += 1

    samout = open(samoutfname,mode="w")

    with open(sam,mode="r") as samfile:
        for line in samfile:
            if line[0] == "@":
                samout.write(line)
            else:
                fields = line.rstrip().split("\t")
                seq = readict[fields[0]]["seq"]
                qual = readict[fields[0]]["qual"]
                if fields[1] == '16': # Aligned to reverse strand
                    seq = rev_comp(seq)
                    qual = qual[::-1]
                fields[9] =  seq
                fields[10] = qual
                samout.write("\t".join(fields)+"\n")
    samout.close()


            
