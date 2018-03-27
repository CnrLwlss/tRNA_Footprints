import subprocess
from itertools import islice
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return(reverse_complement)

def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts

def rebuild(root):
    print("Recovering full polyA tail reads... "+root)
    sam = root + "_unsorted.sam"
    fastqinfname = root + "_trimmed.fastq"

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
                readid = fields[0]
                seq = readict[readid]["seq"]
                qual = readict[readid]["qual"]
                offset = 0
                if fields[1] == '16': # Aligned to reverse strand
                    seq = rev_comp(seq)
                    qual = qual[::-1]
                    offset = -1*(len(seq)-len(fields[9]))
                newstart = int(fields[3]) + offset
                if newstart>0:
                    fields[3] = str(newstart)
                    fields[9] =  seq
                    fields[10] = qual
                    fields[5] = str(len(qual))+"M" # Assume all matches...
                else:
                    print("Warning, read runs off end of region.  Not adding polyA: "+fields[0])
                samout.write("\t".join(fields)+"\n")
    samout.close()

def justpolyA(root):
    print("Keeping reads with polyA tails... "+root)
    fqin = root+"_trimmed.fastq"
    fqout = root+"_filtered.fastq"

    fq = open(fqout,"w")
    counter = 0
    with open(fqin,"r") as reads:
        for line in reads:
            if counter == 0:
                head = line.rstrip()
                readid = head.split(" ")[0]
            elif counter == 1:
                seq = line.rstrip()
            elif counter == 2:
                conc = line.rstrip()
            elif counter == 3:
                qualval = line.rstrip()
                if seq.endswith('A'*5) and not seq.endswith('A'*26):
                    isA = [x is "A" for x in seq]
                    if False in isA:
                        last_nonA = len(isA)-isA[::-1].index(False)
                    else:
                        last_nonA = len(isA)
                    fq.write(head+"\n")
                    fq.write(seq[0:last_nonA]+"\n")
                    fq.write(conc+"\n")
                    fq.write(qualval[0:last_nonA]+"\n")
                counter = -1
            counter += 1

def parseattr(attr):
    at = attr.split(";")
    atdict={}
    for a in at:
        key,val=a.split("=")
        atdict[key]=val
    return(atdict)

def getattr(attr,key="Name"):
    atdict = parseattr(attr)
    return(atdict[key])

def getstop(gff):
    if gff.Strand=="+":
        firstcoord = gff.End-2
    else:
        firstcoord = gff.Start
    return(firstcoord)

if __name__ == "__main__":
    from sys import argv
    myargs = getopts(argv)
    if '-r' in myargs:
        root = myargs['-r']
        rebuild(root)
    if '-f' in myargs:
        root = myargs['-f']
        justpolyA(root)    
    if '-a' in myargs:
        roots = ['ERR2208504','ERR2208505','ERR2208506','ERR2208507','ERR2208508','SRR935452','SRR935453']
        for root in roots:
            justpolyA(root)
            rebuild(root)
            
    roots = ['ERR2208504','ERR2208505','ERR2208506','ERR2208507','ERR2208508','SRR935452','SRR935453']
    gff3 = 'mtDNA.gff3'
    pus = {root:pd.read_csv(root+'_polyA.pileup',sep="\t",header=None,names=["Genome","Coordinate","RefBase","Coverage","ReadQual","AlignQual","Start","End"]) for root in roots}
    gf = pd.read_csv(gff3,skiprows=2,sep="\t",header=None,names=["Sequence","Source","Feature","Start","End","Score","Strand","Phase","Attributes"])
    N = int(gf.End[gf.Feature=="region"])
    coverage = pd.DataFrame(columns=roots)
    
    for root in roots:
        covs = [0 for i in range(0,N-1)]
        for i,g in pus[root].iterrows():
            covs[g.Coordinate]=g.Coverage
        coverage[root]=covs

    totals={}
    for root in roots:
        totals[root] = sum(1 for line in open(root+'_polyA_header.sam'))-3

    gfg = gf[gf.Feature=='gene']
    gfg["Gname"] = [getattr(attr) for attr in gfg.Attributes]
    gfg["FirstStop"] = [getstop(g) for i,g in gfg.iterrows()]

    summary = {}
    for i,g in gfg.iterrows():
        fs = int(g.FirstStop)
        summary[g.Gname] = coverage[fs:(fs+3)].median()

    res = pd.DataFrame(summary)
    perc = 100*res.divide(pd.Series(totals),axis=0)
    res["Total"]=pd.Series(totals)
    res.to_csv("polyA_Counts.txt",sep="\t")
    perc.to_csv("polyA_Percentages.txt",sep="\t")


    






