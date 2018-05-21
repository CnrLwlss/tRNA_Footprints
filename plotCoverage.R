coverage = list.files(".",pattern="*.coverage")
names(coverage) = sapply(coverage, function(x) strsplit(x,"_")[[1]][1])

gff3 = read.delim("mtDNA_stop.gff3",sep="\t",stringsAsFactors=FALSE,skip=2,header=FALSE)
colnames(gff3) = c("seqid","source","type","start","end","score","strand","phase","attributes")
genes = gff3[gff3$type=="gene",]

cols = c("black","black","red","red","red","blue","blue")
typs = c(1,2,1,2,3,1,2)
names(cols) = c("ERR2208504", "ERR2208505", "ERR2208506", "ERR2208507", "ERR2208508", 
"SRR935452", "SRR935453")
names(typs) = names(cols)
pattern=";Name=(.*?);"
genes$gene = sapply(genes$attributes,function(x) regmatches(x,regexec(pattern,x))[[1]][2])

chrom = "NC_012920.1"
N = 16569
counts = c(1393,111,155,122,470,4785,347)
names(counts) = c("ERR2208504", "ERR2208505", "ERR2208506", "ERR2208507", "ERR2208508", "SRR935452", "SRR935453")


df = data.frame(Chromosme = chrom, Coordinate = 1:N)
for (samp in names(coverage)) df[[samp]] = 0

for (c in names(coverage)){
 dfc = read.delim(coverage[c], sep="\t",stringsAsFactors=FALSE,header=FALSE)
 names(dfc) = c("Chromosome","Coordinate","Coverage")
 df[[c]][dfc$Coordinate] = dfc$Coverage
}

# Total number of aligned polyA reads per sample
sums = colSums(df[,3:length(names(df))])
print(sums)

df_frac = df
df_frac[,3:length(names(df))] = sweep(df[,3:length(names(df))],MARGIN=2,FUN="/",STATS=counts)

pdf("GeneCoverages.pdf",width=11.69,height=8.27,pointsize=16)

for(g in 1:length(genes$gene)){

  #g = which(genes$gene=="COX3")
  gene = genes$gene[g]
  start = genes$start[g]
  end = genes$end[g]
  dfg = df_frac[start:end,]
  maxfrac = max(dfg[,3:length(names(df))])

  plot(NULL,xlim=c(start,end),ylim=c(0,maxfrac),xlab="Chromosome coordinate",ylab="Coverage (% aligned poly(A) reads)",main=gene)
  for(c in names(coverage)){
    points(dfg$Coordinate,dfg[[c]],lwd=2,type="l",col = cols[c],lty = typs[c])
  }
  legend("top",col=c("black","red","blue"),legend=c("Control","Mito. disease","Cybrid"),lwd=2)
}

dev.off()




