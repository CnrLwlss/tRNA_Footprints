library(corrplot)
source("cplt.R")

# Correlation plot for examining similarities between profiles
corscat=function(colx,coly,df,nsig=3){
  C=cor(df[colx],df[coly])
  plot(df[[colx]],df[[coly]],pch=16,cex=0.5,main=paste("Correlation: ",formatC(C,nsig)),xlab=colx,ylab=coly,col="grey")
  text(x=df[[colx]],y=df[[coly]],labels=rownames(df),pos=1,cex=0.4,srt=0,offset=0.1)
}

# Draw profiles, highlighting amino acids
profplot=function(countlist,codons,aminodf=c(),main="",ylab="Abundance",lwd=2,type="h",cols=c(),ltys=c(),legpos="",perfcol="",crosses=c(),pluses=c(),crosscol="black",pluscol="black"){
  maxcount=max(unlist(countlist),na.rm=TRUE)
  if(length(cols)==0) cols=rainbow(length(countlist),start=0.7,end=0.1)
  if(length(ltys)==0) ltys=rep(1,length(cols))
  plot(1:length(countlist[[1]]),countlist[[1]],xlab="",ylab=ylab,main=main,type="n",xaxt="n",ylim=c(0,1.05*maxcount))
  codons$Colour="black"
  if(perfcol!="") codons$Colour[codons$Perfect]=perfcol
  Map(function(x,y,z)
    axis(1,at=x,col.axis=y,labels=z,lwd=0,las=2),
    codons$CodonNum,
    codons$Colour,
    codons$Codon
  )
  if(length(crosses)!=0) axis(1,at=codons$CodonNum,labels=ifelse(crosses,rawToChar(as.raw(215)),""),line=1.75,lwd=0,cex=1,col.axis=crosscol)
  if(length(pluses)!=0) axis(1,at=codons$CodonNum,labels=ifelse(pluses,"+",""),line=1.75,lwd=0,cex=1,col.axis=pluscol)
  axis(1,at=codons$CodonNum,labels=FALSE)
  if(length(aminodf)!=0){
    rect(aminodf$Lower,0,aminodf$Upper,1.05*maxcount,col=aminodf$Colour,border=aminodf$Colour)
    text(aminodf$Centre,maxcount,aminodf$AA)
  }
  for(i in 1:length(countlist)){
    cname=names(countlist)[i]
    counts=countlist[[i]]
    points(1:length(counts),counts,type=type,lwd=lwd,col=cols[i],lty=ltys[i])
  }
  if((legpos!="") & (!is.null(names(countlist)))) legend(x=legpos,legend=names(countlist),col=cols,lwd=lwd,lty=ltys,bg="white")
}

visualiseCorrelation = function(df,label=""){
  M=cor(df)
  D=as.dist(1-M)
  opp=par(mfrow=c(1,2))
  plot(hclust(D),xlab="",main=paste(label,"dendogram (1-Correlation)"),sub="")
  cplt(M,method="number",order="hclust",addrect=3,nsig=3,main=paste(label,"correlation matrix"))
  par(opp)
}

#################

gnames=c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6")

pdf("FootprintReport.pdf",width=20,height=8)

#for(fgene in c("",gnames)){
#for(fgene in c(gnames)){
for(fgene in c("")){

# Define how codons code for amino acids, including perfect codon-anticodon matches
codons = data.frame(
AA = c("K", "K", "N", "N", "T", "T", "T", "T",
"M", "M", "I", "I", "Q", "Q", "H", "H", "P", "P", "P", "P", "R",
"R", "R", "R", "L", "L", "L", "L", "L", "L", "E", "E", "D", "D",
"A", "A", "A", "A", "G", "G", "G", "G", "V", "V", "V", "V", "Y",
"Y", "S", "S", "S", "S", "S", "S", "W", "W", "C", "C", "F", "F",
"stop", "stop", "stop", "stop"),
Codon = c("AAA", "AAG", "AAC",
"AAT", "ACA", "ACC", "ACG", "ACT", "ATA", "ATG", "ATC", "ATT",
"CAA", "CAG", "CAC", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA",
"CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "TTA", "TTG",
"GAA", "GAG", "GAC", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA",
"GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT",
"AGC", "AGT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGG", "TGC",
"TGT", "TTC", "TTT", "AGA", "AGG", "TAA", "TAG"),
Perfect = c(T,F,F,F,T,F,F,F,F,T,T,F,T,F,F,F,T,F,F,F,T,F,F,F,T,F,F,F,F,F,T,F,F,F,T,F,F,F,T,F,F,F,T,F,F,F,F,F,T,F,T,F,F,F,T,F,T,F,T,F,F,F,F,F),
stringsAsFactors=FALSE)
codons$CodonNum=1:length(codons$Codon)

#################

# Read in data extracted from Fei's spreadsheet
clines = read.delim("cell_lines.txt",sep="\t",stringsAsFactors=FALSE)
# Discard untranslated regions of genes
clines = clines[!is.na(clines$frame),]
# Name samples systematically
samples = c("Val_1","Val_2","Val_3","Phe_1","Phe_2","Phe_3")
names(samples) = c("SRR935452.sorted", "SRR935453.sorted", "X2934_2_RP2_Control.sorted",
"X2934_1_RP2_.Tet.sorted", "Tet.sorted", "Tet_Plus.sorted")
colnames(clines)[colnames(clines)%in%names(samples)] = samples[colnames(clines)[colnames(clines)%in%names(samples)]]

# Read in some metadata about the sequence, from Fei
meta = read.delim("mitoreference.csv",sep=",",stringsAsFactors=FALSE)

# FILTER BY GENE
if(fgene!=""){
  clines=clines[clines$gene==fgene,]
  meta=meta[meta$gene==fgene,]
  codons=codons[codons$Codon%in%clines$aa_codon,]
  codons$CodonNum=1:length(codons$Codon)
}

# Discard codons which are in untranslated regions of genes
meta = meta[!is.na(meta$frame),]
# Count number of occurances of each codon and number of codons in each gene
abundances = table(meta$aa_codon[meta$frame==1])
codons$Count = abundances[codons$Codon]
codons$Count[is.na(codons$Count)] = 0
gabundances = table(meta$gene[meta$frame==1])

# Generate set of locations on plot for amino acids
aminos=unique(codons$AA)
aminodf=data.frame(
AA=aminos,
Centre=tapply(codons$CodonNum, codons$AA, mean)[aminos],
Lower=tapply(codons$CodonNum, codons$AA, min)[aminos]-0.5,
Upper=tapply(codons$CodonNum, codons$AA, max)[aminos]+0.5,
Colour=rep(c("white","lightgrey"),length(aminos))[1:length(aminos)],
stringsAsFactors=FALSE
)

# Calculate footprint abundances in coding regions by codon
ftcodon = aggregate(clines[,samples],by=list(Codon=clines$aa_codon),sum)
ftcodon$mtDNA = aggregate(clines[,"aa_codon"],by=list(Codon=clines$aa_codon),length)$x/3
ftcodon = ftcodon[!is.na(ftcodon[,1]),]

ftcodon = ftcodon[match(codons$Codon,ftcodon$Codon),]
ftcodon$mtDNA = codons$Count

# Calculate footprint abundances in coding regions by gene
ftgene = aggregate(clines[,samples],by=list(Gene=clines$gene),mean,na.rm=TRUE)
ftgene$mtDNA = aggregate(clines[,"aa_codon"],by=list(Gene=clines$gene),length)$x
ftgene = ftgene[!is.na(ftgene[,1]),]

# Fractional abundances for comparison between mtDNA and samples
ftfcodon=ftcodon
ftfcodon[,c(samples,"mtDNA")]=prop.table(as.matrix(ftfcodon[,c(samples,"mtDNA")]),2)
ftfgene=ftgene
ftfgene[,c(samples,"mtDNA")]=prop.table(as.matrix(ftfgene[,c(samples,"mtDNA")]),2)

# Look at gene profiles for each cell line, comparing with abdundance in mtDNA
genes=data.frame(Codon=unique(ftgene$Gene),Perfect=FALSE,CodonNum=1:length(unique(ftgene$Gene)))

# Compare with Fei's average plots
Vals=c("Val_1","Val_2","Val_3")
Phes=c("Phe_1","Phe_2","Phe_3")
f_ave_frac=data.frame(Val=apply(ftfcodon[,Vals],1,mean), Phe=apply(ftfcodon[,Phes],1,mean))
f_ave_frac$mtDNA=ftfcodon$mtDNA
f_sd_frac=data.frame(Val=apply(ftfcodon[,Vals],1,sd), Phe=apply(ftfcodon[,Phes],1,sd))

mslp=1
mint=0
f_frac=ftfcodon
f_frac$pred=f_frac$mtDNA*mslp+mint
f_frac$Val_diff=f_frac$pred-apply(f_frac[,c("Val_1","Val_2","Val_3")],1,mean)
f_frac$Phe_diff=f_frac$pred-apply(f_frac[,c("Phe_1","Phe_2","Phe_3")],1,mean)
f_frac$Val_p=1
f_frac$Phe_p=1
for(i in 1:length(f_frac$pred)){
  f_frac$Val_p[i]=t.test(as.numeric(f_frac[i,c("Val_1","Val_2","Val_3")]),mu=f_frac$pred[i])$p.value
  f_frac$Phe_p[i]=t.test(as.numeric(f_frac[i,c("Phe_1","Phe_2","Phe_3")]),mu=f_frac$pred[i])$p.value
  f_frac$Val_Phe_p[i]=t.test(f_frac[i,c("Val_1","Val_2","Val_3")],f_frac[i,c("Phe_1","Phe_2","Phe_3")])$p.value
  #f_frac$Val_p[i]=wilcox.test(as.numeric(f_frac[i,c("Val_1","Val_2","Val_3")]),mu=f_frac$pred[i])$p.value
  #f_frac$Phe_p[i]=wilcox.test(as.numeric(f_frac[i,c("Phe_1","Phe_2","Phe_3")]),mu=f_frac$pred[i])$p.value
}
f_frac$Val_q=p.adjust(f_frac$Val_p,method="fdr")
f_frac$Phe_q=p.adjust(f_frac$Phe_p,method="fdr")
f_frac$Val_Phe_q=p.adjust(f_frac$Val_Phe_p,method="fdr")
f_frac$Val_sign=ifelse(sign(f_frac$Val_diff)>0,"+","-")
f_frac$Phe_sign=ifelse(sign(f_frac$Phe_diff)>0,"+","-")
f_frac$Val_p_hits=ifelse(f_frac$Val_p<0.05,"*","")
f_frac$Phe_p_hits=ifelse(f_frac$Phe_p<0.05,"*","")
f_frac$Val_Phe_p_hits=ifelse(f_frac$Val_Phe_p<0.05,"*","")
f_frac$Val_q_hits=ifelse(f_frac$Val_q<0.05,"*","")
f_frac$Phe_q_hits=ifelse(f_frac$Phe_q<0.05,"*","")
f_frac$Val_Phe_q_hits=ifelse(f_frac$Val_Phe_q<0.05,"*","")
f_frac$Matched=ifelse(codons$Perfect,"*","")
f_frac$AA=codons$AA
write.table(file=paste(fgene,"FractionalFootprints.txt",sep="_"),f_frac,sep="\t",quote=FALSE,row.names=FALSE)

# Examine abundance of codons in genome
profplot(countlist=list(mtDNA=codons$Count),aminodf=aminodf,codons=codons,main=paste("Coding mtDNA",fgene),perfcol="red",lwd=3)
if(fgene=="") profplot(countlist=list(Gene=ftgene$mtDNA),aminodf=c(),codons=genes,main=paste("Coding mtDNA",fgene),perfcol="red",lwd=3,ylab="Gene length (b)")

# Look at codon profiles for each cell line, comparing with abundance in mtDNA
profplot(
countlist=as.list(ftfcodon[,c(samples,"mtDNA")]),
aminodf=aminodf,
codons=codons,
main=paste("Ribosome footprint & codon abundance in coding regions of mtDNA (p-value)",fgene),
ylab="Fractional abundance",
type="l",
cols=c("red","red","red","blue","blue","blue","green"),
ltys=c(1,2,3,1,2,3,1),
legpos="right",
perfcol="red",
crosses=f_frac$Val_p_hits=="*",
pluses=f_frac$Phe_p_hits=="*",
crosscol="red",
pluscol="blue"
)

# Look at codon profiles for each cell line, comparing with abundance in mtDNA
profplot(
countlist=as.list(ftfcodon[,c(samples,"mtDNA")]),
aminodf=aminodf,
codons=codons,
main=paste("Ribosome footprint & codon abundance in coding regions of mtDNA (q-value)",fgene),
ylab="Fractional abundance",
type="l",
cols=c("red","red","red","blue","blue","blue","green"),
ltys=c(1,2,3,1,2,3,1),
legpos="right",
perfcol="red",
crosses=f_frac$Val_q_hits=="*",
pluses=f_frac$Phe_q_hits=="*",
crosscol="red",
pluscol="blue"
)

# Look at codon profiles for each cell line, comparing Val with Phe
profplot(
countlist=as.list(ftfcodon[,samples]),
aminodf=aminodf,
codons=codons,
main=paste("Comparing ribosome footprints in coding regions of mtDNA (p-value)",fgene),
ylab="Fractional abundance",
type="l",
cols=c("red","red","red","blue","blue","blue"),
ltys=c(1,2,3,1,2,3),
legpos="right",
perfcol="red",
pluses=f_frac$Val_Phe_p_hits=="*",
pluscol="black"
)

# Look at codon profiles for each cell line, comparing Val with Phe
profplot(
countlist=as.list(ftfcodon[,samples]),
aminodf=aminodf,
codons=codons,
main=paste("Comparing ribosome footprints in coding regions of mtDNA (q-value)",fgene),
ylab="Fractional abundance",
type="l",
cols=c("red","red","red","blue","blue","blue"),
ltys=c(1,2,3,1,2,3),
legpos="right",
perfcol="red",
pluses=f_frac$Val_Phe_q_hits=="*",
pluscol="black"
)

# Look at gene profiles for each cell line, comparing with abundance in mtDNA
if(fgene=="") {
 profplot(
 countlist=as.list(ftfgene[,c(samples,"mtDNA")]),
 aminodf=c(),
 codons=genes,
 main="Ribosome footprint & gene length in coding regions of mtDNA",
 ylab="Fractional abundance/length",
 type="l",
 cols=c("red","red","red","blue","blue","blue","green"),
 ltys=c(1,2,3,1,2,3,1),
 legpos="right",
 perfcol="red"
 )
}


# Check whether cell lines and mtDNA cluster together after unsupervised hierarchical clustering
visualiseCorrelation(ftcodon[,c(samples,"mtDNA")],label=paste("Codon",fgene))
if(fgene=="") visualiseCorrelation(ftgene[,c(samples,"mtDNA")],label="Gene")

profplot(
countlist=as.list(f_ave_frac),
aminodf=aminodf,
codons=codons,
main=paste("Ribosome footprint & codon abundance in coding regions of mtDNA",fgene),
legpos="right",
type="l",
cols=c("red","blue","green"),
ylab="Mean fractional abundance",
perfcol="red"
)

# Linear regression tests:
# Is relationship between mtDNA abundance and mitoribosome occupancy well represented by a straight line?  Yes, high R2, low p-values
# Do slopes go through origins (simple proportionality)?  Often no (p-value for intercept > 0.05).
# However, mtDNA occupancy observations are not evenly spaced, seems likely that fits are affected by spacing: just use proportions
op=par(mfrow=c(2,3))
for(dep in samples){
  corr=cor(ftcodon$mtDNA,ftcodon[[dep]])
  regm=lm(get(dep)~mtDNA+0,data=ftcodon)
  reg=lm(get(dep)~mtDNA+1,data=ftcodon)
  plot(ftcodon$mtDNA,ftcodon[[dep]],xlab=paste("Codon abundance in coding mtDNA",fgene),ylab=paste(dep,"codon footprint reads",fgene),main=paste("Pearson's Correlation:",formatC(corr,4)))
  abline(regm,col="blue",lwd=2)
  abline(reg,col="red",lwd=2)
  print(dep)
  print(summary(reg))
}
legend("bottomright",c("Linear regression","Regression through origin"),col=c("red","blue"),lwd=2)
par(op)

if(fgene=="") {
 op=par(mfrow=c(2,3))
 for(dep in samples){
   corr=cor(ftgene$mtDNA,ftgene[[dep]])
   regm=lm(get(dep)~mtDNA+0,data=ftgene)
   reg=lm(get(dep)~mtDNA+1,data=ftgene)
   plot(ftgene$mtDNA,ftgene[[dep]],xlab="Gene abundance in coding mtDNA",ylab=paste(dep,"gene footprint reads"),main=paste("Pearson's Correlation:",formatC(corr,4)))
   abline(regm,col="blue",lwd=2)
   abline(reg,col="red",lwd=2)
   print(dep)
   print(summary(reg))
 }
 legend("bottomright",c("Linear regression","Regression through origin"),col=c("red","blue"),lwd=2)
 par(op)
}

if(fgene!=""){
 cols=c("red","red","red","blue","blue","blue")
 ltys=c(1,2,3,1,2,3)
 fclines=data.frame(prop.table(as.matrix(clines[,samples]),2),stringsAsFactors=FALSE)
 fclines$pos=clines$pos
 if(length(cols)==0) cols=rainbow(length(countlist),start=0.7,end=0.1)
 if(length(ltys)==0) ltys=rep(1,length(cols))
 plot(fclines$pos,fclines[,samples[1]],type="n",ylim=c(0,0.02),xlab="Chromosome coordinate",ylab="Fractional mitoribosome residency",main=fgene)
 for(i in seq_along(samples)){
   points(fclines$pos,fclines[,samples[i]],type="l",col=cols[i],lty=ltys[i],lwd=2)
 }
 legend("topleft",samples,col=cols,lty=ltys,lwd=2)
}

}

dev.off()
