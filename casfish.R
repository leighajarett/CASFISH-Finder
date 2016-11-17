#CASFISHFind

library(data.table)
library(plyr)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDB.Hsapiens.UCSC.hg19.knownGene)
library(CRISPRseek)
library(org.Hs.eg.db)

#create a new directory for output files
dir.create("casfish")
setwd("casfish")

#user input - chrom, start, end

#load target range
target = subseq(BSgenome.Hsapiens.UCSC.hg19[[chrom]],start = start, end = end)

#find all PAM sequences
p1 = "GG" 
pams = matchPattern(p1,target)

#save a vector of all candidate 20bp guide sequences 
starts = start(pams)
guides - c()
for (i in 1:length(starts)){
	if (starts[i] > 21){
		guides[i] <- toString(subseq(target,start = (starts[i]-21),width=20))
	}
	else{
		guides[i] <- NA
	}
} 
guides <- na.omit(guides)

#create a fasta file of candidate guides
file.create("query.fa")
fileConn <- file("query.fa",'w')
qnames<-c()
for (i in 1:length(guides)){
	qnames[i] <- sprintf(">Query_%s",i)
	write(c(qnames[i],guides[i],fileConn, append = TRUE))
}
close(fileConn)

outputDir <- getwd()
guidesFilePath <- sprintf('%s/query.fa',outputDir)
REpatternFile <- system.file('extdata', 'NEBenzymes.fa', pacakge = 'CRISPRseek')

#perform off target analysis
results <- offTargetAnalysis(inputFilePath = guidesFilePath, findgRNAsWithREcutOnly = FALSE, REpatternFile = REpatternFile, findPairedgRNAOnly = FALSE, findgRNAs = FALSE, BSgenomeName = Hsapiens, chromToSearch = 'all', txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, orgAnn = org.Hs.egSYMBOL, max.mismatch = 4, outputDir = outputDir, overwrite = TRUE)

#determine off target to on target ratio 
Offtarget = read.csv('/Users/leighajarett/Offtarget.csv')
name = Offtarget$name
guide = Offtarget$gRNAPlusPAM
score = Offtarget$score
chrom = Offtarget$chrom
g_start = Offtarget$chromStart
g_end = Offtarget$chromEnd

data = data.table(name, guide, score, chrom, g_start, g_end)

on_tar = c()
off_tar = c()

on_scores = c()
off_scores = c()
gRNA =c()
names = c(levels(data$name))
for (i in names){
	rows = c(which(data$name == i))
	for (j in rows){
		score = data$score[j]
		if ((data$chrom[j] == input$chrom)&(input$start <= data$g_start[j])&(input$end >= data$g_end[j])){
			on_tar = append(on_tar,score)
			off_tar = append(off_tar,0)}
		else {
			off_tar = append(off_tar,score)
			on_tar =  append(on_tar, 0)}}
	on_tar_score = sum(on_tar)
	off_tar_score = sum(off_tar)
	on_scores = append(on_scores, on_tar_score)
	off_scores = append(off_scores, off_tar_score)
	gRNA = append(gRNA, i) 
}

output = data.table(gRNA, on_scores, off_scores)


		
		




