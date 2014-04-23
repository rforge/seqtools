### R code from vignette source 'seqTools.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: seqTools.Rnw:40-41
###################################################
options(width=60)


###################################################
### code chunk number 2: seqTools.Rnw:78-80
###################################################
library(seqTools)
head(getPhredTable())


###################################################
### code chunk number 3: computation
###################################################
basedir<-system.file("extdata",package="seqTools")
filenames<-file.path(basedir,c("g4_l101_n100.fq.gz","g5_l101_n100.fq.gz"))
fq<-fastqq(filenames,k=6,probeLabel=c("g4","g5"))


###################################################
### code chunk number 4: seqTools.Rnw:116-118 (eval = FALSE)
###################################################
## filenames<-dir(path="fastqDir",pattern="*.fastq.gz")
## fq<-fastqq(file.path("fastqDir",filenames),k=6)


###################################################
### code chunk number 5: seqTools.Rnw:123-124
###################################################
fq


###################################################
### code chunk number 6: seqTools.Rnw:191-193
###################################################
fqm<-meltDownK(fq,newK=2)
kmerCount(fqm)[,1]


###################################################
### code chunk number 7: computation
###################################################
files1<-file.path(basedir,c("sfq1_ctrl.fq.gz","sfq2_ctrl.fq.gz"))
files2<-file.path(basedir,c("sfq1_cont.fq.gz","sfq2_cont.fq.gz"))
fq1<-fastqq(files1,k=3,probeLabel=c("ctrl1","ctrl2"))
fq2<-fastqq(files2,k=3,probeLabel=c("cont1","cont2"))


###################################################
### code chunk number 8: seqTools.Rnw:208-212
###################################################
op<-par(mfrow=c(1,2))
plotKmerCount(fq1,2,mxey=9,main="Control")
plotKmerCount(fq2,2,mxey=9,main="Contamination")
par(op)


###################################################
### code chunk number 9: seqTools.Rnw:216-218
###################################################
mrg<-mergeFastqq(fq1,fq2)
mrg


###################################################
### code chunk number 10: seqTools.Rnw:229-230
###################################################
plotNucFreq(fq,1)


###################################################
### code chunk number 11: seqTools.Rnw:234-235
###################################################
plotNucCount(fq)


###################################################
### code chunk number 12: seqTools.Rnw:239-240
###################################################
plotNucCount(fq,c(2,3))


###################################################
### code chunk number 13: seqTools.Rnw:245-246
###################################################
plotGCcontent(fq)


###################################################
### code chunk number 14: seqTools.Rnw:250-251
###################################################
plotPhredQuant(fq,1,"Phred quantiles for 1st file")


###################################################
### code chunk number 15: seqTools.Rnw:256-257
###################################################
plotMergedPhredQuant(fq,main="Phred quantiles for all files")


###################################################
### code chunk number 16: seqTools.Rnw:263-266
###################################################
phred<-phredDist(fq,1)
phred<-phredDist(fq)
head(phred)


###################################################
### code chunk number 17: seqTools.Rnw:271-272
###################################################
plotPhredDist(fq)


###################################################
### code chunk number 18: seqTools.Rnw:280-281
###################################################
getKmerIndex(c("CCC","GGG"))


###################################################
### code chunk number 19: seqTools.Rnw:285-286
###################################################
plotKmerCount(fq,1)


###################################################
### code chunk number 20: seqTools.Rnw:296-298
###################################################
mtx<-cbDistMatrix(mrg)
mtx


###################################################
### code chunk number 21: seqTools.Rnw:302-307
###################################################
hc<-hclust(as.dist(mtx))
hcd<-as.dendrogram(hc,lty=2,lwd=2)
op<-par(mar=c(3,1,1,5))
plot(hcd,horiz=TRUE,las=1,edgePar=list(lwd=2,lty=2,col="blue"))
par(op)


###################################################
### code chunk number 22: seqTools.Rnw:316-318
###################################################
char2ascii("a")
ascii2char(97:99)


###################################################
### code chunk number 23: seqTools.Rnw:322-324 (eval = FALSE)
###################################################
## getPhredTable()
## getPhredTable(20:30)


