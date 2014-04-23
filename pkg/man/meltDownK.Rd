\name{meltDownK}
%
\alias{meltDownK}
\alias{meltDownK-methods}
\alias{meltDownK,Fastqq-method}
%
\title{
meltDownK: Condensing DNA k-mer count data to lower K (i.e. shorter DNA motifs).
}
\description{
Returns a copy of given object where DNA k-mer counts and first DNA k-mer count table are reduced in size. }
\usage{
meltDownK(object,newK)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{Fastqq}: Object from wich data is copied.}
  \item{newK}{\code{integer}: New value for k. Must be >=1 and <= old k.}
}
\details{The function sums all count values which belong the the new motif up. The new motif is the new-k sized prefix of the given k-mer motif.}
\value{S4 Object of class 'Fastqq'.}
\references{
Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM
The sanger fastq file format for sequences with quality scores and the Solexa/Illumina fastq variants.
Nucleic Acids Research 2010 Vol.38 No.6 1767-1771
}
\author{Wolfgang Kaisers}
\note{The static size of the retured k-mer array is 4^k.}
\examples{
basedir<-system.file("extdata",package="seqTools")
setwd(basedir)
fq<-fastqq(c("g4_l101_n100.fq.gz","g5_l101_n100.fq.gz"),k=4,probeLabel=c("g4","g5"))
fqm<-meltDownK(fq,2)
}
\keyword{meltDownK}
\keyword{kmer}
\keyword{fastq}
\keyword{Fastqq}