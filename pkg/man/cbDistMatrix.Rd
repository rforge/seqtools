\name{cbDistMatrix}
%
\alias{cbDistMatrix}
\alias{cbDistMatrix-methods}
\alias{cbDistMatrix,Fastqq-method}
%

\title{
cbDistMatrix function: Calculates pairwise distance matrix from DNA k-mer counts based on a modified Canberra distance.
}
\description{
Calculates pairwise distance matrix from DNA k-mer counts based on a modified Canberra distance. Before calculating canberra distances, sequencing depth artifacts are compensated by saling up read counts in each DNA k-mer count vector so that scaled read counts in each sample are nearly equal.
}
\usage{
cbDistMatrix(object,nReadNorm=max(nReads(object)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{Fastqq}: Object from which DNA k-mer counts are used.}
  \item{nReadNorm}{\code{numeric}: Number of reads per file to wich all contained DNA k-mer counts are scaled. Because it is intended to up-scale the value must be greater than all fastq file read counts (as reported by nReads). Therefore the standard value is chosen to the maximal number of reads recorded in this object. This scaling is necessary to compensate for systematic effects in the canberra distance.}
}
\details{The DNA k-mer counts are shrinked to a k of 6. The distance between two DNA k-mer count vectors is calculated by
\deqn{df_0(X,Y)=\frac{\sum_{i=1}^n cbd(x_i,y_i)}{4^k}}
where cb is given by
\deqn{cbd(x,y)=\frac{|x-y|}{x+y}}.
}
\value{Square \code{matrix}. The number of rows equals the number of files (=\code{nFiles(object)}). }
\references{
Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM
The sanger fastq file format for sequences with quality scores and the Solexa/Illumina fastq variants.
Nucleic Acids Research 2010 Vol.38 No.6 1767-1771
}
\author{
Wolfgang Kaisers
}
\note{The static size of the retured k-mer array is 4^k.}
\examples{
basedir<-system.file("extdata",package="seqTools")
basenames<-c("g4_l101_n100.fq.gz","g5_l101_n100.fq.gz")
filenames<-file.path(basedir,basenames)
fq<-fastqq(filenames,6,c("g4","g5"))
dm<-cbDistMatrix(fq,1e3)
}
\keyword{cbDistMatrix}
\keyword{kmer}
\seealso{hclust}
