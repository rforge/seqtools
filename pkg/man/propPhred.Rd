\name{propPhred}
%
\alias{propPhred}
\alias{propPhred-methods}
\alias{propPhred,Fastqq-method}
%
\title{
propPhred: Lane specific proportion of reads in a specified Phred-region.
}
\description{
The \code{propPhred} function returns a named vector with relative Phred content for all contained lanes. }
\usage{
propPhred(object, greater = 30, less = 93)
}
\arguments{
  \item{object}{\code{Fastqq}: Object which contains collected values from fastq files.}
  \item{greater}{\code{numeric}: Limits the counted proportion of phred to values which are greater than the given value (default: 30).}
  \item{less}{\code{numeric}: Limits the counted proportion of phred to values which are less than the given value (default: 93).}
}
\details{The \code{greater} and \code{less} arguments must be numeric, have length 1 and be >0 and < 94. \code{greater} must be less than \code{less}. With the default settings the reported proportions should be >50 \% for all lanes in order to be acceptable (see 't Hoen et. al.).}
\value{Numeric.}
\references{
Ewing B, Green P
Base-calling of automated sequencer traces using phred. II. Error probabilities.
Genome Research 1998 Vol. 8 No. 3 186-194
}
\author{
Wolfgang Kaisers
}
\examples{
basedir<-system.file("extdata",package="seqTools")
setwd(basedir)
fq<-fastqq(c("g4_l101_n100.fq.gz","g5_l101_n100.fq.gz"),k=4,probeLabel=c("g4","g5"))
# Proportion of phred Values >30
propPhred(fq)
# Proportion of phred Values >10 and < 30
propPhred(fq,greater = 10, less = 30)
}
\keyword{propPhred}