\name{writeSimContFastq}
\alias{writeSimContFastq}
\title{writeSimContFastq: Create fastq files with simulated k-mer sequences}
\description{
Writes compressed fastq files where sequence sections contain concatenated k-mers which are uniformly distributed in the range of k-mers for given k. A fraction of the reads can be contaminated with one or more deterministic k-mers.}
\usage{
writeSimContFastq(k=6,nk=5,nSeq=10,pos=1,kIndex=1,nContam=nSeq,filename="simc.fq.gz")
}
\arguments{
  \item{k}{\code{numeric}. Length of k-mer. Default value is 6.}
  \item{nk}{\code{numeric}. Number of k-mers in each fastq read. Default value is 5.}
  \item{nSeq}{\code{numeric}. Number of simulated reads in fastq-file. Default value is 10.}
  \item{pos}{\code{numeric}. Determines at which position in sequence the k-mer is inserted. 1-based (1=first position).}
  \item{kIndex}{\code{numeric}. k-mer index of inserted k-mer. The k-mer index can be retreaved for a given k-mer with 'getKmerIndex'. }
  \item{nContam}{\code{numeric}. Absolute number of contaminated reads. The k-mer's are inserted at the firsts 'nContam' reads of the sequence array. }
  \item{filename}{\code{character}. Name of written (compressed) fastq file.}
}
\details{The read headers are consequtive numbered. The phred quality values are equally set to 46 (='.') which represents a phred value of 13.}
\value{None.}
\author{Wolfgang Kaisers}
\examples{
writeSimContFastq()
}
\keyword{writeSimContFastq}