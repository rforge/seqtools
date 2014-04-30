\name{writeSimFastq}
\alias{writeSimFastq}
\title{writeSimFastq: Create fastq files with simulated DNA k-mer sequences}
\description{
Writes compressed fastq files where sequence sections contain concatenated k-mers which are uniformly distributed in the range of k-mers for given k.}
\usage{
writeSimFastq(k=6,nk=5,nSeq=10,filename="sim.fq.gz")
}
\arguments{
  \item{k}{\code{numeric}. Length of k-mer. Default value is 6.}
  \item{nk}{\code{numeric}. Number of k-mers in each fastq read. Default value is 5.}
  \item{nSeq}{\code{numeric}. Number of simulated reads in fastq-file. Default value is 10.}
  \item{filename}{\code{character}. Name of written (compressed) fastq file.}
}
\details{The read headers are consequtive numbered. The phred quality values are equally set to 46 (='.') which represents a phred value of 13. This function is not designed for routine use. The random content fastq files can be used in order to measure the separation capabilities of hierarchical clustering mechanisms.}
\value{None.}
\references{
Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM
The sanger fastq file format for sequences with quality scores and the Solexa/Illumina fastq variants.
Nucleic Acids Research 2010 Vol.38 No.6 1767-1771
}
\author{Wolfgang Kaisers}
\examples{
writeSimFastq()
}
\keyword{writeSimFastq}