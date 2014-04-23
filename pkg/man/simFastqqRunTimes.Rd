\name{simFastqqRunTimes}
\alias{simFastqqRunTimes}
\title{simFastqqRunTimes: For given values of k and nSeq the function creates fastq files with simulated data, collects k-mer data with the fastqq function and reports the run times for the data collection.}
\description{
For each combination of the parameters k and nSeq, the function writes one fastq file and collects the data. The fastq files are equally structured: Each read contains 17 randomly selected DNA k-mers of size 6. Therefore the read-length is always 102.}
\usage{
simFastqqRunTimes(k,nSeq,filedir=".")
}
\arguments{
  \item{k}{\code{numeric}. Length of k-mer. Default value is 2:15.}
  \item{nSeq}{\code{numeric}. Number of simulated reads in fastq-file. Default value is (100, 1000,...,10000000).}
  \item{filedir}{\code{character}. The output can be placed in a separate directory. When not existant, the function tries to create 'filedir'.}
}
\details{The fastq files contain the parameter settings inside their filename. The files are created with 'writeSimFastq'.}
\value{\code{data.frame}. The data frame has four columns: id,k,nSeq and runtime.}
\author{Wolfgang Kaisers}
\examples{
simFastqqRunTimes(k=5,nSeq=10)
}
\keyword{simFastqqRunTimes}