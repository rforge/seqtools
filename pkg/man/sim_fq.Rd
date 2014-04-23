\name{sim_fq}
\alias{sim_fq}
\title{sim_fq: Performs a simulation experiment where }
\description{
Writes compressed fastq files where sequence sections contain concatenated k-mers which are uniformly distributed in the range of k-mers for given k. A fraction of the reads can be contaminated with one or more deterministic k-mers.}
\usage{
sim_fq(nRep=2,nContamVec=c(100,1000),grSize=20,nSeq=1e4,k=6,kIndex=1365,pos=20)
}
\arguments{
  \item{nRep}{\code{numeric}. Number of replicates for each combination of each nContamVec value}
  \item{nContamVec}{\code{numeric}. Vector with nContam (absolute number of contaminated reads) values.}
  \item{grSize}{\code{numeric}. Number fastq files in control and contamination group.}
  \item{nSeq}{\code{numeric}. Number of reads per fastq file.}
  \item{k}{\code{numeric}. k value used in fastqq function.}
  \item{kIndex}{\code{numeric}. k-mer index of inserted k-mer(s). The k-mer index can be retreaved for a given k-mer with 'getKmerIndex'. Default value is 1365 (="CCCCCC").}
  \item{pos}{\code{numeric}. Determines at which position in sequence the k-mer is inserted. 1-based (1=first position).}
}
\details{The function is intended to be used as explorative tool (not for routine quality assessment). There are some files written and there will be a lot of output on the terminal. It is therefore recommended to switch to a separate working directory and to run this function on a separate terminal.}
\value{None.}
\author{Wolfgang Kaisers}
\examples{
kMerIndex("CCCCCC")
res<-sim_fq(nRep=2,nContamVec=c(10,100),grSize=4,nSeq=1e2)
}
\keyword{sim_fq}