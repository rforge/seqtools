\name{countGenomeKmers}
\alias{countGenomeKmers}
\title{countGenomeKmers: Counting K-mers in DNA sequences.}
\description{
Counts K-mers in DNA sequences}
\usage{
countGenomeKmers(dna,wSeqid,wStart,wWidth,wStrand,k)
}
\arguments{
  \item{dna}{\code{character}. Vector of DNA sequences (must not contain other characters as "ATCG". Capitalization does not matter.)}
  \item{wSeqid}{\code{numeric}. Vector of (1-based) values depicting one of the given sequences.}
  \item{wStart}{\code{numeric}. Vector of (1-based) start positions for reading windows.}
  \item{wWidth}{\code{numeric}. Vector of window width values.}
  \item{wStrand}{\code{factor or numeric}. First factor level (or numeric: 1) value will be interpreted as (+)-strand. For any other values, the reversed complement sequence will be counted (in left direction from start value).}
  \item{k}{\code{numeric}. Number of nucleotides in tabled DNA motifs. Only a single value is allowed (length(n)=1!)}  
}
\details{The function returns a matrix. Each colum containts the motif-count values for one frame. Each row represents one DNA motif. The DNA sequence of the DNA motif is given as row.name.}
\value{\code{matrix}.}
\author{Wolfgang Kaisers}
\examples{
sq<-"TTTTTCCCCGGGGAAAA"
wSeqid <-as.integer(c(1, 1))
wLstart<-as.integer(c(6,14))
wWidth <-as.integer(c(4, 4))
wStrand<-as.integer(c(1, 0))
wKmer  <-2L
countGenomeKmers(sq,wSeqid,wLstart,wWidth,wStrand,wKmer)
}
\keyword{countGenomeKmers}