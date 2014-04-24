\name{revCountDnaKmers}
\alias{revCountDnaKmers}
\title{revCountDnaKmers: Counting K-mers in DNA sequence.}
\description{
Counts K-mers for reverse complement of given DNA sequence. The k-mers are searched in a set of windows, which are defined by \code{start} and \code{width} parameter. Within the defined windows, the function counts k-mers on the reverse complement of the given DNA sequence.}
\usage{
revCountDnaKmers(dna,k,start,width)
}
\arguments{
  \item{dna}{\code{character}. Vector of DNA sequences (must not contain other characters as "ATCGN". Capitalization does not matter. If any "N"-entry is found during analysis for a DNA k-mer, the k-mer is omitted from counting.)}
  \item{k}{\code{numeric}. Number of nucleotides in tabled DNA motifs.}
  \item{start}{\code{numeric}. Vector of (1-based) start positions for reading frames. Reading frame is counted to the left side of the DNA string.}
  \item{width}{\code{numeric}. Vector of frame width values}
}
\details{

Each colum containts the motif-count values for one frame. Each row represents one DNA motif. The DNA sequence of the DNA motif is given as row.name.}
\value{matrix}
\author{Wolfgang Kaisers}
\examples{
rseq<-"TATTAT"
revCountDnaKmers(rseq,2,6:4,2)
}
\keyword{revCountDnaKmers}