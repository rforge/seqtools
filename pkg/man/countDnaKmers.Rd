\name{countDnaKmers}
\alias{countDnaKmers}
\title{countDnaKmers: Counting k-mers in DNA sequence.}
\description{
Counts occurrence of k-mers in given DNA sequence. The k-mers are searched in a set of windows, which are defined by \code{start} and \code{width} parameter.}
\usage{
countDnaKmers(dna,k,start,width)
}
\arguments{
  \item{dna}{\code{character}. Single DNA sequence (vector of length 1). \code{dna} must not contain other characters as "ATCG". Capitalization does not matter.}
  \item{k}{\code{numeric}. Number of nucleotides in tabled DNA motifs.}
  \item{start}{\code{numeric}. Vector of (1-based) start positions for reading frames.}
  \item{width}{\code{numeric}. Vector of frame width values. When \code{width} has length 1, the value is recycled.}
}
\details{The last start positions for k-mer search are \code{start+width-1}. The last allowed start position counting a k-mer is \code{nchar(dna)-k+1}. The function throws the error \code{'Search region exceeds string end'} when a value \code{start + width + k > nchar(dna) + 2} occurs.
The function returns a matrix. Each colum contains the motif-count values for one frame. The column names are the values in the start vector. Each row represents one DNA motif. The DNA sequence of the DNA motif is given as row.name. }
\value{\code{matrix}.}
\author{Wolfgang Kaisers}
\examples{
seq<-"ATAAATA"
countDnaKmers(seq,2,1:3,3)
}
\seealso{revCountDnaKmers}
\keyword{countDnaKmers}