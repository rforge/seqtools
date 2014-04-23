\name{countSpliceKmers}
\alias{countSpliceKmers}
\title{countSpliceKmers: Counting K-mers in DNA sequences.}
\description{
Counts K-mers in DNA sequences}
\usage{
countSpliceKmers(dna,wSeqid,wLend,wRstart,wWidth,wStrand,k)
}
\arguments{
  \item{dna}{\code{character}. Vector of DNA sequences (must not contain other characters as "ATCG". Capitalization does not matter.)}
  \item{wSeqid}{\code{numeric}. Vector of (1-based) values depicting one of the given sequences.}
  \item{wLend}{\code{numeric}. Vector of (1-based) left-end positions. Will be used as rightmost window positions (over which (n-1) positions overhang will be counted as part of frame).}
  \item{wRstart}{\code{numeric}. Vector of (1-based) right-start positions. Will be used as leftmost window positions (over which(n-1) positions overhang will be counted as part of frame).}
  \item{wWidth}{\code{numeric}. Vector of window width values.}  
  \item{wStrand}{\code{factor or numeric}. First factor level (or numeric: 1) value will be interpreted as (+) strand. For any other values, the reversed complement sequence will be counted (in left direction from start value). For (+) strand, the wLend value will be used as starting position. For (-) strand, the wRstart position will be used as starting positions.}
  \item{k}{\code{numeric}. Number of nucleotides in tabled DNA motifs. Only a single value is allowed (length(n)=1!)}  
}
\details{The function returns a matrix. Each colum containts the motif-count values for one frame. Each row represents one DNA motif. The DNA sequence of the DNA motif is given as row.name.}
\value{\code{matrix}.}
\author{Wolfgang Kaisers}
\examples{
sq1<-"TTTTTCCCCGGGGAAAA"
sq2<-"TTTTTTTCCCCGGGGAAAA"
sq<-c(sq1,sq2)
wSeqid <-c( 1, 1, 2, 2)
wLend  <-c( 9, 9,11,11)
wRstart<-c(14,14,16,16)
wWidth <-c( 4, 4, 4, 4)
wStrand<-c( 1, 0, 1, 0)
wKmer  <-2
countSpliceKmers(sq,wSeqid,wLend,wRstart,wWidth,wStrand,wKmer)
}
\keyword{countSpliceKmers}