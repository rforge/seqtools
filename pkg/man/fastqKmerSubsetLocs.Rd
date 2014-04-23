\name{fastqKmerSubsetLocs}
\alias{fastqKmerSubsetLocs}
\title{
fastqKmerSubsetLocs function: Counts for a given DNA k-mer subset position wise from fastq files.
}
\description{
Reads (compressed) fastq files and counts for given DNA k-mer subset for each position in sequence. The k-mer subset is given by a vector of k-mer indices. k-mer indices can be obtained from DNA k-mers with the function \code{getKmerIndex}.}
\usage{
fastqKmerSubsetLocs(filenames,k=4,kIndex)
}
\arguments{
  \item{filenames}{Character. Vector of fastqKmerSubsetLocs file names. Files can be gz compressed.}
  \item{k}{Integer. Length of counted DNA k-mers.}
  \item{kIndex}{Integer. Numeric values which represent indices of DNA-k mers.}
}
\details{Maximal allowed value for k is 16. maxSeqLen must be >= k.}
\value{\code{list}. Contains for each given file a matrix with 4^k rows and (maxSeqLen-k+1) columns. The matrix contains for each k-mer and k-mer-start position the counted values.}
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
setwd(basedir)
k<-4
kMers<-c("AAAA","AACC","AAGG")
kIdx<-kMerIndex(kMers)
res<-fastqKmerSubsetLocs("test_l6.fq",k,kIdx)
}
\keyword{fastqKmerSubsetLocs}
\keyword{kmer}
