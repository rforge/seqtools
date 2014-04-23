\name{kMerIndex}
\alias{kMerIndex}
\title{
kMerIndex function: Counts for a given DNA k-mer subset position wise from fastq files.
}
\description{
Reads (compressed) fastq files and counts for given DNA k-mer subset for each position in sequence. The k-mer subset is given by a vector of k-mer indices. k-mer indices can be obtained from DNA k-mers with the function \code{kMerIndex}.}
\usage{
kMerIndex(kMers,k=nchar(kMers)[1])
}
\arguments{
  \item{kMers}{\code{character}. Vector of equal sized character strings. The number of characters in each string must be =k (i.e. \code{all(nchar(kMers)==k)})}
  \item{k}{Integer. Length of k-mer.}
}
\details{Maximal allowed value for k is 12.}
\value{\code{integer}.}
\author{
Wolfgang Kaisers
}
\examples{
kMerIndex(c("AACC","ATAA"))
}
\keyword{kMerIndex}