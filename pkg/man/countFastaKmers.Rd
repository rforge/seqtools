\name{countFastaKmers}
\alias{countFastaKmers}
\alias{writeFai}
\title{
countFastaKmers function: Counts DNA k-mers from (compressed) fasta files.
}
\description{
Reads (compressed) fasta files and counts for DNA k-mers in sequence.
}
\usage{
countFastaKmers(filenames,k=4)
}
\arguments{
  \item{filenames}{\code{character}: Vector of fasta file names. Files can be gz compressed.}
  \item{k}{Length of counted DNA k-mers.}
}
\details{Maximal allowed value for k is 16.}
\value{\code{matrix}.}
\references{(Fasta sequence definition)
}
\author{
Wolfgang Kaisers
}
\note{The static size of the retured k-mer array is 4^k.}
\examples{
basedir<-system.file("extdata",package="seqTools")
filename<-file.path(basedir,"small.fa")
#writeFai(filename,"small.fa.fai")
res<-countFastaKmers(filename,k=2)
}
\keyword{countFastaKmers}
\keyword{kmer}
