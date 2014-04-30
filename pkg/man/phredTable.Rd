\name{phredTable}
\alias{phredTable}
\title{phredTable: Returns a data.frame with phred encodings.}
\description{
The function calculates phred values for a given range. As default, a data.frame with all phred values is returned.}
\usage{
phredTable(phred)
}
\arguments{
  \item{phred}{\code{numeric}. Vector with phred values. All values must be in 0:93}
}
\value{\code{data.frame}. The data.frame has three columns: "ascii","phred" and "char"}
\references{
Ewing B, Green P
Base-calling of automated sequencer traces using phred. II. Error probabilities.
Genome Research 1998 Vol. 8 No. 3 186-194
}
\author{Wolfgang Kaisers}
\examples{
phredTable()
}
\keyword{phredTable}
\seealso{char2ascii}