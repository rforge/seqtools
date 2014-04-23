\name{ascii2char}
\alias{ascii2char}
\alias{char2ascii}
\title{ascii2char: Converting ASCII encoded values to character values.}
\description{
The function calculates phred values for a given range. As default, a data.frame with all phred values is returned. \code{char2ascii} does the reverse transformation.}
\usage{
ascii2char(x,multiple = FALSE)
char2ascii(c)
}
\arguments{
  \item{x}{\code{numeric}. Vector with ASCII values. All values must be in 1:255. Other values produce an error.}
  \item{multiple}{\code{logical}. For \code{'FALSE'} (the default), all characters are combined into one single string (i.e. a character vector of length 1). For \code{'TRUE'}, single characters are combined into a vector.}
  \item{c}{\code{character}. Vector of length 1 (Longer vectors will generate Warnings).}
}
\details{The functions are only wrappers for convenience. \code{char2ascii} is defined as \code{strtoi(charToRaw(c), base = 16L)}. \code{ascii2char} is defined as \code{rawToChar(as.raw(x), multiple)}.}
\value{\code{ascii2char} returns a \code{character}. \code{char2ascii} returns \code{integer}.}
\author{Wolfgang Kaisers}
\examples{
ascii2char(97:101,multiple=FALSE)
ascii2char(97:101,multiple=TRUE)
char2ascii("abcde")
ascii2char(char2ascii("abcde"))
}
\seealso{getPhredTable}
\keyword{ascii2char}