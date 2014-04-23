\name{collectDur}
%
\alias{collectTime}
\alias{collectTime-methods}
\alias{collectTime,Fastqq-method}
%
\alias{collectDur}
\alias{collectDur-methods}
\alias{collectDur,Fastqq-method}
%
\title{
collectDur: Returning elapsed time (in seconds) for collection of data from fastq files.
}
\description{Objects of class \code{Fastqq} are created by reading data from fastq-files using the function \code{fastqq}. The \code{fastqq} function calls \code{Sys.time()} before and after execution of the core collecting routine. \code{collectDur} returns the number of seconds between these two times (as numeric value). \code{collectTime} returns the two timestamps inside a \code{list}.}
\usage{
collectDur(object)
collectTime(object)
}
\arguments{
  \item{object}{\code{Fastqq}. Object from which collection duration (or times) is returned.}
}
\details{Maximal allowed value for k is 12.}
\value{\code{collectTime} returns \code{numeric}. \code{collectTime} returns \code{list}.}
\author{
Wolfgang Kaisers
}
\examples{
basedir<-system.file("extdata",package="seqTools")
setwd(basedir)
fq<-fastqq(c("g4_l101_n100.fq.gz","g5_l101_n100.fq.gz"),k=4,probeLabel=c("g4","g5"))
collectTime(fq)
collectDur(fq)
}
\keyword{collectDur}
\seealso{fastqq}