\name{mergeFastqq}
%
\alias{mergeFastqq}
\alias{mergeFastqq-methods}
\alias{mergeFastqq,Fastqq-method}
%
\title{mergeFastqq: Merges two Fastqq object into one.}
\description{The contents of two given Fastqq objects are merged together into one resulting Fastqq object.}
\usage{
mergeFastqq(lhs,rhs)
}
\arguments{
  \item{lhs}{\code{Fastqq}.}
  \item{rhs}{\code{Fastqq}.}
}
\details{The function returns a data frame. Each colum containts the motif-count values for one frame. Each row represents one DNA motif. The DNA sequence of the DNA motif is given as row.name.}
\value{S4 Object of class 'Fastqq'.}
\author{Wolfgang Kaisers}
\examples{
basedir<-system.file("extdata",package="seqTools")
setwd(basedir)
lhs<-fastqq("g4_l101_n100.fq.gz",k=4,"g4")
rhs<-fastqq("g5_l101_n100.fq.gz",k=4,"g5")
fq<-mergeFastqq(lhs,rhs)
}
\keyword{mergeFastqq}