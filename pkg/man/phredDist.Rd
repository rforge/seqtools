\name{phredDist}
%
\alias{phredDist}
\alias{phredDist-methods}
\alias{phredDist,Fastqq-method}
%
\alias{plotPhredDist}
\alias{plotPhredDist-methods}
\alias{plotPhredDist,Fastqq-method}
%
\title{
phredDist: Global relative content of Phred values in Fastqq objects (or subsets).
}
\description{
The \code{phredDist} function returns a named vector with relative Phred content from the whole \code{Fastqq} object or a subset which is denoted by a index \code{i}. The \code{plotPhredDist} function produces a plot of the \code{phredDist} values.}
\usage{
phredDist(object,i)
plotPhredDist(object, i, maxp=45, col, ...)
}
\arguments{
  \item{object}{\code{Fastqq}: Object which contains collected values from \code{nFiles} fastq files.}
  \item{i}{\code{integer}(optional): Index of fastq file(s) from which Phred values are counted. When value is missing, Phred counts for all contained data is returned.}
  \item{maxp}{\code{numeric}(optional): Value of maximal plotted phred value (right limit of x-axis).}
  \item{col}{Colour encoding for plotted lines.}
  \item{...}{Additional values passed to plot function.}
}
\details{i must be a numerical vector with values in 1,...,nFiles. The \code{plotPhredDist} function is also prepared for additional arguments: The maxp value denotes the maximal Phred value until which the Phred values are plottted (evtl. shrinks the x-Axis). The standard line color is \code{topo.colors(10)[3]}. Additional arguments (e.g. main="") can be passed to the plot function.}
\value{\code{phredDist} returns \code{numeric}. \code{plotPhredDist} returns nothing.}
\references{
Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM
The sanger fastq file format for sequences with quality scores and the Solexa/Illumina fastq variants.
Nucleic Acids Research 2010 Vol.38 No.6 1767-1771
}
\author{
Wolfgang Kaisers
}
\examples{
basedir<-system.file("extdata",package="seqTools")
setwd(basedir)
fq<-fastqq(c("g4_l101_n100.fq.gz","g5_l101_n100.fq.gz"),k=4,probeLabel=c("g4","g5"))
#
phredDist(fq)
plotPhredDist(fq,main="g4 and g5")
#
}
\keyword{phredDist}
\keyword{plotPhredDist}
