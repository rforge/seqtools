\name{trimFastq}
\alias{trimFastq}
\title{trimFastq: Performs removing sequences, fixed and quality base trimming and nucleotide masking on fastq files.}
\description{
Writes compressed fastq files where sequence sections contain concatenated k-mers which are uniformly distributed in the range of k-mers for given k. A fraction of the reads can be contaminated with one or more deterministic k-mers.}
\usage{
trimFastq(infile,outfile="keep.fq.gz",discard="disc.fq.gz",qualDiscard=0,qualMask=0,
                    fixTrimLeft=0,fixTrimRight=0,qualTrimLeft=0,qualTrimRight=0,
                    qualMaskValue=78,minSeqLen=0)
}
\arguments{
  \item{infile}{\code{character}. Input fastq file.}
  \item{outfile}{\code{character}. Output fastq file.}
  \item{discard}{\code{character}. Output file in which discarded reads are written.}
  \item{qualDiscard}{\code{numeric}. All reads which contain one or more phred scores <qualDiscard will be discarded (i.e. output to discard).}
  \item{qualMask}{\code{numeric}. All nucleotides for which phred score<qualMask will be overwritten with \code{qualMaskValue}.}
  \item{fixTrimLeft}{\code{numeric}. Prefix of this size will be trimmed.}
  \item{fixTrimRight}{\code{numeric}. Suffix of this size will be trimmed.}
  \item{qualMaskValue}{\code{numeric}. ASCII replace value for masked nucleotides}
  \item{qualTrimLeft}{\code{numeric}. Prefix where all phred scores are < qualTrimLeft will be trimmed.}
  \item{qualTrimRight}{\code{numeric}. Suffix where all phred scores are < qualTrimRight will be trimmed.}
  \item{minSeqLen}{\code{numeric}. All reads where sequence length after (fixed and quality based) trimming is <minSeqLen will be descarded (i.e. output to discard).}
}
\details{After trim operations, the function checks for remaining read length. When read length is smaller than minSeqLen, the read will be discarded (i.e. the read will be written to the discard file).}
\value{Numeric. A vector of length 2 which contains the number of reads which are written to output and to discard}
\author{Wolfgang Kaisers}
\examples{
basedir<-system.file("extdata",package="seqTools")
setwd(basedir)
trimFastq("sim.fq.gz",qualDiscard=10,qualMask=15,fixTrimLeft=2,
    fixTrimRight=2,qualTrimLeft=28,qualTrimRight=30,minSeqLen=5)
}
\keyword{trimFastq}