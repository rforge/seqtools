
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
##  Project   :   seqTools                                                   ##
##  Created   :   26.August.2013                                             ##
##  Author    :   W. Kaisers                                                 ##
##  File      :   kMer.r                                                     ##
##  Content   :   Functionality for counting DNA k-mers                      ##
##                (independent of fastq or fasta files)                      ##
##                countKmers, countDnaKmers, revCountDnaKmers,               ##
##                countGenomeKmers, countSpliceKmers                         ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Counts DNA k-mers on specified regions inside single (character) sequence
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
countDnaKmers <- function(dna, k, start=1, width=nchar(dna) - k + 1)
{
  if(!is.character(dna))
    stop("'dna' must be character.")
  if(length(dna) != 1)
    stop("'dna' must have length 1.")
  if(is.numeric(start))
    start <- as.integer(start)
  if(is.numeric(width))
    width <- as.integer(width)
  if(length(width) == 1)
    width <- rep(width, length(start))
  if(is.numeric(k))
    k <- as.integer(k)
  if(length(k)!=1)
    stop("'k' must have length 1.")
  if(k < 1)
    stop("'k' must be positive.")
  if(k > max_k)
    stop("'k' must not exceed", max_k, ".")
  
  nc<-nchar(dna)
  if(k>nc)
    stop("'k' must be <= nchar(dna).")
  if(any(start+width+k>nc+2))
    stop("Search region exceeds string end.")
  
  
  # Counts N's
  # ToDo: Return value
  nn <- integer(length(start))
  return(.Call("count_dna_Kmers", dna, start, width, k, nn, 
               PACKAGE = "seqTools"))
}

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Counts DNA k-mers on specified regions inside single (character) sequence 
## in reverse direction
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
revCountDnaKmers <- function(dna, k, start, width)
{
  if(!is.character(dna))
    stop("'dna' must be character.")
  if(length(dna) != 1)
    stop("'dna' must have length 1.")
  if(is.numeric(start))
    start <- as.integer(start)
  if(is.numeric(width))
    width <- as.integer(width)
  if(length(width) == 1)
    width <- rep(width, length(start))
  if(is.numeric(k))
    k <- as.integer(k)  
  if(any(width + k > start))
    stop("'width' must be <=  'start' - 'k'.")
  
  # Counts N's
  # ToDo: Return value
  nn <- integer(length(start))
  return(.Call("rev_count_dna_Kmers", dna, start, width, k, nn,
               PACKAGE = "seqTools"))
}

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Counts DNA k-mers on specified regions inside multiple (character) sequences
##  in possibly reversed direction (depending on wStrand)
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
countGenomeKmers <- function(dna, wSeqid, wStart, wWidth, wStrand, k)
{
  if(!is.character(dna))
    stop("'dna' must be character.")
  if(!is.numeric(wSeqid))
    wSeqid <- as.integer(wSeqid)
  rg <- range(wSeqid)
  if(rg[1] < 0)
    stop("Negative seqid's are not allowed.")
  if(rg[2] > length(dna))
    stop("Out of range seqid's.")
  
  if(!is.numeric(wStart))
    stop("'wStart' must be numeric.")
  wStart <- as.integer(wStart)
  if(!is.numeric(wWidth))
    stop("'wWidth' must be numeric")
  wWidth <- as.integer(wWidth)
  
  if(is.factor(wStrand))
    wStrand <- as.integer(wStrand)
  else
  {
    if(!is.numeric("wStrand"))
      wStrand <- as.integer(wStrand)    
  }
  
  nStart <- length(wStart)
  if( (length(wSeqid) != nStart) | (length(wWidth) != nStart) | 
        (length(wStrand) != nStart) )
    stop("'wSeqid', 'wStart', 'wWidth' and 'wStrand' must have same length.")
  if(length(k) != 1)
    stop("'k' must be a single value.")
  if(k>max_k)
    stop("'k' must not exceed", max_k, ".")
  
  # Counts N's
  # ToDo: Return value
  nn <- integer(length(wStart))
  
  return(.Call("count_genome_Kmers", dna, wSeqid, wStart, wWidth, 
               wStrand, k, nn, PACKAGE = "seqTools")) 
}


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Counts DNA k-mers on each border of a splice-site defined by wLend and 
## wRstart in range of size wWidth
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
countSpliceKmers <- function(dna, wSeqid, wLend, wRstart, wWidth, wStrand, k)
{
  if(!is.character(dna))
    stop("'dna' must be character.")
  if(!is.numeric(wSeqid))
    stop("'wSeqid' must be numeric.")
  wSeqid <- as.integer(wSeqid)
  
  rg <- range(wSeqid)
  if(rg[1] < 0)
    stop("Negative seqid's are not allowed.")
  if(rg[2] > length(dna))
    stop("Out of range seqid's.")
  
  if(!is.numeric(wLend))
    stop("'wLend' must be numeric.")
  wLend <- as.integer(wLend)
  if(!is.numeric(wRstart))
    stop("'wRstart' must be numeric.")
  wRstart <- as.integer(wRstart)
  
  if(!is.numeric(wWidth))
    stop("'wWidth' must be numeric.")
  wWidth <- as.integer(wWidth)
  
  if(is.factor(wStrand))
    wStrand <- as.integer(wStrand)
  else
  {
    if(!is.numeric("wStrand"))
      wStrand <- as.integer(wStrand)    
  }
  
  nStart <- length(wLend)
  if(length(wSeqid) != nStart | length(wRstart) != nStart | 
       length(wWidth) != nStart | length(wStrand) != nStart)
    stop(
"'wSeqid', 'wLend', 'wRstart', 'wWidth' and 'wStrand' must have equal length.")

  if(!is.numeric(k))
    stop("'k' must be numeric.")
  k <- as.integer(k)
  if(length(k) != 1)
    stop("'k' must be a single value.")
  if(k > max_k)
    stop("'k' must not exceed", max_k, ".")
  
  # Counts N's
  # ToDo: Return value
  nn <- integer(length = nStart)
  
  return(.Call("count_splice_Kmers", dna, wSeqid, wLend, wRstart, wWidth, 
               wStrand, k, nn, PACKAGE = "seqTools")) 
}

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## END OF FILE
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##