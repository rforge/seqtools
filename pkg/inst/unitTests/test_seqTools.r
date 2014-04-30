

test_countDnaKmers<-function()
{
  ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
  res<-countDnaKmers("ACGT", k=1, start=3:1, width=1)
  checkEquals(res[,1],c(0,0,1,0))
  checkEquals(res[,2],c(0,1,0,0))
  checkEquals(res[,3],c(1,0,0,0))
  
  ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
  res<-countDnaKmers("ACGT", k=1, start=3, width=1)
  
  checkEquals(res,c(0,0,1,0))

  res<-countDnaKmers("ATTNAC", k=2, start=1:3, width=1)     
  checkEquals(res[,1],c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0))
  checkEquals(res[,2],c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1))
  checkEquals(res[,3],c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  
  
}

test_kMerIndex<-function()
{
  checkEquals(kMerIndex(c("AACC","ATAA")),c(6,49))
}


test_countSpliceKmers<-function()
{
  dna<-"atcgGTccAGatcg"
  mt<-countSpliceKmers(dna,seqid=1,lEnd=4,rStart=11,width=2,strand=1,k=3)
  
  checkEquals(sum(mt),2)
  checkTrue(all(mt[kMerIndex(c("atc","tcg"))]==1),"Check: countSpliceKmers")  
}