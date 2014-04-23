

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
  checkEquals(kMerIndex(c("AACC","ATAA")),c(5,48))
}
