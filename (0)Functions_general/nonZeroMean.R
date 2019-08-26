nonZeroMean<-function(x){
  x[x==0]<-NA
  thisMean<-mean(x,na.rm=TRUE)
  return(thisMean)
}