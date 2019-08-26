calcR<-function(SSB,a,b,SSBnumerator=NULL){
  if(is.null(SSBnumerator)){SSBnumerator<-SSB}
  R<-(a*SSBnumerator)/(b+SSB)
  return(R)
}