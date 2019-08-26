calcTLindex<-function(biomassByGroup, TLbyGroup){
  ##biomassByGroup is an array of dimensions ntracers * ntime. TLbyGroup must match order of tracers in biomassByGroup
  A<-rep(0, dim(biomassByGroup)[2]); B<-A #will use A for numerator, B for denominator
  for(t in 1:dim(biomassByGroup)[1]){
    thisTL<-TLbyGroup[t]
    if(thisTL>0){
      thisA<-biomassByGroup[t,]*thisTL; thisB<-biomassByGroup[t,]
      A<-A+thisA; B<-B+thisB
    }
  }
  TLindex<-A/B
  return(TLindex)
}