getICNumbers<-function(name,cohort,ThisIC.nc){
  #assumes has availabilt init which is the biol initial conditions file lines
  thisVar<-paste(name,cohort,"_Nums",sep="")
  xx<-ncvar_get(ThisIC.nc,thisVar)
  if(length(dim(xx))==3){xx<-xx[,,1]}
  totalNums<-sum(xx,na.rm=TRUE)
  return(totalNums)
}